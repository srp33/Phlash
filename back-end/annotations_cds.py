"""Contains the functions for the CDS page.

Returns CDS data.
Updates CDS data.
Queries BLAST results.

Attributes:
    response_object:
        The dictionary that is returned by the main functions.
"""
import models
from models import *
import helper
import pandas as pd
import re
import json
from datetime import datetime

response_object = {}

# ------------------------------ MAIN FUNCTIONS ------------------------------
def annotate_cds(request, cds_id, UPLOAD_FOLDER):
    """Updates a CDS given the data and the ID.

    Updates the start, stop, and function of a given CDS.
    If CDS does not exist returns an error message.
    Updates the CDS status.

    Args:
        request:
            A dictionary containing the new CDS data.
        cds_id:
            The ID of the CDS to be updated.

    Returns:
        A dictionary containing a pass or fail message.
    """
    put_data = request.get_json()
    cds = DNAMaster.query.filter_by(id=cds_id).first()
    if cds:
        cds.strand = put_data.get('strand')
        cds.start = put_data.get('start')
        cds.stop = put_data.get('stop')
        cds.function = put_data.get('function')
        cds.notes = put_data.get('notes')
        print(cds.notes)
        response_object['message'] = 'CDS updated!'
        print(cds)
    else:
        response_object['message'] = 'CDS did not update.'
    coding_potential = {}
    genemark_gdata_file = helper.get_file_path("gdata", UPLOAD_FOLDER)
    gdata_df = pd.read_csv(genemark_gdata_file, sep='\t', skiprows=16)
    gdata_df.columns = ['Base', '1', '2', '3', '4', '5', '6']
    coding_potential['x_data'] = gdata_df["Base"].to_list()
    coding_potential['y_data_1'] = gdata_df["1"].to_list()
    coding_potential['y_data_2'] = gdata_df["2"].to_list()
    coding_potential['y_data_3'] = gdata_df["3"].to_list()
    coding_potential['y_data_4'] = gdata_df["4"].to_list()
    coding_potential['y_data_5'] = gdata_df["5"].to_list()
    coding_potential['y_data_6'] = gdata_df["6"].to_list()
    
    if (cds.status == "trnaDELETED" or cds.status == "tRNA"):
        cds.status = put_data.get('status')
        cds.frame = put_data.get('frame')
    else:
        frame, status = helper.get_frame_and_status(cds.start, cds.stop, cds.strand, coding_potential)
        cds.status = status
        cds.frame = frame
    db.session.commit()

    return response_object

def get_cds_data(UPLOAD_FOLDER, cds_id):
    """Queries and returns all of the data for a CDS given the ID.

    Gets start, stop, function, and status for the CDS.
    Gets the stop of the previous CDS and the start of the next CDS.
    Gets the blast results for the CDS.
    Gets the genemark coding potential data.
    Gets the next non-updated CDS ID.

    Args:
        UPLOAD_FOLDER:
            The folder containing all of the uploaded files.
        cds_id:
            The ID of the requested CDS.
    Returns:
        A dictionary containing the CDS data including the blast results and coding potential.
        
    """
    num_begins = cds_id.rfind('_') + 1
    index = float(cds_id[num_begins:])
    index = int(index)
    prev_id = cds_id[:num_begins] + str(index - 1)
    next_id = cds_id[:num_begins] + str(index + 1)
    cds = DNAMaster.query.filter_by(id=cds_id).first()
    prev_cds = DNAMaster.query.filter_by(id=prev_id).first()
    next_cds = DNAMaster.query.filter_by(id=next_id).first()
    if prev_cds is not None:
        response_object['prev_stop'] = prev_cds.stop
    else:
        response_object['prev_stop'] = 0
    if next_cds is not None:
        response_object['next_start'] = next_cds.start
    else:
        response_object['next_start'] = cds.stop
    response_object['cds'] = {'id': cds.id,
                                'start': cds.start,
                                'stop': cds.stop,
                                'strand': cds.strand,
                                'function': cds.function,
                                'status': cds.status,
                                'frame': cds.frame,
                                'notes': cds.notes}
    
    starts, stops = get_blasts(cds.start)

    genemark_gdata_file = helper.get_file_path("gdata", UPLOAD_FOLDER)
    gdata_df = pd.read_csv(genemark_gdata_file, sep='\t', skiprows=16)
    gdata_df.columns = ['Base', '1', '2', '3', '4', '5', '6']
    try:
        gdata_df = gdata_df[gdata_df.Base.isin(
            range(min(starts) - 100, max(stops) + 100))]
    except:
        response_object['message'] = "Not finished parsing"
        return response_object
    response_object['message'] = "Finished"
    response_object['x_data'] = gdata_df["Base"].to_list()
    response_object['y_data_1'] = gdata_df["1"].to_list()
    response_object['y_data_2'] = gdata_df["2"].to_list()
    response_object['y_data_3'] = gdata_df["3"].to_list()
    response_object['y_data_4'] = gdata_df["4"].to_list()
    response_object['y_data_5'] = gdata_df["5"].to_list()
    response_object['y_data_6'] = gdata_df["6"].to_list()

    reached_CDS = False
    response_object['nextCDS'] = 'undefined'
    for cds in db.session.query(DNAMaster).order_by(DNAMaster.start):
        if reached_CDS and cds.function == "None selected":
            response_object['nextCDS'] = cds.id
            break
        elif cds.id == cds_id:
            reached_CDS = True
    response_object['glimmer'] = Gene_Calls.query.filter_by(id='Glimmer').first().calls.split(',')
    response_object['genemark'] = Gene_Calls.query.filter_by(id='GeneMark').first().calls.split(',')
    response_object['phanotate'] = Gene_Calls.query.filter_by(id='Phanotate').first().calls.split(',')


    return response_object

# ---------- BLAST HELPER FUNCTIONS ----------
def get_blasts(start):
    """Queries and returns all of the data for a CDS given the current start position.

    Gets alternate starts and stops within range defined in settings
    Gets the blast results for all alternate CDS.

    Args:
        start:
            The start position of the current CDS.
    Returns:
        Lists of the alternate starts and stops.
        
    """
    dir_blasts = {}
    comp_blasts = {}
    dir_starts = []
    dir_stops = []
    comp_starts = []
    comp_stops = []
    starts = []
    stops = []
    strands = []
    setting = db.session.query(Settings).order_by(Settings.back_start_range).first()
    minimum = start - setting.back_start_range
    maximum = start + setting.forward_start_range
    for blast in db.session.query(Blast_Results).filter_by(strand='+').order_by(Blast_Results.start):
        if blast.start > minimum and blast.start < maximum:
            starts.append(blast.start)
            stops.append(blast.stop)
            dir_starts.append(blast.start)
            dir_stops.append(blast.stop)
            dir_blasts[str(blast.start) + '-' + str(blast.stop) + '  ' + blast.strand] = eval(blast.results)
    for blast in db.session.query(Blast_Results).filter_by(strand='-').order_by(Blast_Results.start):
        if blast.start > minimum and blast.start < maximum:
            starts.append(blast.start)
            stops.append(blast.stop)
            comp_starts.append(blast.start)
            comp_stops.append(blast.stop)
            comp_blasts[str(blast.start) + '-' + str(blast.stop) + '  ' + blast.strand] = eval(blast.results)
    print(dir_blasts)
    response_object['comp_start_options'] = comp_starts
    response_object['comp_stop_options'] = comp_stops
    response_object['dir_start_options'] = dir_starts
    response_object['dir_stop_options'] = dir_stops
    response_object['dir_blast'] = dir_blasts
    response_object['comp_blast'] = comp_blasts
    return starts, stops