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
def annotate_cds(phage_id, request, cds_id, UPLOAD_FOLDER):
    """Updates a CDS given the data and the ID.

    Updates the left, right, and function of a given CDS.
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
    cds = Annotations.query.filter_by(phage_id=phage_id).filter_by(id=cds_id).first()
    if cds:
        if put_data.get('strand'):
            cds.strand = put_data.get('strand')
            cds.left = put_data.get('left')
            cds.right = put_data.get('right')
            cds.function = put_data.get('function')
            cds.notes = put_data.get('notes')
        else:
            cds.notes = put_data.get('notes')
        response_object['message'] = 'CDS updated!'
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
        frame, status = helper.get_frame_and_status(cds.left, cds.right, cds.strand, coding_potential)
        cds.status = status
        cds.frame = frame
    db.session.commit()

    return response_object

def get_cds_data(phage_id, UPLOAD_FOLDER, cds_id):
    """Queries and returns all of the data for a CDS given the ID.

    Gets left, right, function, and status for the CDS.
    Gets the right of the previous CDS and the left of the next CDS.
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
    cds = Annotations.query.filter_by(phage_id=phage_id).filter_by(id=cds_id).first()
    prev_cds = Annotations.query.filter_by(phage_id=phage_id).filter_by(id=prev_id).first()
    next_cds = Annotations.query.filter_by(phage_id=phage_id).filter_by(id=next_id).first()
    if prev_cds is not None:
        response_object['prevCDS'] = prev_id
        response_object['prev_right'] = prev_cds.right
    else:
        response_object['prevCDS'] = 'undefined'
        response_object['prev_right'] = 0
    if next_cds is not None:
        response_object['next_left'] = next_cds.left
    else:
        response_object['next_left'] = cds.right
    response_object['cds'] = {'id': cds.id,
                                'left': cds.left,
                                'right': cds.right,
                                'strand': cds.strand,
                                'function': cds.function,
                                'status': cds.status,
                                'frame': cds.frame,
                                'notes': cds.notes}
    
    left_positions, right_positions = get_blasts(phage_id, cds.left)

    genemark_gdata_file = helper.get_file_path("gdata", UPLOAD_FOLDER)
    gdata_df = pd.read_csv(genemark_gdata_file, sep='\t', skiprows=16)
    gdata_df.columns = ['Base', '1', '2', '3', '4', '5', '6']
    try:
        gdata_df = gdata_df[gdata_df.Base.isin(
            range(min(left_positions) - 100, max(right_positions) + 100))]
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
    for cds in db.session.query(Annotations).filter_by(phage_id=phage_id).order_by(Annotations.left):
        if reached_CDS and cds.function != "@DELETED" and cds.status != "tRNA":
            response_object['nextCDS'] = cds.id
            break
        elif cds.id == cds_id:
            reached_CDS = True

    reached_CDS = False
    response_object['prevCDS'] = 'undefined'
    for cds in db.session.query(Annotations).filter_by(phage_id=phage_id).order_by(Annotations.left.desc()):
        if reached_CDS and cds.function != "@DELETED" and cds.status != "tRNA":
            response_object['prevCDS'] = cds.id
            break
        elif cds.id == cds_id:
            reached_CDS = True
    response_object['glimmer'] = Gene_Calls.query.filter_by(phage_id=phage_id).filter_by(id='Glimmer').first().calls.split(',')
    response_object['genemark'] = Gene_Calls.query.filter_by(phage_id=phage_id).filter_by(id='GeneMark').first().calls.split(',')
    response_object['phanotate'] = Gene_Calls.query.filter_by(phage_id=phage_id).filter_by(id='Phanotate').first().calls.split(',')


    return response_object

# ---------- BLAST HELPER FUNCTIONS ----------
def get_blasts(phage_id, left):
    """Queries and returns all of the data for a CDS given the current left position.

    Gets alternate lefts and rights within range defined in settings
    Gets the blast results for all alternate CDS.

    Args:
        left:
            The left position of the current CDS.
    Returns:
        Lists of the alternate lefts and rights.
        
    """
    dir_blasts = {}
    comp_blasts = {}
    all_blasts = {}
    dir_lefts = []
    dir_rights = []
    comp_lefts = []
    comp_rights = []
    lefts = []
    rights = []
    setting = db.session.query(Settings).filter_by(phage_id=phage_id).first()
    minimum = left - setting.back_left_range
    maximum = left + setting.forward_left_range
    for blast in db.session.query(Blast_Results).filter_by(phage_id=phage_id).filter_by(strand='+').order_by(Blast_Results.left):
        if blast.left > minimum and blast.left < maximum:
            lefts.append(blast.left)
            rights.append(blast.right)
            dir_lefts.append(blast.left)
            dir_rights.append(blast.right)
            dir_blasts[str(blast.left) + '-' + str(blast.right) + '  ' + blast.strand] = eval(blast.results)
            all_blasts[str(blast.left) + '-' + str(blast.right) + '  ' + blast.strand] = eval(blast.results)
    for blast in db.session.query(Blast_Results).filter_by(phage_id=phage_id).filter_by(strand='-').order_by(Blast_Results.left):
        if blast.left > minimum and blast.left < maximum:
            lefts.append(blast.left)
            rights.append(blast.right)
            comp_lefts.append(blast.left)
            comp_rights.append(blast.right)
            comp_blasts[str(blast.left) + '-' + str(blast.right) + '  ' + blast.strand] = eval(blast.results)
            all_blasts[str(blast.left) + '-' + str(blast.right) + '  ' + blast.strand] = eval(blast.results)
    response_object['comp_left_options'] = comp_lefts
    response_object['comp_right_options'] = comp_rights
    response_object['dir_left_options'] = dir_lefts
    response_object['dir_right_options'] = dir_rights
    response_object['dir_blast'] = dir_blasts
    response_object['comp_blast'] = comp_blasts
    response_object['all_blast'] = all_blasts
    return lefts, rights