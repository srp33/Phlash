"""Contains the functions for the CDS page.

Returns CDS data.
Updates CDS data.
Parses BLAST results.

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

response_object = {}

# ------------------------------ MAIN FUNCTIONS ------------------------------
def annotate_cds(request, cds_id):
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
    genemark_cds = GeneMark.query.filter_by(stop=cds.stop).first()
    if cds:
        cds.start = put_data.get('start')
        cds.stop = put_data.get('stop')
        cds.function = put_data.get('function')
        response_object['message'] = 'CDS updated!'
    else:
        response_object['message'] = 'CDS did not update.'
    if genemark_cds:
        orf_length_dnamaster = cds.stop - cds.start
        orf_length_genemark = genemark_cds.stop - genemark_cds.stop
        if orf_length_dnamaster < orf_length_genemark:
            cds.status = "Fail"
        else:
            cds.status = "Pass"
    else:
        cds.status = "Undetermined"
        
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
                                'status': cds.status}
    
    starts = [int(start) for start in cds.start_options.split(",")]
    stops = [int(stop) for stop in cds.stop_options.split(",")]
    response_object['start_options'] = starts
    response_object['stop_options'] = stops
    print(stops)
    print(starts)
    blast_files = helper.get_file_path("blast", UPLOAD_FOLDER)
    E_VALUE_THRESH = 1e-7
    blast_results = parse_blast_results(blast_files, cds.id, E_VALUE_THRESH)
    response_object['blast'] = blast_results

    genemark_gdata_file = helper.get_file_path("gdata", UPLOAD_FOLDER)
    gdata_df = pd.read_csv(genemark_gdata_file, sep='\t', skiprows=16)
    gdata_df.columns = ['Base', '1', '2', '3', '4', '5', '6']
    gdata_df = gdata_df[gdata_df.Base.isin(
        range(min(starts) - 100, cds.stop + 100))]
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
            print(cds.id)
            break
        elif cds.id == cds_id:
            reached_CDS = True

    return response_object

# ---------- BLAST HELPER FUNCTIONS ----------
def parse_blast_results(blast_files, cds_id, e_value_thresh):
    """Parses through the blast results and returns a dictionary containing data.

    Removes instances of CREATE_VIEW created by BLAST.
    Finds the associated blast results for the CDS.

    Args:
        blast_files:
            The files containing the blast results.
        cds_id:
            The ID of the CDS to be found.
        e_value_thresh:
            The blast similarity threshold.

    Returns:
        A dictionary containing all of the blast data for the CDS.
    """
    blast_results = {}
    for blast_file in blast_files:
        newLines = []
        with open(blast_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line != "CREATE_VIEW\n":
                    newLines.append(line)
        with open(blast_file, 'w') as f:
            f.writelines(newLines)    
            
        with open(blast_file) as f:
            blasts = json.load(f)["BlastOutput2"]
            for blast in blasts:
                search = blast["report"]["results"]["search"]
                title = re.search(
                    "([A-Z]+_\d+_*\d*), (\d+-\d+)", search["query_title"])
                if title:
                    curr_id = title.group(1)
                    curr_start = title.group(2)
                    cds_id_ = cds_id + "_"
                    if cds_id == curr_id or cds_id_ in curr_id:
                        blast_results[curr_start] = []
                        hits = search["hits"]
                        for hit in hits:
                            hsps = hit["hsps"][0]
                            if hsps["evalue"] <= e_value_thresh:
                                alignment = {}
                                description = hit["description"][0]
                                alignment['accession'] = description["accession"]
                                alignment["title"] = description["title"]
                                alignment["evalue"] = '{:0.2e}'.format(
                                    hsps["evalue"])
                                alignment["query_from"] = hsps["query_from"]
                                alignment["query_to"] = hsps["query_to"]
                                alignment["hit_from"] = hsps["hit_from"]
                                alignment["hit_to"] = hsps["hit_to"]
                                alignment["percent_identity"] = round(
                                    hsps["identity"] / hsps["align_len"] * 100, 2)
                                blast_results[curr_start].append(alignment)

    return blast_results