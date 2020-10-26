"""
Contains the methods for the CDS page.
"""
import models
from models import *
import helper
import pandas as pd
import re
import json

response_object = {}

def annotate_cds(request, cds_id):
    '''
    Updates a CDS given the data and the ID.
    '''
    put_data = request.get_json()
    cds = DNAMaster.query.filter_by(id=cds_id).first()
    if cds:
        cds.start = put_data.get('start')
        cds.function = put_data.get('function')
        cds.status = put_data.get('status')
        db.session.commit()
        response_object['message'] = 'CDS updated!'
    else:
        response_object['message'] = 'CDS did not update.'

    return response_object

def get_cds_data(UPLOAD_FOLDER, cds_id):
    '''
    Queries and returns all of the data for a CDS given the ID
    '''
    cds = DNAMaster.query.filter_by(id=cds_id).first()
    response_object['cds'] = {'id': cds.id,
                                'start': cds.start,
                                'stop': cds.stop,
                                'strand': cds.strand,
                                'function': cds.function,
                                'status': cds.status}

    print(cds)
    
    starts = [int(start) for start in cds.start_options.split(",")]
    stops = [int(stop) for stop in cds.stop_options.split(",")]
    print(starts)
    print(stops)
    response_object['start_options'] = starts
    response_object['stop_options'] = stops
    blast_files = helper.get_file_path("blast", UPLOAD_FOLDER)
    E_VALUE_THRESH = 1e-7
    print("parsing blast")
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

    dnamaster = []
    reachedCDS = False
    for cds in db.session.query(DNAMaster).order_by(DNAMaster.start):
        if reachedCDS and cds.function == "None selected":
            response_object['nextCDS'] = cds.id
            print(cds.id)
            break
        elif cds.id == cds_id:
            reachedCDS = True

    return response_object

def parse_blast_results(blast_files, cds_id, e_value_thresh):
    '''
    Parses through the blast results and returns dictionary containing data.
    '''
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
                    "([A-Z]+_\d+_*\d*), (\d+)-\d+", search["query_title"])
                if title:
                    curr_id = title.group(1)
                    curr_start = int(title.group(2))
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
