"""Contains the functions for the Annotations page.

Returns autoannotation data.
Modifies and returns the GenBank file.
Adds a CDS.
Parses BLAST data.

Attributes:
    response_object:
        The dictionary that is returned by the main functions.
"""
import helper
from Bio import SeqIO, Seq, SeqFeature
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from models import *
import re
import json
import pandas as pd
from datetime import datetime
import os
import time

response_object = {}

# ------------------------------ MAIN FUNCTIONS ------------------------------
def get_annotations_data(phage_id):
    """Queries and returns all CDS data created by Annotations.

    Returns:
        A dictionary containing the Annotations data.
    """
    setting = db.session.query(Settings).filter_by(phage_id=phage_id).order_by(Settings.back_left_range).first()
    response_object['gap'] = setting.gap
    response_object['opposite_gap'] = setting.opposite_gap
    response_object['overlap'] = setting.overlap
    response_object['short'] = setting.short
    annotations = []
    for cds in db.session.query(Annotations).filter_by(phage_id=phage_id).order_by(Annotations.left):
        function = cds.function[:cds.function.find('#')]
        if cds.function == "None selected":
            function = "None selected"
        annotations.append({'id': cds.id,
                            'left': cds.left,
                            'right': cds.right,
                            'strand': cds.strand,
                            'function': function,
                            'status': cds.status,
                            'frame': cds.frame})
    response_object['annotations'] = annotations

    return response_object

def parse_blast(phage_id, UPLOAD_FOLDER):
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
    message = "success"
    db.session.query(Blast_Results).filter_by(phage_id=phage_id).delete()
    print(datetime.now())
    blast_files = helper.get_file_path("blast", UPLOAD_FOLDER)
    E_VALUE_THRESH = 1e-7
    counter = 0
    try:
        for blast_file in blast_files:
            print(blast_file)
            newLines = []
            with open(blast_file, 'r') as f:
                lines = f.readlines()
                skip = 0
                for i, line in enumerate(lines):
                    if skip <= 0:
                        if line.endswith("CREATE_VIEW\n"):
                            line = line[0:-12] + lines[i + 3]
                            newLines.append(line)
                            skip = 4
                        else:
                            newLines.append(line)
                    skip -= 1
            with open(blast_file, 'w') as f:
                f.writelines(newLines)    
                
            with open(blast_file) as f:
                try:
                    blasts = json.load(f)["BlastOutput2"]
                except:
                    message = "error"
                    print("Not in correct json format.")
                    continue
                for blast in blasts:
                    search = blast["report"]["results"]["search"]
                    title = re.search(
                        "(.), (\d+)-(\d+)", search["query_title"])
                    if title:
                        curr_strand = title.group(1)
                        curr_left = title.group(2)
                        curr_right = title.group(3)
                        results = []
                        hits = search["hits"]
                        for hit in hits:
                            hsps = hit["hsps"][0]
                            if hsps["evalue"] <= E_VALUE_THRESH:
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
                                results.append(alignment)
                        counter += 1
                        blast_result = Blast_Results(phage_id = phage_id,
                                                    id = counter,
                                                    left = curr_left,
                                                    right = curr_right,
                                                    strand = curr_strand,
                                                    results = str(results))
                        try:
                            db.session.add(blast_result)
                            db.session.commit()
                        except:
                            print("An error occured while adding a blast result. Probably a duplicate.")
                            db.session.query(Blast_Results).filter_by(phage_id=phage_id).delete()
                            db.session.commit()
                            return("error")
    except:
        print("An error occured while parsing blast results.")
        db.session.query(Blast_Results).filter_by(phage_id=phage_id).delete()
        db.session.commit()
        return("error")
        
    for filename in os.listdir(UPLOAD_FOLDER):
        if filename.endswith('.json'):
            os.remove(os.path.join(UPLOAD_FOLDER, filename))
            
    print(datetime.now())
    return(message)

def add_cds(request, UPLOAD_FOLDER, phage_id):
    """Adds a new CDS to the database.

    Checks to see if the CDS is an ORF.

    Args:
        request:
            The data sent from the front-end.
        UPLOAD_FOLDER:
            The folder containing all of the uploaded files.
        phage_id:
            The current user's ID.
    
    Returns:
        A dictionary containing a success or fail message.
    
    
    """
    new_cds_data = request.get_json()
    force = new_cds_data.get('force')
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
    frame, status = helper.get_frame_and_status(int(new_cds_data.get('left')), int(new_cds_data.get('right')), new_cds_data.get('strand'), coding_potential)
    cds = Annotations(phage_id = phage_id,
                    id = new_cds_data.get('id'),
                    left = new_cds_data.get('left'),
                    right = new_cds_data.get('right'),
                    strand = new_cds_data.get('strand'),
                    function = "None selected",
                    status = status,
                    frame = frame)

    exists = Annotations.query.filter_by(phage_id=phage_id).filter_by(left=new_cds_data.get('left'), right=new_cds_data.get('right'), strand=new_cds_data.get('strand')).first()
    orf = Blast_Results.query.filter_by(phage_id=phage_id).filter_by(left=new_cds_data.get('left'), right=new_cds_data.get('right'), strand=new_cds_data.get('strand')).first()
    if force:
        db.session.add(cds)
        db.session.commit()
        response_object['message'] = "Added succesfully."
        id_index = 0
        for cds in db.session.query(Annotations).filter_by(phage_id=phage_id).order_by(Annotations.left):
            id_index += 1
            cds.id = str(id_index)
        db.session.commit()
        id_index = 0
        for cds in db.session.query(Annotations).filter_by(phage_id=phage_id).order_by(Annotations.left):
            id_index += 1
            cds.id = phage_id + '_' + str(id_index)
        db.session.commit()
        return response_object
    if exists:
        response_object['message'] = "ID already exists."
    elif not orf:
        response_object['message'] = "Not orf."
    else:
        db.session.add(cds)
        db.session.commit()
        response_object['message'] = "Added succesfully."
        id_index = 0
        for cds in db.session.query(Annotations).filter_by(phage_id=phage_id).order_by(Annotations.left):
            id_index += 1
            cds.id = str(id_index)
        db.session.commit()
        id_index = 0
        for cds in db.session.query(Annotations).filter_by(phage_id=phage_id).order_by(Annotations.left):
            id_index += 1
            cds.id = phage_id + '_' + str(id_index)
        db.session.commit()
    return response_object