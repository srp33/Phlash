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

def add_blast_task(phage_id, UPLOAD_FOLDER):
    """Adds task to database to be executed.

    Args:
        phage_id:
            The ID of the current user.
        UPLOAD_FOLDER:
            The folder containing all the uploaded files.
    """
    if db.session.query(Blast_Results).filter_by(phage_id=phage_id).first() is None and db.session.query(Tasks).filter_by(phage_id=phage_id).filter_by(function="parse_blast").filter_by(complete=False).first() is None:
        args = phage_id + " " + UPLOAD_FOLDER
        task = Tasks(phage_id=phage_id,
                        function="parse_blast",
                        arguments=args,
                        complete=False,
                        result="waiting",
                        time=datetime.now())
        try:
            db.session.add(task)
            db.session.commit()
        except:
            return "Error in adding task to queue"
        return "empty"
    if db.session.query(Tasks).filter_by(phage_id=phage_id).filter_by(function="parse_blast").first() is not None:
        return "empty"
    else:
        return "not empty"

def check_blast_task(phage_id):
    """Checks the task's status.

    Args:
        phage_id:
            The ID of the current user.
    """
    curr_tasks = db.session.query(Tasks).filter_by(complete=False).order_by(Tasks.time)
    counter = 0
    for curr_task in curr_tasks:
        if curr_task.phage_id == phage_id:
            break
        counter += 1
    task = db.session.query(Tasks).filter_by(phage_id=phage_id).filter_by(function="parse_blast").first()
    if task is not None and task.complete:
        result = task.result
        db.session.delete(task)
        db.session.commit()
        return result
    elif task is None:
        return "complete"
    else:
        return str(counter)

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
        # id_index = 0
        # for cds in db.session.query(Annotations).filter_by(phage_id=phage_id).order_by(Annotations.left):
        #     id_index += 1
        #     cds.id = str(id_index)
        # db.session.commit()
        id_index = 0
        phage_name = db.session.query(Users).filter_by(id=phage_id).first().phage_id
        for cds in db.session.query(Annotations).filter_by(phage_id=phage_id).order_by(Annotations.left):
            id_index += 1
            cds.id = phage_name + '_' + str(id_index)
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
        # id_index = 0
        # for cds in db.session.query(Annotations).filter_by(phage_id=phage_id).order_by(Annotations.left):
        #     id_index += 1
        #     cds.id = str(id_index)
        # db.session.commit()
        id_index = 0
        phage_name = db.session.query(Users).filter_by(id=phage_id).first().phage_id
        for cds in db.session.query(Annotations).filter_by(phage_id=phage_id).order_by(Annotations.left):
            id_index += 1
            cds.id = phage_name + '_' + str(id_index)
        db.session.commit()
    return response_object