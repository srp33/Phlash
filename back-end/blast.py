"""Contains the functions for the Blast page.

Returns CDS data.
Updates CDS data.
Parses BLAST results.

Attributes:
    response_object:
        The dictionary that is returned by the main functions.
    ROOT:
        The root directory.
"""
from werkzeug.utils import secure_filename
from contextlib import closing
from zipfile import ZipFile
import time
from Bio import SeqIO, Seq, SeqFeature
from collections import OrderedDict
import subprocess
from models import *
import json
import os
import pandas as pd
import re
import zipfile
from sys import getsizeof
from datetime import datetime

ROOT = os.path.dirname(os.path.abspath(__file__))

# ------------------------------ MAIN FUNCTIONS ------------------------------
def find_blast_zip(phage_id):
    """Finds if the blast zip exists.

    Args:
        phage_id:
            The current user ID.
    
    Returns:
        A dictionary containing download boolean indicator.
    """
    response_object = {}
    response_object["blast_downloaded"] = False
    response_object["uploaded"] = True
    response_object["annotated"] = False
    response_object["annotation_in_progress"] = False
    response_object["blast_input_in_progress"] = False
    if db.session.query(Tasks).filter_by(phage_id=phage_id).filter_by(function="auto_annotate").first() is not None:
        response_object["annotation_in_progress"] = True
    for filename in os.listdir(os.path.join(ROOT, 'users', phage_id, 'uploads')):
        if (filename.endswith('.gdata')):
            response_object["annotated"] = True
            break
    if (db.session.query(Blast_Results).filter_by(phage_id=phage_id).first() is None):
        response_object["uploaded"] = False
    if db.session.query(Tasks).filter_by(phage_id=phage_id).filter_by(function="blast_input").first() is not None:
        response_object["blast_input_in_progress"] = True
    for filename in os.listdir(os.path.join(ROOT, 'users', phage_id)):
        if filename.endswith('.zip'):
            response_object["blast_downloaded"] = True
            break
    return response_object

def download_blast_input(phage_id):
    """Returns the blast input zip folder.

    Args:
        phage_id:
            The current user ID.

    Returns:
        The blast input files in a zip folder.
    """
    f = open(os.path.join(ROOT, 'users', phage_id, f"{phage_id}_blast.zip"), "rb")

    return f.read()

def dropzone(phage_id, UPLOAD_FOLDER, request):
    """Adds the blast output file to the upload directory if of type json.

    Args:
        phage_id:
            The ID of the current user.
        UPLOAD_FOLDER:
            The folder containing all the uploaded files.
        request:
            A dictionary containing the files to be uploaded.
    """
    file = request.files['file']
    contents = str(file.read(), 'utf-8')
    print(request.files)
    if file:
        file_name = secure_filename(file.filename)
        print(file_name)
        found = False
        for existing_file in os.listdir(UPLOAD_FOLDER):
            if existing_file.endswith(file_name):
                found = True
                with open(os.path.join(UPLOAD_FOLDER, existing_file), 'a+') as f:
                    f.write(contents)
                    if contents[-6:] == "\n]\n}\n\n":
                        print("yes")
                        file_data = db.session.query(Files).filter_by(phage_id=phage_id).filter_by(name=file_name).first()
                        file_data.complete = True
                        db.session.commit()
        if not found:
            with open(os.path.join(UPLOAD_FOLDER, file_name), 'w') as f:
                f.write(contents)
            if contents[-6:] == "\n]\n}\n\n":
                print("yes")
                file_data = db.session.query(Files).filter_by(phage_id=phage_id).filter_by(name=file_name).first()
                file_data.complete = True
                db.session.commit()

def get_blast_output_names(phage_id, UPLOAD_FOLDER, type_of_call):
    """Gets the names of all the files of type json in the upload directory.

    Args:
        phage_id:
            The ID of the current user.
        UPLOAD_FOLDER:
            The folder containing all the uploaded files.
        type_of_call:
            A string indicating if the if this function is called from vue or from dropzone.

    Returns:
        A dictionary containing a list of the blast output file names.
    """
    response_object = {}
    file_names = []
    file_sizes = []
    file_mods = []
    bad_files = []
    for file in os.listdir(UPLOAD_FOLDER):
        if file.endswith(".json"):
            file_data = db.session.query(Files).filter_by(phage_id=phage_id).filter_by(name=file).first()
            if file_data != None and file_data.complete:
                file_mods.append(file_data.date)
                file_sizes.append(file_data.size)
                file_names.append(file_data.name)
            elif file_data != None:
                print(file)
                bad_files.append(file_data.name)
            else:
                os.remove(os.path.join(UPLOAD_FOLDER, file))
    if type_of_call == "refresh":
        for file_name in bad_files:
            os.remove(os.path.join(UPLOAD_FOLDER, file_name))
            delete_file = db.session.query(Files).filter_by(phage_id=phage_id).filter_by(name=file_name).first()
            db.session.delete(delete_file)
    response_object["bad_files"] = bad_files
    response_object["file_names"] = file_names
    response_object["file_sizes"] = file_sizes
    response_object["file_mods"] = file_mods
    response_object["in_process"] = False
    response_object["position"] = -1
    response_object["result"] = "not complete"
    response_object["blast_input_complete"] = False
    task = db.session.query(Tasks).filter_by(phage_id=phage_id).filter_by(function="auto_annotate").first()
    if (task is not None):
        curr_tasks = db.session.query(Tasks).filter_by(complete=False).order_by(Tasks.time)
        counter = 0
        for curr_task in curr_tasks:
            if curr_task.phage_id == phage_id:
                break
            counter += 1
        response_object["position"] = counter
        response_object["in_process"] = True
        if (task.complete):
            response_object["in_process"] = False
            response_object["result"] = task.result
            db.session.delete(task)
            db.session.commit()
    task = db.session.query(Tasks).filter_by(phage_id=phage_id).filter_by(function="blast_input").first()
    if (task and task.complete):
        response_object["blast_input_complete"] = True
        response_object["num_files"] = task.result
        db.session.delete(task)
        db.session.commit()
    return response_object

def delete_blast_output(phage_id, UPLOAD_FOLDER, file_path):
    """Removes a file from the upload directory given the file path.
    
    Args:
        phage_id:
            The ID of the current user.
        UPLOAD_FOLDER:
            The folder containing all the uploaded files.
        file_path:
            The path of the file to be removed.

    Returns:
        A dictionary containing a success message.
    """
    response_object = {}
    try:
        os.remove(os.path.join(UPLOAD_FOLDER, file_path))
        delete_file = db.session.query(Files).filter_by(phage_id=phage_id).filter_by(name=file_path).first()
        db.session.delete(delete_file)
        response_object["status"] = "success"
        db.session.query(Blast_Results).filter_by(phage_id=phage_id).delete()
        try:
            db.session.commit()
        except:
            print("error in clearing table")
    except:
        print("error in deleting " + file_path)
        response_object["status"] = "error in deleting " + file_path

    return response_object

def delete_all_blast(phage_id, UPLOAD_FOLDER):
    """Deletes all data associated with the BLAST results.

    Args:
        phage_id:
            The ID of the current user.
        UPLOAD_FOLDER:
            The folder containing all the uploaded files.
    """
    if db.session.query(Tasks).filter_by(phage_id=phage_id).filter_by(function="parse_blast").first() is None:
        db.session.query(Blast_Results).filter_by(phage_id=phage_id).delete()
        db.session.query(Files).filter_by(phage_id=phage_id).delete()
        db.session.commit()
        for filename in os.listdir(UPLOAD_FOLDER):
            if filename.endswith('.json'):
                os.remove(os.path.join(UPLOAD_FOLDER, filename))
        return "success"
    else:
        return "fail"

def new_file(phage_id, file_path, file_method):
    """Adds the new file information to database.

    Args:
        phage_id:
            The ID of the current user.
    """
    index = file_path.find(".json")
    name = secure_filename(file_path[0:index + 5])
    size = file_path[index + 5:]
    file_data = Files(phage_id=phage_id,
                    name=name,
                    date=file_method,
                    size=size,
                    complete=False)
    try:
        db.session.add(file_data)
        db.session.commit()
    except:
        return "already added"
    return "success"

def add_annotation_task(phage_id, UPLOAD_FOLDER):
    """Adds task to database to be executed.

    Args:
        phage_id:
            The ID of the current user.
        UPLOAD_FOLDER:
            The folder containing all the uploaded files.
    """
    args = UPLOAD_FOLDER + " " + phage_id
    task = Tasks(phage_id=phage_id,
                function="auto_annotate",
                arguments=args,
                complete=False,
                result="waiting",
                time=str(datetime.now()))
    try:
        db.session.add(task)
        db.session.commit()
    except:
        return "Error in adding task to queue"
    return "success"

def add_blast_input_task(UPLOAD_FOLDER, phage_id):
    """Adds task to database to be executed.

    Args:
        phage_id:
            The ID of the current user.
        UPLOAD_FOLDER:
            The folder containing all the uploaded files.
    """
    args = UPLOAD_FOLDER + " " + phage_id
    task = Tasks(phage_id=phage_id,
                function="blast_input",
                arguments=args,
                complete=False,
                result="waiting",
                time='0' + str(datetime.now()))
    try:
        db.session.add(task)
        db.session.commit()
    except:
        return "Error in adding task to queue"
    return "success"


def get_num_blast_files(phage_id):
    """Gets the number of Blast input files in the zip folder.

    Args:
        phage_id:
            The ID of the current user.
    
    Returns:
        A string containing the number of Blast files or 'None' if not found.
    """
    for filename in os.listdir(os.path.join(ROOT, 'users', phage_id)):
        if filename.endswith('.zip'):
            with closing(ZipFile(os.path.join(ROOT, 'users', phage_id, filename))) as archive:
                num_blast_files = len(archive.infolist())
                return str(num_blast_files)
    return "None"