"""Contains the functions for the Upload page.

Attributes:
    response_object:
        The dictionary that the main functions return.
    ROOT:
        The root directory.
    FASTA_EXTENSIONS:
        The allowed fasta extensions.
    GENBANK_EXTENSIONS:
        The allowed genbank extensions.
    GDATA_EXTENSIONS:
        The allowed gdata extensions.
    LDATA_EXTENSIONS:
        The allowed ldata extensions.
"""
from werkzeug.utils import secure_filename
from Bio import SeqIO, Seq
from models import *
import os
from builtins import FileNotFoundError
import subprocess
import models
import helper
import shutil
import pandas as pd
import re

response_object = {}
ROOT = os.path.dirname(os.path.abspath(__file__))
FASTA_EXTENSIONS = set(['.fasta', '.fna', '.fa'])
GENBANK_EXTENSIONS = set(['.gb', '.gbk', '.gbf'])
GDATA_EXTENSIONS = set(['.gdata'])
LDATA_EXTENSIONS = set(['.ldata'])

# ------------------------------ MAIN FUNCTIONS ------------------------------
def check_uploaded_files(UPLOAD_FOLDER, phage_id):
    """Checks if fasta file is uploaded and if blast results have been added to the data base.

    Args:
        UPLOAD_FOLDER:
            The directory containing all of the uploaded files.

    Returns:
        A dictionary containing a boolean indicating if the files are uploaded.
    """
    response_object = {} 
    response_object["fasta"] = False
    for filename in os.listdir(UPLOAD_FOLDER):
        ext = os.path.splitext(filename)[1].lower()
        if ext in FASTA_EXTENSIONS:
            response_object["fasta"] = handle_fasta(UPLOAD_FOLDER)

    response_object["blast_completed"] = False if db.session.query(Blast_Results).filter_by(phage_id=phage_id).first() is None else True

    return response_object

def display_files(UPLOAD_FOLDER):
    """Finds the files (genbank or fasta) that have already been uploaded and returns their names. 

    Args:
        UPLOAD_FOLDER:
            The directory containing all of the uploaded files.

    Returns:
        A dictionary containing a success message and the genbank or fasta file names.
    """
    response_object["fasta_file"] = "Not found"
    for file in os.listdir(UPLOAD_FOLDER):
        if (file.endswith(".fasta") or file.endswith(".fna") or file.endswith(".fa")):
            response_object["fasta_file"] = file
            response_object["fasta_file_size"] = os.path.getsize(os.path.join(UPLOAD_FOLDER, file))
    return response_object

def delete_file(phage_id, file_path, UPLOAD_FOLDER):
    """Deletes a file given the file_path.

    Removes all data that has been saved associated with that file.

    Args:
        file_path:
            The path to the file to be deleted.
        UPLOAD_FOLDER:
            The directory containing all of the uploaded files.

    Returns:
        A dictionary containing a success message.
    """
    try:
        if (db.session.query(Tasks).filter_by(phage_id=phage_id).first() is None):
            if (file_path.endswith(".fasta") or file_path.endswith(".fna") or file_path.endswith(".fa")):
                db.session.query(Files).filter_by(phage_id=phage_id).delete()
                db.session.query(Annotations).filter_by(phage_id=phage_id).delete()
                db.session.query(Blast_Results).filter_by(phage_id=phage_id).delete()
                db.session.query(Gene_Calls).filter_by(phage_id=phage_id).delete()
                db.session.query(Tasks).filter_by(phage_id=phage_id).delete()
                for file in os.listdir(UPLOAD_FOLDER):
                    os.remove(os.path.join(UPLOAD_FOLDER, file))
    except:
        response_object["status"] = "error in deleting files"
    
    response_object['status'] = helper.delete_blast_zip(UPLOAD_FOLDER)

    try:
        db.session.commit()
    except:
        print("error in clearing tables")

    return response_object

def dropzone_fasta(UPLOAD_FOLDER, request, phage_id):
    """Uploads the FASTA file. 

    Args:
        UPLOAD_FOLDER:
            The directory containing all of the uploaded files.
        request:
            The file name and data from the front end.
        phage_id:
            The ID of the current user.
    """
    if (db.session.query(Tasks).filter_by(phage_id=phage_id).first() is None):
        file = request.files['file']
        contents = str(file.read(), 'utf-8')
        if file:
            with open(os.path.join(UPLOAD_FOLDER, phage_id + ".fasta"), 'a') as f:
                f.write(contents)

# ---------- FASTA FILE HELPER FUNCTIONS ----------
def handle_fasta(UPLOAD_FOLDER):
    """Ensures that the file is in FastA format.

    Args:
        UPLOAD_FOLDER:
            The directory containing all of the uploaded files.
    """
    fasta_file = helper.get_file_path("fasta", UPLOAD_FOLDER)
    with open(fasta_file, "r") as handle:
        lines = handle.readlines()
        correct_start = False
        not_empty = False
        for line in lines:
            if line.startswith('>') and not correct_start:
                correct_start = True
            elif not bool(re.match('^[ACTG\n]+$', line.upper())):
                return False
            elif len(line) > 3:
                not_empty = True
        return correct_start and not_empty