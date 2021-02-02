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

response_object = {}
ROOT = os.path.dirname(os.path.abspath(__file__))
FASTA_EXTENSIONS = set(['.fasta', '.fna'])
GENBANK_EXTENSIONS = set(['.gb', '.gbk', '.gbf'])
GDATA_EXTENSIONS = set(['.gdata'])
LDATA_EXTENSIONS = set(['.ldata'])

# ------------------------------ MAIN FUNCTIONS ------------------------------
def check_uploaded_files(UPLOAD_FOLDER):
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
            response_object["fasta"] = True

    response_object["blast_completed"] = False if db.session.query(Blast_Results).first() is None else True

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
        if (file.endswith(".fasta") or file.endswith(".fna")):
            response_object["fasta_file"] = file
            response_object["fasta_file_size"] = os.path.getsize(os.path.join(UPLOAD_FOLDER, file))
    return response_object

def delete_file(file_path, UPLOAD_FOLDER):
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
        if (file_path.endswith(".fasta") or file_path.endswith(".fna")):
            db.session.query(Files).delete()
            db.session.query(DNAMaster).delete()
            db.session.query(Blast_Results).delete()
            db.session.query(Gene_Calls).delete()
            for file in os.listdir(UPLOAD_FOLDER):
                os.remove(os.path.join(UPLOAD_FOLDER, file))
    except:
        print("error")
        response_object["status"] = "error in deleting files"
    
    response_object['status'] = helper.delete_blast_zip(UPLOAD_FOLDER)

    try:
        db.session.commit()
    except:
        print("error in clearing tables")

    return response_object

def dropzone_fasta(UPLOAD_FOLDER, request, current_user):
    file = request.files['file']
    contents = str(file.read(), 'utf-8')
    print(request.files)
    if file:
        with open(os.path.join(UPLOAD_FOLDER, current_user + ".fasta"), 'w') as f:
            f.write(contents)
        handle_fasta(UPLOAD_FOLDER)

# ------------------------------ UPLOAD HELPER FUNCTIONS ------------------------------
def overwrite_files(file_type, UPLOAD_FOLDER):
    """Overwrites existing files with specific extension with newly uploaded files.

    Args:
        file_type:
            The type of the file to be overwritten.
        UPLOAD_FOLDER:
            The directory containing all of the uploaded files.
    """
    for existing_file in os.listdir(UPLOAD_FOLDER):
        ext = os.path.splitext(existing_file)[1].lower()
        if (file_type == "fasta" and ext in ['.fasta', '.fna', '.ldata', '.gdata', '.lst', '.ps']) or \
            (file_type == "genbank" and ext in GENBANK_EXTENSIONS):
            os.remove(os.path.join(UPLOAD_FOLDER, existing_file))
            print(f" * removed {existing_file}")

# ---------- FASTA FILE HELPER FUNCTIONS ----------
def handle_fasta(UPLOAD_FOLDER):
    """Runs GeneMark and parses through ldata file.

    Args:
        UPLOAD_FOLDER:
            The directory containing all of the uploaded files.
    """
    fasta_file = helper.get_file_path("fasta", UPLOAD_FOLDER)
    # run_genemark(fasta_file)
    try:
        ldata_file = helper.get_file_path("ldata", UPLOAD_FOLDER)
    except FileNotFoundError:
        print("GeneMark did not save ldata file properly.")
        print("Check that your GeneMark key is in the correct location.")
        return("error")

def run_genemark(fasta_file_path):
    """Invokes the GeneMarkS utilities and adds .gdata and .ldata files.

    Args:
        fasta_file_path:
            The path to the fasta file.
    """
    result = subprocess.run(["/genemark_suite_linux_64/gmsuite/gc", fasta_file_path], stdout=subprocess.PIPE)
    gc_percent = result.stdout.decode("utf-8").split(" ")[3]
    gc_percent = "{:d}".format(round(float(gc_percent)))

    subprocess.run(["/genemark_suite_linux_64/gmsuite/gm", "-D", "-g", "0", "-s", "1", "-m", "/genemark_suite_linux_64/gmsuite/heuristic_mat/heu_11_{}.mat".format(gc_percent), "-v", fasta_file_path])

    gdata_file_path = "{}.gdata".format(fasta_file_path)
    ldata_file_path = "{}.ldata".format(fasta_file_path)

    glimmer_file = fasta_file_path[0 : -6] + "_glimmer"
    subprocess.run(["/opt/glimmer3.02/scripts/g3-from-scratch.csh", fasta_file_path, glimmer_file])

    aragorn_file = fasta_file_path[0 : -6] + "_aragorn.txt"
    subprocess.run(["aragorn", fasta_file_path, aragorn_file])

def parse_genemark_ldata(gm_file):
    """Parses through GeneMark ldata file to gather CDS calls. 

    Adds each CDS to GeneMark table in user's database. 

    Args:
        gm_file:
            The GeneMark ldata file.
    """
    genemark_all_genes = dict()
    with open(gm_file, "r") as handle:
        for line in handle:
            if line == '\n':
                break
            if line[0] != '#':
                column = line.strip().split()
                start = int(column[0])
                stop = int(column[1])
                if stop not in genemark_all_genes:
                    genemark_all_genes[stop] = [start]
                else:
                    if start not in genemark_all_genes[stop]:
                        genemark_all_genes[stop].append(start)
                curr_keys = get_keys_by_value(genemark_all_genes, start)
                if len(curr_keys) > 1:
                    max_right = max(curr_keys)
                    for key in curr_keys:
                        if key != max_right:
                            del genemark_all_genes[key]

    genemark_cdss = dict()
    with open(gm_file, "r") as handle:
        for line in handle:
            if line == '\n':
                break
            if line[0] != '#':
                column = line.strip().split()
                curr_start = int(column[0])
                curr_stop = int(column[1])
                frame = int(column[2])
                if 1 <= frame <= 3:
                    curr_frame = "+"
                elif 4 <= frame <= 6:
                    curr_frame = "-"
                id_number = 1
                for stop in genemark_all_genes:
                    min_start = min(genemark_all_genes[stop])
                    if min_start == curr_start and stop == curr_stop:
                        if (min_start, stop) not in genemark_cdss:
                            id = "genemark_" + str(id_number)
                            genemark_cdss[(min_start, stop)] = curr_frame
                            cds = GeneMark(id=id,
                                           start=min_start,
                                           stop=stop,
                                           strand=curr_frame)
                            exists = GeneMark.query.filter_by(id=id).first()
                            if not exists:
                                db.session.add(cds)
                                db.session.commit()
                    id_number += 1

def get_keys_by_value(dict, value_to_find):
    """Finds the stop location (key) given a start location (value).

    Args:
        dict:
            The dictionary containing keys and values to find.
        value_to_find:
            The value to be found in the dictionary.

    Returns:
        A list of the keys associated with the value.
    """
    keys = list()
    items = dict.items()
    for item in items:
        if value_to_find in item[1]:
            keys.append(item[0])
    return keys