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
import helper
from datetime import datetime

response_object = {}
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
    response_object["blast_downloaded"] = False
    response_object["uploaded"] = True
    response_object["annotated"] = False
    response_object["annotation_in_progress"] = False
    if db.session.query(Tasks).filter_by(phage_id=phage_id).filter_by(function="auto_annotate").first() is not None:
        response_object["annotation_in_progress"] = True
    for filename in os.listdir(os.path.join(ROOT, 'users', phage_id, 'uploads')):
        if (filename.endswith('.gdata')):
            response_object["annotated"] = True
            break
    if (db.session.query(Blast_Results).filter_by(phage_id=phage_id).first() is None):
        response_object["uploaded"] = False
    for filename in os.listdir(os.path.join(ROOT, 'users', phage_id)):
        if filename.endswith('.zip'):
            response_object["blast_downloaded"] = True
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

def create_blast_input(UPLOAD_FOLDER, phage_id):
    """Creates fasta file(s) for BLAST input.

    If more than one file is created, then each file should be 30 kb.

    Args:
        UPLOAD_FOLDER:
            The path to the directory containing files for the current user.
        phage_id:
            The ID of the current user.

    Returns:
        The Number of blast fasta files created.
    """
    fasta_file = helper.get_file_path("fasta", UPLOAD_FOLDER)
    filename = re.search('(.*/users/.*)/uploads/.*.\w*', fasta_file)
    genome = SeqIO.read(fasta_file, "fasta").seq
    output = ""
    blast_file_count = 1  # keep track of num blast files created
    out_file = f"{str(filename.group(1))}/{phage_id}_blast_{blast_file_count}.fasta"
    files_to_zip = [out_file]
    lefts, rights = get_starts_stops('+', genome)
    for i in range(len(lefts)):
        output += f">+, {lefts[i]}-{rights[i]}\n"
        output += f"{Seq.translate(sequence=helper.get_sequence(genome, '+', lefts[i]-1, rights[i]), table=11)}\n"
        if getsizeof(output) > 30000:  # only file size to reach 30 kb, else you get a CPU limit from NCBI blast
            with open(out_file, "w") as f:
                f.write(output)
            output = ""
            blast_file_count += 1
            out_file = f"{str(filename.group(1))}/{phage_id}_blast_{blast_file_count}.fasta"
            files_to_zip.append(out_file)
    lefts, rights = get_starts_stops('-', genome)
    for i in range(len(lefts)):
        output += f">-, {lefts[i]}-{rights[i]}\n"
        left = len(genome) - rights[i]
        right = len(genome) - lefts[i] + 1
        output += f"{Seq.translate(sequence=helper.get_sequence(genome, '-', left, right), table=11)}\n"
        if getsizeof(output) > 30000:  # only file size to reach 30 kb, else you get a CPU limit from NCBI blast
            with open(out_file, "w") as f:
                f.write(output)
            output = ""
            blast_file_count += 1
            out_file = f"{str(filename.group(1))}/{phage_id}_blast_{blast_file_count}.fasta"
            files_to_zip.append(out_file)

    with open(out_file, "w") as f:
        f.write(output)
    
    # zip all out_files together.
    zip_file = zipfile.ZipFile(f"{str(filename.group(1))}/{phage_id}_blast.zip", 'w', zipfile.ZIP_DEFLATED)
    for filename in files_to_zip:
        arcname = filename.rsplit('/', 1)[-1].lower()
        zip_file.write(filename, arcname)
    zip_file.close()
    
    # delete files that are not in zip folder.
    for filename in os.listdir(os.path.join(ROOT, 'users', phage_id)):
        if filename.endswith(".fasta"):
            os.remove(os.path.join(ROOT, 'users', phage_id, filename))

    return blast_file_count

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
                time=datetime.now())
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

# ----- BLAST CREATION HELPER FUNCTIONS ------
def get_stop_options(genome, start, strand):
    """Finds the next stop codon given a start codon index.

    Args:
        genome:
            The genome of the phage.
        start:
            The location of the start codon.
        strand:
            Complimentary or direct strand.

    Returns:
        The stop codon index.
    """
    
    bacteria_stop_codons = ["TAG", "TAA", "TGA", "tag", "taa", "tga"]
    start -= 1
    gene = ""
    if (strand != "-"):
        gene = genome[start:]
    else:
        gene = genome.reverse_complement()[start:]
    for index in range(0, len(gene), 3):
        codon = gene[index:index + 3]
        if codon in bacteria_stop_codons:
            return (start + index + 3)

def get_start_options(genome, maximum, strand, minimum):
    """Finds all the start codons within a range of DNA indicated by the minimum and maximum parameters.

    Args:
        genome:
            The genome of the phage.
        maximum:
            The max index to search for start codons.
        strand:
            Complimentary or direct strand.
        minimum:
            The minimum index to search for start codons.
    """

    bacteria_start_codons = ["ATG", "GTG", "TTG", "atg", "gtg", "ttg"]
    start_options = []
    maximum += 3
    gene = ""
    if (strand != "-"):
        gene = genome[minimum:maximum]
    else:
        gene = genome.reverse_complement()[minimum:maximum]
    for index in range(0, len(gene)):
        codon = gene[index:index + 3]
        if codon in bacteria_start_codons:
            start_options.append(minimum + index + 1)
    return start_options

def get_starts_stops(strand, genome):
    """Finds alternative, possible start and stop positions for a given CDS. 
    
    Args:
        cds_id:
            The ID of the CDS.
        genome:
            The genome of the phage.

    Returns:
        start_options, stop_options: 
            list of alternative starts and an associated list of alternative stops.
    """

    genome_length = len(genome)
    possible_start_options = get_start_options(genome, genome_length - 1, strand, 0)

    start_options = []
    stop_options = []
    for start in possible_start_options:
        stop_options.append(get_stop_options(genome, start, strand))
        if stop_options[-1] is None:
            stop_options.pop()
        elif stop_options[-1] - start < 30:
            stop_options.pop()
        else:
            if strand == '+':
                start_options.append(start)
            else:
                start_options.append(genome_length - start + 1)
                stop = stop_options.pop()
                stop_options.append(genome_length - stop + 1)

    if strand == '+':
        return start_options, stop_options
    else:
        return stop_options, start_options