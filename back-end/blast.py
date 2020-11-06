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
from datetime import datetime
from Bio import SeqIO, Seq, SeqFeature
from collections import OrderedDict
from models import *
import json
import os
import pandas as pd
import re
import zipfile
from sys import getsizeof
import helper

response_object = {}
ROOT = os.path.dirname(os.path.abspath(__file__))

# ------------------------------ MAIN FUNCTIONS ------------------------------
def find_blast_zip(current_user):
    """Finds if the blast zip exists.

    Args:
        current_user:
            The current user ID.
    
    Returns:
        A dictionary containing download boolean indicator.
    """
    response_object["blast_downloaded"] = False

    for filename in os.listdir(os.path.join(ROOT, 'users', current_user)):
        if filename.endswith('.zip'):
            response_object["blast_downloaded"] = True

    return response_object

def download_blast_input(UPLOAD_FOLDER, current_user):
    """Creates and returns the blast input zip folder.

    Args:
        UPLOAD_FOLDER:
            The folder containing all the uploaded files.
        current_user:
            The current user ID.

    Returns:
        The blast input files in a zip folder.
    """
    print("Starting comparisons")
    compare()
    fasta_file = helper.get_file_path("fasta", UPLOAD_FOLDER)
    genemark_gdata_file = helper.get_file_path("gdata", UPLOAD_FOLDER)
    print(datetime.now())
    num_files_downloaded = create_blast_fasta(current_user, fasta_file, genemark_gdata_file)
    print(datetime.now())
    f = open(os.path.join(ROOT, 'users', current_user, f"{current_user}_blast.zip"), "rb")

    return f.read()

def upload_blast_output(UPLOAD_FOLDER, request):
    """Adds the blast output file to the upload directory if of type json.

    Args:
        UPLOAD_FOLDER:
            The folder containing all the uploaded files.
        request:
            A dictionary containing the files to be uploaded.

    Returns:
        A dictionary containing a success or fail message.
    """
    response_object = {}
    if 'file' not in request.files:
        print(request.files)
        response_object["status"] = "'file' not in request.files"
        print("in fail")
    else:
        print("in success")
        file = request.files['file']
        if file:
            file_name = secure_filename(file.filename)
            if helper.allowed_file(file_name, set(['.json'])):

                file_ext = file_name.rsplit('.', 1)[1].lower()
                for existing_file in os.listdir(UPLOAD_FOLDER):
                    if existing_file.endswith(file_name):
                        os.remove(os.path.join(
                            UPLOAD_FOLDER, existing_file))

                file.save(os.path.join(UPLOAD_FOLDER, file_name))
                response_object["uploaded"] = file_name
                print(' * uploaded', file_name)
            else:
                response_object["not_allowed"] = file.filename
        else:
            response_object["status"] = "error in uploading file" 

    return response_object

def get_blast_output_names(UPLOAD_FOLDER):
    """Gets the names of all the files of type json in the upload directory.

    Args:
        UPLOAD_FOLDER:
            The folder containing all the uploaded files.

    Returns:
        A dictionary containing a list of the blast output file names.
    """
    file_names = []
    for file in os.listdir(UPLOAD_FOLDER):
        if file.endswith(".json"):
            file_names.append(file)
    response_object["file_names"] = file_names

    return response_object

def download_blast_output(UPLOAD_FOLDER, file_path):
    """Returns the contents of a file given the file path.

    Returns fail message if file not found.

    Args:
        UPLOAD_FOLDER:
            The folder containing all the uploaded files.
        file_path:
            The path to the file to be downloaded.

    Returns:
        A dictionary containing the contents of the file and a success message.
    """
    try:
        response_object["file_data"] = open(os.path.join(UPLOAD_FOLDER, file_path)).read()
        response_object["status"] = "success"
    except:
        print("error in downloading " + file_path)
        response_object["status"] = "error in downloading " + file_path

    return response_object

def delete_blast_output(UPLOAD_FOLDER, file_path):
    """Removes a file from the upload directory given the file path.
    
    Args:
        UPLOAD_FOLDER:
            The folder containing all the uploaded files.
        file_path:
            The path of the file to be removed.

    Returns:
        A dictionary containing a success message.
    """
    try:
        os.remove(os.path.join(UPLOAD_FOLDER, file_path))
        response_object["status"] = "success"
    except:
        print("error in deleting " + file_path)
        response_object["status"] = "error in deleting " + file_path

    return response_object

def get_num_blast_files(current_user):
    """Gets the number of Blast input files in the zip folder.

    Args:
        current_user:
            The ID of the current user.
    
    Returns:
        A string containing the number of Blast files or 'None' if not found.
    """
    for filename in os.listdir(os.path.join(ROOT, 'users', current_user)):
        if filename.endswith('.zip'):
            with closing(ZipFile(os.path.join(ROOT, 'users', current_user, filename))) as archive:
                num_blast_files = len(archive.infolist())
                return str(num_blast_files)
    return "None"

# ---------- BLAST HELPER FUNCTIONS ----------
def compare():
    """Compares the GeneMark and the DNAMaster CDS calls according to length. 
    """
    for cds in DNAMaster.query.all():
        dnamaster_cds = DNAMaster.query.filter_by(stop=cds.stop).first()
        genemark_cds = GeneMark.query.filter_by(stop=cds.stop).first()
        if dnamaster_cds and genemark_cds:
            if dnamaster_cds.start <= genemark_cds.start:
                dnamaster_cds.status = "Pass"
            else:
                dnamaster_cds.status = "Fail"
        elif not genemark_cds:
            dnamaster_cds.status = "Undetermined"
        db.session.commit()

def create_blast_fasta(current_user, fasta_file, gdata_file):
    """Creates fasta file(s) for BLAST input.

    If more than one file is created, then each file should be 15 kb.

    Args:
        current_user:
            The ID of the current user.
        fasta_file:
            The file containing the DNA sequence of the phage.
        gdata_file:
            The file containing the GeneMark gene calls.

    Returns:
        The Number of blast fasta files created.
    """
    filename = re.search('(.*/users/.*)/uploads/.*.\w*', fasta_file)
    genome = SeqIO.read(fasta_file, "fasta").seq
    output = ""
    blast_file_count = 1  # keep track of num blast files created
    out_file = f"{str(filename.group(1))}/{current_user}_blast_{blast_file_count}.fasta"
    files_to_zip = [out_file]
    for cds in db.session.query(DNAMaster).order_by(DNAMaster.start):
        starts, stops = get_starts_stops(cds.id, genome, gdata_file)
        cds.start_options = ", ".join([str(start) for start in starts])
        cds.stop_options = ", ".join([str(stop) for stop in stops])
        db.session.commit()
        for i in range(len(starts)):
            if cds.strand == '+':
                output += f">{cds.id}_{i + 1}, {starts[i]}-{stops[i]}\n"
                output += f"{Seq.translate(sequence=helper.get_sequence(genome, cds.strand, starts[i]-1, stops[i]), table=11)}\n"
            else:
                output += f">{cds.id}_{i + 1}, {starts[i]}-{stops[i]}\n"
                start = len(genome) - stops[i]
                stop = len(genome) - starts[i] + 1
                output += f"{Seq.translate(sequence=helper.get_sequence(genome, cds.strand, start, stop), table=11)}\n"
            if getsizeof(output) > 15000:  # only file size to reach 15 kb, else you get a CPU limit from NCBI blast
                with open(out_file, "w") as f:
                    f.write(output)
                output = ""
                blast_file_count += 1
                out_file = f"{str(filename.group(1))}/{current_user}_blast_{blast_file_count}.fasta"
                files_to_zip.append(out_file)

    with open(out_file, "w") as f:
        f.write(output)
    
    # zip all out_files together
    zip_file = zipfile.ZipFile(f"{str(filename.group(1))}/{current_user}_blast.zip", 'w', zipfile.ZIP_DEFLATED)
    for filename in files_to_zip:
        arcname = filename.rsplit('/', 1)[-1].lower()
        zip_file.write(filename, arcname)
    zip_file.close()
    
    return blast_file_count

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
    bacteria_stop_codons = ["TAG", "TAA", "TGA"]
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
    print(minimum)
    print(maximum)
    print(strand)
    bacteria_start_codons = ["ATG", "GTG", "TTG"]
    start_options = []
    maximum += 3
    gene = ""
    if (strand != "-"):
        gene = genome[minimum:maximum]
        print(gene)
    else:
        gene = genome.reverse_complement()[minimum:maximum]
    for index in range(0, len(gene)):
        codon = gene[index:index + 3]
        if codon in bacteria_start_codons:
            start_options.append(minimum + index + 1)
    return start_options

def get_starts_stops(cds_id, genome, genemark_gdata_file):
    """Finds alternative, possible start and stop positions for a given CDS. 
    
    Args:
        cds_id:
            The ID of the CDS.
        genome:
            The genome of the phage.
        genemark_gdata_file:
            The file containing the genemark gene calls.

    Returns:
        start_options, stop_options: 
            list of alternative starts and an associated list of alternative stops.
    """
    gdata_df = pd.read_csv(genemark_gdata_file, sep='\t', skiprows=16)
    gdata_df.columns = ['Base', '1', '2', '3', '4', '5', '6']
    gdata_df = gdata_df.set_index('Base')

    dnamaster_cds = DNAMaster.query.filter_by(id=cds_id).first()
    genemark_cds = GeneMark.query.filter_by(stop=dnamaster_cds.stop).first()

    minimum = -100 # 100 base pairs lower than the start position
    maximum = 100 # 100 base pairs higher than the start position
    genome_length = len(genome)
    if dnamaster_cds.strand == '+':
        if genemark_cds is not None:
            if (genemark_cds.start < dnamaster_cds.start):
                minimum += genemark_cds.start
                if minimum < 0:
                    minimum = 0
                maximum += dnamaster_cds.start
                if maximum >= genome_length:
                    maximum = genome_length - 1
            else:
                minimum += dnamaster_cds.start
                if minimum < 0:
                    minimum = 0
                maximum += genemark_cds.start
                if maximum >= genome_length:
                    maximum = genome_length - 1
        else:
            minimum += dnamaster_cds.start
            if minimum < 0:
                minimum = 0
            maximum += dnamaster_cds.start
            if maximum >= genome_length:
                maximum = genome_length - 1
    else:
        d_start = genome_length - dnamaster_cds.stop
        if genemark_cds is not None:
            g_start = genome_length - genemark_cds.stop
            if (g_start < d_start):
                minimum += g_start
                if minimum < 0:
                    minimum = 0
                maximum += d_start
                if maximum >= genome_length:
                    maximum = genome_length - 1
            else:
                minimum += d_start
                if minimum < 0:
                    minimum = 0
                maximum += g_start
                if maximum >= genome_length:
                    maximum = genome_length - 1
        else:
            minimum += d_start
            if minimum < 0:
                minimum = 0
            maximum += d_start
            if maximum >= genome_length:
                maximum = genome_length - 1

    possible_start_options = get_start_options(genome, maximum, dnamaster_cds.strand, minimum)

    start_options = []
    stop_options = []
    for start in possible_start_options:
        stop_options.append(get_stop_options(genome, start, dnamaster_cds.strand))
        if (stop_options[-1] - start < 75 and start != dnamaster_cds.start):
            stop_options.pop()
        else:
            if dnamaster_cds.strand == '+':
                start_options.append(start)
            else:
                start_options.append(genome_length - start + 1)
                stop = stop_options.pop()
                stop_options.append(genome_length - stop + 1)

    if len(start_options) == 0:
        start_options.append(dnamaster_cds.start)
        stop_options.append(dnamaster_cds.stop)

    print(start_options)
    print(stop_options)

    if dnamaster_cds.strand == '+':
        return start_options, stop_options
    else:
        return stop_options, start_options

