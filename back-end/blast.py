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
    response_object["uploaded"] = True
    response_object["annotated"] = False
    response_object["annotation_in_progress"] = False
    if db.session.query(Tasks).filter_by(phage_id=current_user).filter_by(function="auto_annotate").first() is not None:
        response_object["annotation_in_progress"] = True
    for filename in os.listdir(os.path.join(ROOT, 'users', current_user, 'uploads')):
        if (filename.endswith('.gdata')):
            response_object["annotated"] = True
            break
    if (db.session.query(Blast_Results).filter_by(phage_id=current_user).first() is None):
        response_object["uploaded"] = False
    for filename in os.listdir(os.path.join(ROOT, 'users', current_user)):
        if filename.endswith('.zip'):
            response_object["blast_downloaded"] = True
    return response_object

def auto_annotate(UPLOAD_FOLDER, current_user):
    """Auto annotates the entire genome.

    Args:
        UPLOAD_FOLDER:
            The path to the folder containing all files associated with the current user.
        current_user:
            The ID of the current user.
    """

    fasta_file_path = helper.get_file_path("fasta", UPLOAD_FOLDER)
    
    result = subprocess.run(["/genemark_suite_linux_64/gmsuite/gc", fasta_file_path], stdout=subprocess.PIPE)
    gc_percent = result.stdout.decode("utf-8").split(" ")[3]
    gc_percent = "{:d}".format(round(float(gc_percent)))

    print("\n\nGeneMark\n")
    subprocess.run(["/genemark_suite_linux_64/gmsuite/gm", "-D", "-g", "0", "-s", "1", "-m", "/genemark_suite_linux_64/gmsuite/heuristic_mat/heu_11_{}.mat".format(gc_percent), "-v", fasta_file_path])

    gdata_file_path = "{}.gdata".format(fasta_file_path)
    ldata_file_path = "{}.ldata".format(fasta_file_path)

    print("\n\nGlimmer\n")
    glimmer_file = os.path.join(UPLOAD_FOLDER, current_user + "_glimmer")
    subprocess.run(["/opt/glimmer3.02/scripts/g3-from-scratch.csh", fasta_file_path, glimmer_file])

    print("\n\nAragorn\n")
    aragorn_file = os.path.join(UPLOAD_FOLDER, current_user + "_aragorn.txt")
    subprocess.run(["aragorn", fasta_file_path, aragorn_file])
    
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
    Annotations.query.filter_by(phage_id=current_user).delete()
    db.session.commit()

    print("\n\nPhanotate\n")
    phanotate_file = os.path.join(UPLOAD_FOLDER, current_user + "_phanotate.txt")
    subprocess.run(["/usr/local/lib/python3.7/site-packages/phanotate.py", fasta_file_path, "-o" + phanotate_file], timeout=1000)

    id_index = 0
    glimmer_calls = ""
    with open(os.path.join(UPLOAD_FOLDER, current_user + "_glimmer.predict"), 'r') as glimmer:
        for line in glimmer:
            if line[0] == '>':
                continue
            id_index += 1
            column = line.strip().split()
            frame = int(column[3])
            left = int(column[1])
            right = int(column[2])
            strand = '+'
            if frame < 0:
                right = int(column[1])
                left = int(column[2])
                strand = '-'
            frame, status = helper.get_frame_and_status(left, right, strand, coding_potential)
            cds = Annotations(phage_id = current_user,
                            id = 'Glimmer_' + str(id_index),
                            left = left,
                            right = right,
                            strand = strand,
                            function = "None selected",
                            status = status,
                            frame = frame)
            db.session.add(cds)
            glimmer_calls += (str(left) + '-' + str(right) + ' ' + strand + ',')
    calls = Gene_Calls(phage_id = current_user,
                        id = 'Glimmer',
                        calls = glimmer_calls)
    db.session.add(calls)
    db.session.commit()
    with open(os.path.join(UPLOAD_FOLDER, current_user + "_aragorn.txt"), 'r') as aragorn:
        pattern1 = re.compile("tRNA-.*")
        pattern2 = re.compile("Sequence (.*)\[(.*),(.*)\]")
        function = ""
        right = ""
        left = ""
        strand = ""
        id_index = 0
        for line in aragorn:
            line = line.strip()
            if re.match(pattern1, line) != None:
                function = line
            if re.match(pattern2, line) != None:
                id_index += 1
                orf = re.search(pattern2, line)
                left = orf.group(2)
                right = orf.group(3)
                if (orf.group(1)) == 'c':
                    strand = '-'
                else:
                    strand = '+'
                cds = Annotations(phage_id = current_user,
                                id = "Aragorn_" + str(id_index),
                                left = int(left),
                                right = int(right),
                                strand = strand,
                                function = function,
                                status = "tRNA",
                                frame = 0)
                db.session.add(cds)
    db.session.commit()
    genemark_calls = ""
    with open(os.path.join(UPLOAD_FOLDER, current_user + ".fasta.ldata"), 'r') as genemark:
        left = 0
        right = 0
        frame = 0
        prob = 0.0
        id_index = 0
        first = True
        for line in genemark:
            if line == '\n':
                break
            if line[0] != '#':
                column = line.strip().split()
                new_prob = 0.0
                try:
                    new_prob = (float(column[3]) + float(column[4])) / 2
                except:
                    new_prob = 0.0
                if frame > 3 and left == int(column[0]) and new_prob < prob:
                    continue
                if frame < 3 and right == int(column[1]) and new_prob < prob:
                    continue
                if (first) or (frame > 3 and left == int(column[0])) or (frame < 3 and right == int(column[1])):
                    first = False
                    left = int(column[0])
                    right = int(column[1])
                    frame = int(column[2])
                    prob = new_prob
                else:
                    found = False
                    if frame > 3:
                        genemark_calls += (str(left) + '-' + str(right) + ' ' + '-' + ',')
                        cds = Annotations.query.filter_by(phage_id=current_user).filter_by(left=left).first()
                        if cds == None or cds.status == 'Fail':
                            frame, status = helper.get_frame_and_status(left, right, '-', coding_potential)
                            if cds == None or (cds.status == 'Fail' and status == 'Pass'):
                                if cds != None:
                                    db.session.delete(cds)
                                id_index += 1
                                cds = Annotations(phage_id = current_user,
                                                id = "Genemark_" + str(id_index),
                                                left = int(left),
                                                right = int(right),
                                                strand = '-',
                                                function = "None selected",
                                                status = status,
                                                frame = frame)
                            db.session.add(cds)
                            db.session.commit()
                    else:
                        genemark_calls += (str(left) + '-' + str(right) + ' ' + '+' + ',')
                        cds = Annotations.query.filter_by(phage_id=current_user).filter_by(right=right).first()
                        if cds == None or cds.status == 'Fail':
                            frame, status = helper.get_frame_and_status(left, right, '+', coding_potential)
                            if cds == None or (cds.status == 'Fail' and status == 'Pass'):
                                if cds != None:
                                    db.session.delete(cds)
                                id_index += 1
                                cds = Annotations(phage_id = current_user,
                                                id = "Genemark_" + str(id_index),
                                                left = int(left),
                                                right = int(right),
                                                strand = '+',
                                                function = "None selected",
                                                status = status,
                                                frame = frame)
                                db.session.add(cds)
                                db.session.commit()
                    left = int(column[0])
                    right = int(column[1])
                    frame = int(column[2])
                    prob = new_prob
    calls = Gene_Calls(phage_id = current_user,
                       id = 'GeneMark',
                       calls = genemark_calls)
    db.session.add(calls)
    db.session.commit()

    phanotate_calls = ""
    id_index = 0
    with open(phanotate_file, 'r') as phanotate:
        for line in phanotate:
            if line[0] == '#':
                continue
            id_index += 1
            column = line.strip().split()
            strand = str(column[2])
            left = int(column[1])
            right = int(column[0])
            cds = Annotations.query.filter_by(phage_id=current_user).filter_by(left=left).first()
            if strand == '+':
                left = int(column[0])
                right = int(column[1])
                cds = Annotations.query.filter_by(phage_id=current_user).filter_by(right=right).first()
            phanotate_calls += (str(left) + '-' + str(right) + ' ' + strand + ',')
            if cds == None or cds.status == 'Fail':
                frame, status = helper.get_frame_and_status(left, right, strand, coding_potential)
                if cds == None or (cds.status == 'Fail' and status == 'Pass'):
                    if cds != None:
                        db.session.delete(cds)
                    cds = Annotations(phage_id = current_user,
                                    id = 'Phanotate_' + str(id_index),
                                    left = left,
                                    right = right,
                                    strand = strand,
                                    function = "None selected",
                                    status = status,
                                    frame = frame)
                    db.session.add(cds)
    calls = Gene_Calls(phage_id = current_user,
                       id = 'Phanotate',
                       calls = phanotate_calls)
    db.session.add(calls)
    db.session.commit()
    id_index = 0

    for cds in db.session.query(Annotations).filter_by(phage_id=current_user).order_by(Annotations.left):
        id_index += 1
        cds.id = current_user + '_' + str(id_index)
        if (cds.right < cds.left) or (Blast_Results.query.filter_by(phage_id=current_user).filter_by(left=cds.left, right=cds.right, strand=cds.strand).first()):
            db.session.delete(cds)
    db.session.commit()

    # Remove files created for autoannotation.
    for filename in os.listdir(UPLOAD_FOLDER):
        if not filename.endswith(".fasta") and not filename.endswith(".gdata") and not filename.endswith(".json"):
            os.remove(os.path.join(UPLOAD_FOLDER, filename))

def download_blast_input(current_user):
    """Returns the blast input zip folder.

    Args:
        current_user:
            The current user ID.

    Returns:
        The blast input files in a zip folder.
    """
    f = open(os.path.join(ROOT, 'users', current_user, f"{current_user}_blast.zip"), "rb")

    return f.read()

def create_blast_input(UPLOAD_FOLDER, current_user):
    """Creates fasta file(s) for BLAST input.

    If more than one file is created, then each file should be 30 kb.

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
    fasta_file = helper.get_file_path("fasta", UPLOAD_FOLDER)
    filename = re.search('(.*/users/.*)/uploads/.*.\w*', fasta_file)
    genome = SeqIO.read(fasta_file, "fasta").seq
    output = ""
    blast_file_count = 1  # keep track of num blast files created
    out_file = f"{str(filename.group(1))}/{current_user}_blast_{blast_file_count}.fasta"
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
            out_file = f"{str(filename.group(1))}/{current_user}_blast_{blast_file_count}.fasta"
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
            out_file = f"{str(filename.group(1))}/{current_user}_blast_{blast_file_count}.fasta"
            files_to_zip.append(out_file)

    with open(out_file, "w") as f:
        f.write(output)
    
    # zip all out_files together.
    zip_file = zipfile.ZipFile(f"{str(filename.group(1))}/{current_user}_blast.zip", 'w', zipfile.ZIP_DEFLATED)
    for filename in files_to_zip:
        arcname = filename.rsplit('/', 1)[-1].lower()
        zip_file.write(filename, arcname)
    zip_file.close()
    
    # delete files that are not in zip folder.
    for filename in os.listdir(os.path.join(ROOT, 'users', current_user)):
        if filename.endswith(".fasta"):
            os.remove(os.path.join(ROOT, 'users', current_user, filename))

    return blast_file_count

def dropzone(current_user, UPLOAD_FOLDER, request):
    """Adds the blast output file to the upload directory if of type json.

    Args:
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
                        file_data = db.session.query(Files).filter_by(phage_id=current_user).filter_by(name=file_name).first()
                        file_data.complete = True
                        db.session.commit()
        if not found:
            with open(os.path.join(UPLOAD_FOLDER, file_name), 'w') as f:
                f.write(contents)
            if contents[-6:] == "\n]\n}\n\n":
                print("yes")
                file_data = db.session.query(Files).filter_by(phage_id=current_user).filter_by(name=file_name).first()
                file_data.complete = True
                db.session.commit()

def get_blast_output_names(current_user, UPLOAD_FOLDER, type_of_call):
    """Gets the names of all the files of type json in the upload directory.

    Args:
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
            file_data = db.session.query(Files).filter_by(phage_id=current_user).filter_by(name=file).first()
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
            delete_file = db.session.query(Files).filter_by(phage_id=current_user).filter_by(name=file_name).first()
            db.session.delete(delete_file)
    response_object["bad_files"] = bad_files
    response_object["file_names"] = file_names
    response_object["file_sizes"] = file_sizes
    response_object["file_mods"] = file_mods
    response_object["in_process"] = False
    response_object["position"] = -1
    task = db.session.query(Tasks).filter_by(phage_id=current_user).filter_by(function="auto_annotate").first()
    if (task is not None):
        curr_tasks = db.session.query(Tasks).filter_by(complete=False).order_by(Tasks.time)
        counter = 0
        for curr_task in curr_tasks:
            if curr_task.phage_id == current_user:
                break
            counter += 1
        response_object["position"] = counter
        response_object["in_process"] = True
        if (task.complete):
            response_object["in_process"] = False
            db.session.delete(task)
            db.session.commit()

    return response_object

def delete_blast_output(current_user, UPLOAD_FOLDER, file_path):
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
        delete_file = db.session.query(Files).filter_by(phage_id=current_user).filter_by(name=file_path).first()
        db.session.delete(delete_file)
        response_object["status"] = "success"
        db.session.query(Blast_Results).filter_by(phage_id=current_user).delete()
        try:
            db.session.commit()
        except:
            print("error in clearing table")
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