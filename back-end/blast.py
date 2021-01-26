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
    response_object["annotated"] = True
    if db.session.query(DNAMaster).first() is None:
        response_object["annotated"] = False
    if (db.session.query(Blast_Results).first() is None):
        response_object["uploaded"] = False
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
    # compare()
    fasta_file = helper.get_file_path("fasta", UPLOAD_FOLDER)
    # genemark_gdata_file = helper.get_file_path("gdata", UPLOAD_FOLDER)
    print(datetime.now())
    num_files_downloaded = create_blast_fastas(current_user, fasta_file)
    print(datetime.now())
    f = open(os.path.join(ROOT, 'users', current_user, f"{current_user}_blast.zip"), "rb")

    return f.read()

def dropzone(UPLOAD_FOLDER, request):
    file = request.files['file']
    contents = str(file.read(), 'utf-8')
    print(request.files)
    if file:
        file_name = secure_filename(file.filename)
        found = False
        for existing_file in os.listdir(UPLOAD_FOLDER):
            if existing_file.endswith(file_name):
                found = True
                with open(os.path.join(UPLOAD_FOLDER, existing_file), 'a+') as f:
                    f.write(contents)
        if not found:
            #file.save(os.path.join(UPLOAD_FOLDER, file_name))
            with open(os.path.join(UPLOAD_FOLDER, file_name), 'w') as f:
                f.write(contents)

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
            print(list(bytes((file.read()[0:6]))))
            print(list(bytes(b'{\n"Bla')))
            print(1)
            file_name = secure_filename(file.filename)
            print(file_name)
            if compare(list(bytes(file.read()[0:6])), list(bytes(b'{\n"Bla'))):
                print(2)
                if helper.allowed_file(file_name, set(['.json'])):
                    print(3)
                    file_ext = file_name.rsplit('.', 1)[1].lower()
                    for existing_file in os.listdir(UPLOAD_FOLDER):
                        if existing_file.endswith(file_name):
                            os.remove(os.path.join(
                                UPLOAD_FOLDER, existing_file))

                    file.save(os.path.join(UPLOAD_FOLDER, file_name))
                    print(' * uploaded', file_name)
                else:
                    response_object["not_allowed"] = file.filename
            else:
                print(4)
                for existing_file in os.listdir(UPLOAD_FOLDER):
                    if existing_file.endswith(file_name[1:]):
                        print(5)
                        with open(existing_file, 'a+') as f1:
                            with open(file, 'r') as f2:
                                print(f1.read())
                                print(f2.read())
                                f1.write(f2.read())
                                print(f1.read())
            response_object["uploaded"] = file_name

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
    file_sizes = []
    for file in os.listdir(UPLOAD_FOLDER):
        if file.endswith(".json"):
            file_names.append(file)
            file_sizes.append(os.path.getsize(os.path.join(UPLOAD_FOLDER, file)))
    response_object["file_names"] = file_names
    response_object["file_sizes"] = file_sizes

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
        db.session.query(Blast_Results).delete()
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

# ---------- BLAST HELPER FUNCTIONS ----------

def auto_annotate(UPLOAD_FOLDER, current_user):
    fasta_file_path = helper.get_file_path("fasta", UPLOAD_FOLDER)
    
    result = subprocess.run(["/genemark_suite_linux_64/gmsuite/gc", fasta_file_path], stdout=subprocess.PIPE)
    gc_percent = result.stdout.decode("utf-8").split(" ")[3]
    gc_percent = "{:d}".format(round(float(gc_percent)))

    subprocess.run(["/genemark_suite_linux_64/gmsuite/gm", "-D", "-g", "0", "-s", "1", "-m", "/genemark_suite_linux_64/gmsuite/heuristic_mat/heu_11_{}.mat".format(gc_percent), "-v", fasta_file_path])

    gdata_file_path = "{}.gdata".format(fasta_file_path)
    ldata_file_path = "{}.ldata".format(fasta_file_path)

    glimmer_file = os.path.join(UPLOAD_FOLDER, current_user + "_glimmer")
    subprocess.run(["/opt/glimmer3.02/scripts/g3-from-scratch.csh", fasta_file_path, glimmer_file])

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
    DNAMaster.query.delete()
    db.session.commit()

    phanotate_file = os.path.join(UPLOAD_FOLDER, current_user + "_phanotate.txt")
    subprocess.run(["/usr/local/lib/python3.7/site-packages/phanotate.py", fasta_file_path, "-o" + phanotate_file])

    id_index = 0
    glimmer_calls = ""
    with open(os.path.join(UPLOAD_FOLDER, current_user + "_glimmer.predict"), 'r') as glimmer:
        for line in glimmer:
            if line[0] == '>':
                continue
            id_index += 1
            column = line.strip().split()
            frame = int(column[3])
            start = int(column[1])
            stop = int(column[2])
            strand = '+'
            if frame < 0:
                stop = int(column[1])
                start = int(column[2])
                strand = '-'
            frame, status = helper.get_frame_and_status(start, stop, strand, coding_potential)
            cds = DNAMaster(id = 'Glimmer_' + str(id_index),
                            start = start,
                            stop = stop,
                            strand = strand,
                            function = "None selected",
                            status = status,
                            frame = frame)
            db.session.add(cds)
            glimmer_calls += (str(start) + '-' + str(stop) + ' ' + strand + ',')
    calls = Gene_Calls(id = 'Glimmer',
                       calls = glimmer_calls)
    db.session.add(calls)
    db.session.commit()
    with open(os.path.join(UPLOAD_FOLDER, current_user + "_aragorn.txt"), 'r') as aragorn:
        pattern1 = re.compile("tRNA-.*")
        pattern2 = re.compile("Sequence (.*)\[(.*),(.*)\]")
        function = ""
        stop = ""
        start = ""
        strand = ""
        id_index = 0
        for line in aragorn:
            line = line.strip()
            if re.match(pattern1, line) != None:
                function = line
            if re.match(pattern2, line) != None:
                id_index += 1
                orf = re.search(pattern2, line)
                start = orf.group(2)
                stop = orf.group(3)
                if (orf.group(1)) == 'c':
                    strand = '-'
                else:
                    strand = '+'
                cds = DNAMaster(id = "Aragorn_" + str(id_index),
                                start = int(start),
                                stop = int(stop),
                                strand = strand,
                                function = function,
                                status = "tRNA",
                                frame = 0)
                db.session.add(cds)
    db.session.commit()
    genemark_calls = ""
    with open(os.path.join(UPLOAD_FOLDER, current_user + ".fasta.ldata"), 'r') as genemark:
        start = 0
        stop = 0
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
                if frame > 3 and start == int(column[0]) and new_prob < prob:
                    continue
                if frame < 3 and stop == int(column[1]) and new_prob < prob:
                    continue
                if (first) or (frame > 3 and start == int(column[0])) or (frame < 3 and stop == int(column[1])):
                    first = False
                    start = int(column[0])
                    stop = int(column[1])
                    frame = int(column[2])
                    prob = new_prob
                else:
                    found = False
                    if frame > 3:
                        genemark_calls += (str(start) + '-' + str(stop) + ' ' + '-' + ',')
                        cds = DNAMaster.query.filter_by(start=start).first()
                        if cds == None or cds.status == 'Fail':
                            frame, status = helper.get_frame_and_status(start, stop, '-', coding_potential)
                            if cds == None or (cds.status == 'Fail' and status == 'Pass'):
                                if cds != None:
                                    db.session.delete(cds)
                                id_index += 1
                                cds = DNAMaster(id = "Genemark_" + str(id_index),
                                                start = int(start),
                                                stop = int(stop),
                                                strand = '-',
                                                function = "None selected",
                                                status = status,
                                                frame = frame)
                            db.session.add(cds)
                            db.session.commit()
                    else:
                        genemark_calls += (str(start) + '-' + str(stop) + ' ' + '+' + ',')
                        cds = DNAMaster.query.filter_by(stop=stop).first()
                        if cds == None or cds.status == 'Fail':
                            frame, status = helper.get_frame_and_status(start, stop, '+', coding_potential)
                            if cds == None or (cds.status == 'Fail' and status == 'Pass'):
                                if cds != None:
                                    db.session.delete(cds)
                                id_index += 1
                                cds = DNAMaster(id = "Genemark_" + str(id_index),
                                                start = int(start),
                                                stop = int(stop),
                                                strand = '+',
                                                function = "None selected",
                                                status = status,
                                                frame = frame)
                                db.session.add(cds)
                                db.session.commit()
                    start = int(column[0])
                    stop = int(column[1])
                    frame = int(column[2])
                    prob = new_prob
    calls = Gene_Calls(id = 'GeneMark',
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
            start = int(column[1])
            stop = int(column[0])
            cds = DNAMaster.query.filter_by(start=start).first()
            if strand == '+':
                start = int(column[0])
                stop = int(column[1])
                cds = DNAMaster.query.filter_by(stop=stop).first()
            phanotate_calls += (str(start) + '-' + str(stop) + ' ' + strand + ',')
            if cds == None or cds.status == 'Fail':
                frame, status = helper.get_frame_and_status(start, stop, strand, coding_potential)
                if cds == None or (cds.status == 'Fail' and status == 'Pass'):
                    if cds != None:
                        db.session.delete(cds)
                    cds = DNAMaster(id = 'Phanotate_' + str(id_index),
                                    start = start,
                                    stop = stop,
                                    strand = strand,
                                    function = "None selected",
                                    status = status,
                                    frame = frame)
                    db.session.add(cds)
    calls = Gene_Calls(id = 'Phanotate',
                       calls = phanotate_calls)
    db.session.add(calls)
    db.session.commit()
    print(phanotate_calls)
    print(genemark_calls)
    print(glimmer_calls)
    id_index = 0
    for cds in db.session.query(DNAMaster).order_by(DNAMaster.start):
        id_index += 1
        cds.id = current_user + '_' + str(id_index)
        if (cds.stop < cds.start) or (Blast_Results.query.filter_by(start=cds.start, stop=cds.stop, strand=cds.strand).first()):
            db.session.delete(cds)
    db.session.commit()


# def create_blast_fasta(current_user, fasta_file, gdata_file):
#     """Creates fasta file(s) for BLAST input.

#     If more than one file is created, then each file should be 30 kb.

#     Args:
#         current_user:
#             The ID of the current user.
#         fasta_file:
#             The file containing the DNA sequence of the phage.
#         gdata_file:
#             The file containing the GeneMark gene calls.

#     Returns:
#         The Number of blast fasta files created.
#     """
#     filename = re.search('(.*/users/.*)/uploads/.*.\w*', fasta_file)
#     genome = SeqIO.read(fasta_file, "fasta").seq
#     output = ""
#     blast_file_count = 1  # keep track of num blast files created
#     out_file = f"{str(filename.group(1))}/{current_user}_blast_{blast_file_count}.fasta"
#     files_to_zip = [out_file]
#     for cds in db.session.query(DNAMaster).order_by(DNAMaster.start):
#         starts, stops = get_starts_stops(cds.id, genome, gdata_file)
#         cds.start_options = ", ".join([str(start) for start in starts])
#         cds.stop_options = ", ".join([str(stop) for stop in stops])
#         db.session.commit()
#         for i in range(len(starts)):
#             if cds.strand == '+':
#                 output += f">{cds.id}_{i + 1}, {starts[i]}-{stops[i]}\n"
#                 output += f"{Seq.translate(sequence=helper.get_sequence(genome, cds.strand, starts[i]-1, stops[i]), table=11)}\n"
#             else:
#                 output += f">{cds.id}_{i + 1}, {starts[i]}-{stops[i]}\n"
#                 start = len(genome) - stops[i]
#                 stop = len(genome) - starts[i] + 1
#                 output += f"{Seq.translate(sequence=helper.get_sequence(genome, cds.strand, start, stop), table=11)}\n"
#             if getsizeof(output) > 30000:  # only file size to reach 30 kb, else you get a CPU limit from NCBI blast
#                 with open(out_file, "w") as f:
#                     f.write(output)
#                 output = ""
#                 blast_file_count += 1
#                 out_file = f"{str(filename.group(1))}/{current_user}_blast_{blast_file_count}.fasta"
#                 files_to_zip.append(out_file)

#     with open(out_file, "w") as f:
#         f.write(output)
    
#     # zip all out_files together
#     zip_file = zipfile.ZipFile(f"{str(filename.group(1))}/{current_user}_blast.zip", 'w', zipfile.ZIP_DEFLATED)
#     for filename in files_to_zip:
#         arcname = filename.rsplit('/', 1)[-1].lower()
#         zip_file.write(filename, arcname)
#     zip_file.close()
    
#     return blast_file_count

def create_blast_fastas(current_user, fasta_file):
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
    filename = re.search('(.*/users/.*)/uploads/.*.\w*', fasta_file)
    genome = SeqIO.read(fasta_file, "fasta").seq
    output = ""
    blast_file_count = 1  # keep track of num blast files created
    out_file = f"{str(filename.group(1))}/{current_user}_blast_{blast_file_count}.fasta"
    files_to_zip = [out_file]
    starts, stops = get_start_stop('+', genome)
    for i in range(len(starts)):
        output += f">+, {starts[i]}-{stops[i]}\n"
        output += f"{Seq.translate(sequence=helper.get_sequence(genome, '+', starts[i]-1, stops[i]), table=11)}\n"
        if getsizeof(output) > 30000:  # only file size to reach 30 kb, else you get a CPU limit from NCBI blast
            with open(out_file, "w") as f:
                f.write(output)
            output = ""
            blast_file_count += 1
            out_file = f"{str(filename.group(1))}/{current_user}_blast_{blast_file_count}.fasta"
            files_to_zip.append(out_file)
    starts, stops = get_start_stop('-', genome)
    for i in range(len(starts)):
        output += f">-, {starts[i]}-{stops[i]}\n"
        start = len(genome) - stops[i]
        stop = len(genome) - starts[i] + 1
        output += f"{Seq.translate(sequence=helper.get_sequence(genome, '-', start, stop), table=11)}\n"
        if getsizeof(output) > 30000:  # only file size to reach 30 kb, else you get a CPU limit from NCBI blast
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
    bacteria_start_codons = ["ATG", "GTG", "TTG"]
    start_options = []
    maximum += 3
    gene = ""
    if (strand != "-"):
        gene = genome[minimum:maximum]
    else:
        print("complimentary")
        gene = genome.reverse_complement()[minimum:maximum]
    for index in range(0, len(gene)):
        codon = gene[index:index + 3]
        if (index == 3 or index == 217 or index == 218 or index == 219):
            print(codon)
        if codon in bacteria_start_codons:
            start_options.append(minimum + index + 1)
    return start_options

def get_starts_stops(cds_id, genome):
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
    # gdata_df = pd.read_csv(genemark_gdata_file, sep='\t', skiprows=16)
    # gdata_df.columns = ['Base', '1', '2', '3', '4', '5', '6']
    # gdata_df = gdata_df.set_index('Base')

    dnamaster_cds = DNAMaster.query.filter_by(id=cds_id).first()
    genemark_cds = GeneMark.query.filter_by(stop=dnamaster_cds.stop).first()

    minimum = -300 # 100 base pairs lower than the start position
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
        if (stop_options[-1] - start < 30 and start != dnamaster_cds.start):
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

    if dnamaster_cds.strand == '+':
        return start_options, stop_options
    else:
        return stop_options, start_options

def get_start_stop(strand, genome):
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
    # gdata_df = pd.read_csv(genemark_gdata_file, sep='\t', skiprows=16)
    # gdata_df.columns = ['Base', '1', '2', '3', '4', '5', '6']
    # gdata_df = gdata_df.set_index('Base')

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

    if len(start_options) == 0:
        start_options.append(dnamaster_cds.start)
        stop_options.append(dnamaster_cds.stop)

    if strand == '+':
        return start_options, stop_options
    else:
        return stop_options, start_options