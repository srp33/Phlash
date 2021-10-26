import subprocess
import pandas as pd
import re
import os
from models import *
import models
from datetime import datetime
import json
from Bio import SeqIO, Seq, SeqFeature
from sys import getsizeof
import zipfile

def run_genemark(UPLOAD_FOLDER, phage_id):
    fasta_file_path = os.path.join(UPLOAD_FOLDER, phage_id + ".fasta")
    result = subprocess.run(["/genemark_suite_linux_64/gmsuite/gc", fasta_file_path], stdout=subprocess.PIPE)
    gc_percent = result.stdout.decode("utf-8").split(" ")[3]
    gc_percent = "{:d}".format(round(float(gc_percent)))
    subprocess.run(["/genemark_suite_linux_64/gmsuite/gm", "-D", "-g", "0", "-s", "1", "-m", "/genemark_suite_linux_64/gmsuite/heuristic_mat/heu_11_{}.mat".format(gc_percent), "-v", fasta_file_path])

def run_glimmer(UPLOAD_FOLDER, phage_id):
    fasta_file_path = os.path.join(UPLOAD_FOLDER, phage_id + ".fasta")
    glimmer_file = os.path.join(UPLOAD_FOLDER, phage_id + "_glimmer")
    subprocess.run(["/opt/glimmer3.02/scripts/g3-from-scratch.csh", fasta_file_path, glimmer_file])

def run_aragorn(UPLOAD_FOLDER, phage_id):
    fasta_file_path = os.path.join(UPLOAD_FOLDER, phage_id + ".fasta")
    aragorn_file = os.path.join(UPLOAD_FOLDER, phage_id + "_aragorn.txt")
    subprocess.run(["aragorn", fasta_file_path, aragorn_file])

def run_phanotate(UPLOAD_FOLDER, phage_id):
    fasta_file_path = os.path.join(UPLOAD_FOLDER, phage_id + ".fasta")
    phanotate_file = os.path.join(UPLOAD_FOLDER, phage_id + "_phanotate.txt")
    subprocess.run(["/usr/local/lib/python3.7/site-packages/phanotate.py", fasta_file_path, "-o" + phanotate_file], timeout=300)

def run_prodigal(UPLOAD_FOLDER, phage_id):
    fasta_file_path = os.path.join(UPLOAD_FOLDER, phage_id + ".fasta")
    prodigal_file = os.path.join(UPLOAD_FOLDER, phage_id + "_prodigal.txt")
    subprocess.run(["/usr/local/bin/prodigal", "-i", fasta_file_path, "-o", prodigal_file, "-f", "sco"])

def auto_annotate(UPLOAD_FOLDER, phage_id, session):
    """Auto annotates the entire genome.

    Args:
        UPLOAD_FOLDER:
            The path to the folder containing all files associated with the current user.
        phage_id:
            The ID of the current user.
    """
    fasta_file_path = os.path.join(UPLOAD_FOLDER, phage_id + ".fasta")
    gdata_file_path = "{}.gdata".format(fasta_file_path)
    ldata_file_path = "{}.ldata".format(fasta_file_path)
    phanotate_file = os.path.join(UPLOAD_FOLDER, phage_id + "_phanotate.txt")
    prodigal_file = os.path.join(UPLOAD_FOLDER, phage_id + "_prodigal.txt")

    coding_potential = {}
    gdata_df = pd.read_csv(gdata_file_path, sep='\t', skiprows=16)
    gdata_df.columns = ['Base', '1', '2', '3', '4', '5', '6']
    coding_potential['x_data'] = gdata_df["Base"].to_list()
    coding_potential['y_data_1'] = gdata_df["1"].to_list()
    coding_potential['y_data_2'] = gdata_df["2"].to_list()
    coding_potential['y_data_3'] = gdata_df["3"].to_list()
    coding_potential['y_data_4'] = gdata_df["4"].to_list()
    coding_potential['y_data_5'] = gdata_df["5"].to_list()
    coding_potential['y_data_6'] = gdata_df["6"].to_list()

    setting = session.query(Settings).filter_by(phage_id=phage_id).order_by(Settings.back_left_range).first()

    # Add glimmer calls.
    if setting.glimmer:
        id_index = 0
        glimmer_calls = ""
        with open(os.path.join(UPLOAD_FOLDER, phage_id + "_glimmer.predict"), 'r') as glimmer:
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
                frame, status = get_frame_and_status(left, right, strand, coding_potential)
                cds = Annotations(phage_id = phage_id,
                                id = 'Glimmer_' + str(id_index),
                                left = left,
                                right = right,
                                strand = strand,
                                function = "None selected",
                                status = status,
                                frame = frame)
                session.add(cds)
                glimmer_calls += (str(left) + '-' + str(right) + ' ' + strand + ',')
        calls = Gene_Calls(phage_id = phage_id,
                            id = 'Glimmer',
                            calls = glimmer_calls)
        session.add(calls)
        session.commit()

    # Add aragorn calls.
    if setting.aragorn:
        with open(os.path.join(UPLOAD_FOLDER, phage_id + "_aragorn.txt"), 'r') as aragorn:
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
                    cds = Annotations(phage_id = phage_id,
                                    id = "Aragorn_" + str(id_index),
                                    left = int(left),
                                    right = int(right),
                                    strand = strand,
                                    function = function,
                                    status = "tRNA",
                                    frame = 0)
                    session.add(cds)
        session.commit()

    # Add genemark calls.
    if setting.genemark:
        genemark_calls = ""
        with open(ldata_file_path) as genemark:
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
                            cds = session.query(Annotations).filter_by(phage_id=phage_id).filter_by(left=left).first()
                            if cds == None or cds.status == 'Fail':
                                frame, status = get_frame_and_status(left, right, '-', coding_potential)
                                if cds == None or (cds.status == 'Fail' and status == 'Pass'):
                                    if cds != None:
                                        session.delete(cds)
                                    id_index += 1
                                    cds = Annotations(phage_id = phage_id,
                                                    id = "Genemark_" + str(id_index),
                                                    left = int(left),
                                                    right = int(right),
                                                    strand = '-',
                                                    function = "None selected",
                                                    status = status,
                                                    frame = frame)
                                session.add(cds)
                                session.commit()
                        else:
                            genemark_calls += (str(left) + '-' + str(right) + ' ' + '+' + ',')
                            cds = session.query(Annotations).filter_by(phage_id=phage_id).filter_by(right=right).first()
                            if cds == None or cds.status == 'Fail':
                                frame, status = get_frame_and_status(left, right, '+', coding_potential)
                                if cds == None or (cds.status == 'Fail' and status == 'Pass'):
                                    if cds != None:
                                        session.delete(cds)
                                    id_index += 1
                                    cds = Annotations(phage_id = phage_id,
                                                    id = "Genemark_" + str(id_index),
                                                    left = int(left),
                                                    right = int(right),
                                                    strand = '+',
                                                    function = "None selected",
                                                    status = status,
                                                    frame = frame)
                                    session.add(cds)
                                    session.commit()
                        left = int(column[0])
                        right = int(column[1])
                        frame = int(column[2])
                        prob = new_prob
        calls = Gene_Calls(phage_id = phage_id,
                        id = 'GeneMark',
                        calls = genemark_calls)
        session.add(calls)
        session.commit()

    # Add prodigal calls.
    if setting.prodigal:
        prodigal_calls = ""
        id_index = 0
        with open(prodigal_file, 'r') as prodigal:
            for line in prodigal:
                if line[0] != '>':
                    continue
                id_index += 1
                column = line.strip().split('_')
                strand = str(column[3])
                left = int(column[2])
                right = int(column[1])
                cds = session.query(Annotations).filter_by(phage_id=phage_id).filter_by(left=left).first()
                if strand == '+':
                    left = int(column[1])
                    right = int(column[2])
                    cds = session.query(Annotations).filter_by(phage_id=phage_id).filter_by(right=right).first()
                prodigal_calls += (str(left) + '-' + str(right) + ' ' + strand + ',')
                if cds == None or cds.status == 'Fail':
                    frame, status = get_frame_and_status(left, right, strand, coding_potential)
                    if cds == None or (cds.status == 'Fail' and status == 'Pass'):
                        if cds != None:
                            session.delete(cds)
                        cds = Annotations(phage_id = phage_id,
                                        id = 'Prodigal_' + str(id_index),
                                        left = left,
                                        right = right,
                                        strand = strand,
                                        function = "None selected",
                                        status = status,
                                        frame = frame)
                        session.add(cds)
        calls = Gene_Calls(phage_id = phage_id,
                        id = 'Prodigal',
                        calls = prodigal_calls)
        session.add(calls)
        session.commit()

    # Add phanotate calls.
    if setting.phanotate:
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
                cds = session.query(Annotations).filter_by(phage_id=phage_id).filter_by(left=left).first()
                if strand == '+':
                    left = int(column[0])
                    right = int(column[1])
                    cds = session.query(Annotations).filter_by(phage_id=phage_id).filter_by(right=right).first()
                phanotate_calls += (str(left) + '-' + str(right) + ' ' + strand + ',')
                if cds == None or cds.status == 'Fail':
                    frame, status = get_frame_and_status(left, right, strand, coding_potential)
                    if cds == None or (cds.status == 'Fail' and status == 'Pass'):
                        if cds != None:
                            session.delete(cds)
                        cds = Annotations(phage_id = phage_id,
                                        id = 'Phanotate_' + str(id_index),
                                        left = left,
                                        right = right,
                                        strand = strand,
                                        function = "None selected",
                                        status = status,
                                        frame = frame)
                        session.add(cds)
        calls = Gene_Calls(phage_id = phage_id,
                        id = 'Phanotate',
                        calls = phanotate_calls)
        session.add(calls)
        session.commit()
    id_index = 1

    # Rename gene calls.
    user = session.query(Users).filter_by(id=phage_id).first()
    for cds in session.query(Annotations).filter_by(phage_id=phage_id).order_by(Annotations.left):
        cds.id = user.phage_id + '_' + str(id_index)
        if (cds.right < cds.left) or (session.query(Blast_Results).filter_by(phage_id=phage_id).filter_by(left=cds.left, right=cds.right, strand=cds.strand).first()):
            session.delete(cds)
        else:
            id_index += 1
    session.commit()

    # Remove files created for auto-annotation.
    for filename in os.listdir(UPLOAD_FOLDER):
        if not filename.endswith(".fasta") and not filename.endswith(".gdata") and not filename.endswith(".json"):
            os.remove(os.path.join(UPLOAD_FOLDER, filename))

    return("success")

def get_frame_and_status(left, right, strand, coding_potential):
    """Finds the frame that a given CDS is on and whether it covers the coding potential.

    Args:
        left:
            The left position for a given cds.
        right:
            The right postition for a given cds.
        strand:
            The strand of a given cds.
        coding_potential:
            The coding potential created by genemark.
    
    Returns:
        The frame and status of a given CDS.
    """
    left_index = 0
    find_base = left
    found = False
    while not found:
        if find_base < 1 or find_base > len(coding_potential['x_data']):
            left_index = 0
            break
        found = True
        try:
            left_index = coding_potential['x_data'].index(find_base)
        except:
            found = False
            find_base -= 1
    if left_index > 0:
        left_index -= 1
    right_index = 0
    find_base = right
    found = False
    while not found:
        if find_base < 1 or find_base > len(coding_potential['x_data']):
            right_index = len(coding_potential['x_data']) - 1
            break
        found = True
        try:
            right_index = coding_potential['x_data'].index(find_base)
        except:
            found = False
            find_base += 1
    if right_index < len(coding_potential['x_data']) - 1:
        left_index += 1
    frame = 0
    if strand == "-":
        frame = ((right + 2) % 3) + 4
    else:
        frame = ((left + 2) % 3) + 1

    y_key = 'y_data_' + str(frame)
    status = "Pass"
    if coding_potential[y_key][left_index] >= .5 or coding_potential[y_key][right_index] >= .5:
        status = "Fail"
    return frame, status

def parse_blast(phage_id, UPLOAD_FOLDER, session):
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
    session.query(Blast_Results).filter_by(phage_id=phage_id).delete()
    blast_files = []
    for filename in os.listdir(UPLOAD_FOLDER):
        if filename.endswith('.json'):
            blast_files.append(os.path.join(UPLOAD_FOLDER, filename))
    E_VALUE_THRESH = 1e-7
    counter = 0
    try:
        for blast_file in blast_files:
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
                            session.add(blast_result)
                            session.commit()
                        except:
                            print("An error occured while adding a blast result. Probably a duplicate.")
                            session.query(Blast_Results).filter_by(phage_id=phage_id).delete()
                            session.commit()
                            return("error")
    except:
        print("An error occured while parsing blast results.")
        session.query(Blast_Results).filter_by(phage_id=phage_id).delete()
        session.commit()
        return("error")
        
    for filename in os.listdir(UPLOAD_FOLDER):
        if filename.endswith('.json'):
            os.remove(os.path.join(UPLOAD_FOLDER, filename))
            
    return(message)

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
    fasta_file = os.path.join(UPLOAD_FOLDER, phage_id + ".fasta")
    filename = re.search('(.*/users/.*)/uploads/.*.\w*', fasta_file)
    genome = SeqIO.read(fasta_file, "fasta").seq
    output = ""
    blast_file_count = 1  # keep track of num blast files created
    out_file = f"{str(filename.group(1))}/{phage_id}_blast_{blast_file_count}.fasta"
    files_to_zip = [out_file]
    lefts, rights = get_starts_stops('+', genome)
    for i in range(len(lefts)):
        output += f">+, {lefts[i]}-{rights[i]}\n"
        output += f"{Seq.translate(sequence=get_sequence(genome, '+', lefts[i]-1, rights[i]), table=11)}\n"
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
        output += f"{Seq.translate(sequence=get_sequence(genome, '-', left, right), table=11)}\n"
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
    USER_FOLDER = UPLOAD_FOLDER[:-8]
    for filename in os.listdir(os.path.join(USER_FOLDER)):
        if filename.endswith(".fasta"):
            os.remove(os.path.join(USER_FOLDER, filename))

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

def get_sequence(genome, strand, left, right):
    """Gets sequence (ranging from left to right) of the direct or complementary strand of genome.

    Args:
        genome:
            The genome of the phage.
        strand:
            Direct or complimentary strand.
        left:
            The index of the left codon.
        right:
            The index of the right codon.
    
    Returns:
        The sequence.
    """
    if strand == '-':
        return genome.reverse_complement()[left : right]
    else:
        return genome[left : right]