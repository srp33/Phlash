"""Contains functions that are used in multiple pages.

Gets the path of a file.
Finds out if a file is allowed.
Gets the sequence of a genome.
Deletes the Blast zip file.

Attributes:
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
    BLAST_EXTENSIONS:
        The allowed blast extensions.
"""
import os
from Bio import SeqIO, Seq, SeqFeature
from builtins import FileNotFoundError

ROOT = os.path.dirname(os.path.abspath(__file__))
FASTA_EXTENSIONS = set(['.fasta', '.fna', '.fa'])
GDATA_EXTENSIONS = set(['.gdata'])

def allowed_file(filename, allowed_extensions):
    """Checks if file extension is acceptable.

    Args:
        filename:
            The name of the file.
        allowed_extensions:
            The extensions that are allowed.

    """
    return '.' in filename and os.path.splitext(filename)[1].lower() in allowed_extensions

def get_file_path(preference, upload_directory):
    """Gets path of required file.

    Args:
        preference:
            The type of file being searched for.
        upload_directory:
            The directory containing all of the uploaded files.

    Returns:
        The file path.
    """
    blast_files = []
    for filename in os.listdir(upload_directory):
        file_ext = os.path.splitext(filename)[1].lower()
        if preference == "fasta":
            if file_ext in FASTA_EXTENSIONS:
                return os.path.join(upload_directory, filename)
        elif preference == "gdata":
            if file_ext in GDATA_EXTENSIONS:
                return os.path.join(upload_directory, filename)
        else:
            return("Couldn't find file.")
    return blast_files

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

def delete_blast_zip(UPLOAD_FOLDER):
    """Deletes the Blast zip so that any updates that are made are reflected in the blast files.

    Args:
        UPLOAD_FOLDER:
            The directory containing all of the uploaded files.

    Returns:
        A success message.
    """
    try:
        USER_FOLDER = UPLOAD_FOLDER[:-8]
        for file in os.listdir(USER_FOLDER):
            if file.endswith(".zip") or file.endswith(".fasta"):
                os.remove(os.path.join(USER_FOLDER, file))
    except FileNotFoundError:
        return "Blast zip file has not been created yet."
    return "success"

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
        print(coding_potential[y_key][left_index])
        print(left_index)
        print(coding_potential[y_key][right_index])
        print(right_index)
        status = "Fail"
    return frame, status