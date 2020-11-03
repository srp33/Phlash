"""
Contains methods that are used in multiple pages.
"""
import os
from Bio import SeqIO, Seq, SeqFeature
from builtins import FileNotFoundError

ROOT = os.path.dirname(os.path.abspath(__file__))
FASTA_EXTENSIONS = set(['.fasta', '.fna'])
GENBANK_EXTENSIONS = set(['.gb', '.gbk', '.gbf'])
GDATA_EXTENSIONS = set(['.gdata'])
LDATA_EXTENSIONS = set(['.ldata'])
BLAST_EXTENSIONS = set(['.json'])

def allowed_file(filename, allowed_extensions):
    """
    Checks if file extension is acceptable.
    """
    return '.' in filename and os.path.splitext(filename)[1].lower() in allowed_extensions

def get_file_path(preference, upload_directory):
    """
    Gets path of required file.
    """
    blast_files = []
    for filename in os.listdir(upload_directory):
        file_ext = os.path.splitext(filename)[1].lower()
        if preference == "fasta":
            if file_ext in FASTA_EXTENSIONS:
                return os.path.join(upload_directory, filename)
        elif preference == "genbank":
            if file_ext in GENBANK_EXTENSIONS:
                return os.path.join(upload_directory, filename)
        elif preference == "ldata":
            if file_ext in LDATA_EXTENSIONS:
                return os.path.join(upload_directory, filename)
        elif preference == "gdata":
            if file_ext in GDATA_EXTENSIONS:
                return os.path.join(upload_directory, filename)
        elif preference == "blast":
            if file_ext in BLAST_EXTENSIONS:
                blast_files.append(os.path.join(upload_directory, filename))
        else:
            return("Couldn't find file.")
    return blast_files

def get_sequence(genome, strand, start, stop):
    '''
    Gets sequence (ranging from start to stop) of the direct or complementary strand of genome.
    '''
    if strand == '-':
        return genome.reverse_complement()[start : stop]
    else:
        return genome[start : stop]

def delete_blast_zip(UPLOAD_FOLDER):
    '''
    Deletes the Blast zip so that any updates that are made are reflected in the blast files.
    '''
    try:
        USER_FOLDER = UPLOAD_FOLDER[:-8]
        for file in os.listdir(USER_FOLDER):
            if file.endswith(".zip") or file.endswith(".fasta"):
                os.remove(os.path.join(USER_FOLDER, file))
    except FileNotFoundError:
        return "Blast zip file has not been created yet."
    return "success"