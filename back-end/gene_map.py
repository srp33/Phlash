"""Contains the functions for the Genome Map page.

Creates and returns the genome map.
"""
from dna_features_viewer import GraphicFeature, GraphicRecord
from Bio import SeqIO, Seq, SeqFeature
from werkzeug.utils import secure_filename
from contextlib import closing
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
import base64

# ------------------------------ MAIN FUNCTIONS ------------------------------
def get_map(UPLOAD_FOLDER):
    """Creates and returns a map of the genome.

    Args:
        UPLOAD_FOLDER:
            The folder containing all of the uploaded files.

    Returns:
        A dictionary containing an image of the genome map.

    """
    features = []
    for cds in db.session.query(DNAMaster).order_by(DNAMaster.start):
        if cds.function != 'DELETED':
            if cds.strand == '+':
                features.append(GraphicFeature(start=cds.start, end=cds.stop, strand=+1, color="#add8e6", label=cds.id))
            else:
                features.append(GraphicFeature(start=cds.start, end=cds.stop, strand=-1, color="#fed8b1", label=cds.id))
    fasta_file = helper.get_file_path("fasta", UPLOAD_FOLDER)
    genome = SeqIO.read(fasta_file, "fasta").seq
    sequence = str(genome)
    record = GraphicRecord(sequence_length=len(sequence), features=features)
    ax, _ = record.plot(figure_width=len(sequence)/1000)
    ax.figure.savefig(os.path.join(UPLOAD_FOLDER, 'sequence_and_translation.png'), bbox_inches='tight')
    image_byte_string = ""
    with open(os.path.join(UPLOAD_FOLDER, 'sequence_and_translation.png'), "rb") as image_file:
        image_byte_string = base64.b64encode(image_file.read())
    response_object = {}
    response_object['status'] = "success"
    response_object['image'] = str(image_byte_string)
    return response_object