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


def get_map(phage_id, UPLOAD_FOLDER):
    """Creates and returns a map of the genome.

    Args:
        UPLOAD_FOLDER:
            The folder containing all of the uploaded files.

    Returns:
        A dictionary containing an image of the genome map.

    """
    features = []
    for cds in db.session.query(Annotations).filter_by(phage_id=phage_id).order_by(Annotations.left):
        if cds.function != '@DELETED' and cds.status != 'trnaDELETED':
            if cds.strand == '+':
                if cds.status == "tRNA":
                    features.append(GraphicFeature(
                        start=cds.left, end=cds.right, strand=+1, color="#7570b3", label=cds.id))
                else:
                    features.append(GraphicFeature(
                        start=cds.left, end=cds.right, strand=+1, color="#1b9e77", label=cds.id))
            else:
                if cds.status == "tRNA":
                    features.append(GraphicFeature(
                        start=cds.left, end=cds.right, strand=-1, color="#7570b3", label=cds.id))
                else:
                    features.append(GraphicFeature(
                        start=cds.left, end=cds.right, strand=-1, color="#d95f02", label=cds.id))

    fasta_file = helper.get_file_path("fasta", UPLOAD_FOLDER)
    genome = SeqIO.read(fasta_file, "fasta").seq
    sequence = str(genome)
    record = GraphicRecord(sequence_length=len(sequence), features=features)
    ax, _ = record.plot(figure_width=len(sequence)/1000)
    ax.figure.savefig(os.path.join(
        UPLOAD_FOLDER, 'sequence_and_translation.png'), bbox_inches='tight')
    image_byte_string = ""

    with open(os.path.join(UPLOAD_FOLDER, 'sequence_and_translation.png'), "rb") as image_file:
        image_byte_string = base64.b64encode(image_file.read())

    response_object = {}
    response_object['status'] = "success"
    response_object['image'] = str(image_byte_string)
    return response_object
