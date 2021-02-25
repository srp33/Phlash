"""Contains the functions for the Genbank page.

Modifies and returns the GenBank file.

Attributes:
    response_object:
        The dictionary that is returned by the main functions.
"""
import helper
from Bio import SeqIO, Seq, SeqFeature
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from models import *
import re
import json
import pandas as pd
from datetime import datetime
import os

response_object = {}

def get_genbank(UPLOAD_FOLDER, current_user, payload):
    """Calls function to create the genbank file and returns it.

    Args:
        UPLOAD_FOLDER:
            The folder containing all of the uploaded files.
        current_user:
            The current user's ID.
        payload:
            The data sent from the front-end.
    
    Returns:
        The GenBank file.
    """
    gb_file = helper.get_file_path("genbank", UPLOAD_FOLDER)
    fasta_file = helper.get_file_path("fasta", UPLOAD_FOLDER)
    gb_file = create_genbank(fasta_file, UPLOAD_FOLDER, current_user, payload)
    f = open(gb_file, "r")
    return f.read()

def create_genbank(fasta_file, UPLOAD_FOLDER, current_user, payload):
    """Creates and returns the genbank file.
    Args:
        fasta_file:
            The file containing the phage's dna sequence.
        UPLOAD_FOLDER:
            The folder containing all of the uploaded files.
        current_user:
            The current user's ID.
        payload:
            The data sent from the front-end.
    
    Returns:
        The GenBank file.

    """
    headers = payload.get_json()
    print(headers)
    gb_file = os.path.join(UPLOAD_FOLDER, current_user + ".gb")
    genome = SeqIO.read(fasta_file, "fasta").seq
    genome = Seq(str(genome), IUPAC.unambiguous_dna)
    record = SeqRecord(genome, id='', name=headers["phageName"], description=headers["source"])
    record.annotations["AUTHORS"] = "Becker, L.W."
    record.annotations["Reference"] = "whole thing"

    qualifiers = {}
    qualifiers["organism"] = headers["source"]
    qualifiers["mol_type"] = headers["molType"]
    qualifiers["isolation_source"] = headers["isolationSource"]
    qualifiers["lab_host"] = headers["labHost"]
    qualifiers["country"] = headers["country"]
    qualifiers["identified_by"] = headers["identifiedBy"]
    qualifiers["note"] = headers["notes"]
    feature = SeqFeature(FeatureLocation(start=0, end=len(genome)), type='source', qualifiers=qualifiers)
    record.features.append(feature)

    idNumber = 0
    for cds in DNAMaster.query.order_by(DNAMaster.start).all():
        if (cds.function == "DELETED" or cds.status == "trnaDELETED"):
            continue
        idNumber += 1
        if cds.strand == '-':
            qualifiers = {}
            qualifiers["gene"] = str(idNumber)
            qualifiers["locus_tag"] = cds.id[0:-1] + str(idNumber)
            if headers["includeNotes"]:
                qualifiers["note"] = cds.notes
            feature = SeqFeature(FeatureLocation(start=cds.start - 1, end=cds.stop, strand=-1), id=cds.id[0:-1] + str(idNumber), type='gene', qualifiers=qualifiers)
            record.features.append(feature)
            if cds.status == "tRNA":
                qualifiers = {}
                qualifiers["gene"] = str(idNumber)
                qualifiers["locus_tag"] = headers["phageName"] + '_' + str(idNumber)
                qualifiers["note"] = cds.function
                feature = SeqFeature(FeatureLocation(start=cds.start - 1, end=cds.stop, strand=-1), id=cds.id[0:-1] + str(idNumber), type='tRNA', qualifiers=qualifiers)
                record.features.append(feature)
            else:
                qualifiers = {}
                qualifiers["gene"] = str(idNumber)
                qualifiers["locus_tag"] = headers["phageName"] + '_' + str(idNumber)
                qualifiers["codon_start"] = [1]
                qualifiers["transl_table"] = [11]
                pattern = re.compile("@*(.*)##(.*)")
                matches = pattern.search(cds.function)
                if matches:
                    qualifiers["product"] = matches.group(1)
                    qualifiers["protein_id"] = matches.group(2)
                    print(qualifiers)
                else:
                    qualifiers["product"] = "Hypothetical Protein"
                    qualifiers["protein_id"] = "unknown:" + qualifiers["locus_tag"]
                start = len(genome) - cds.stop
                stop = len(genome) - cds.start + 1
                qualifiers["translation"] = Seq.translate(helper.get_sequence(genome, cds.strand, start, stop), table=11)[0:-1]
                feature = SeqFeature(FeatureLocation(start=cds.start - 1, end=cds.stop, strand=-1), id=cds.id[0:-1] + str(idNumber), type='CDS', qualifiers=qualifiers)
                record.features.append(feature)
        else:
            qualifiers = {}
            qualifiers["gene"] = str(idNumber)
            qualifiers["locus_tag"] = headers["phageName"] + '_' + str(idNumber)
            if headers["includeNotes"]:
                qualifiers["note"] = cds.notes
            feature = SeqFeature(FeatureLocation(start=cds.start - 1, end=cds.stop), id=cds.id[0:-1] + str(idNumber), type='gene', qualifiers=qualifiers)
            record.features.append(feature)
            if cds.status == "tRNA":
                qualifiers = {}
                qualifiers["gene"] = str(idNumber)
                qualifiers["locus_tag"] = headers["phageName"] + '_' + str(idNumber)
                qualifiers["note"] = cds.function
                feature = SeqFeature(FeatureLocation(start=cds.start - 1, end=cds.stop), id=cds.id[0:-1] + str(idNumber), type='tRNA', qualifiers=qualifiers)
                record.features.append(feature)
            else:
                qualifiers = {}
                qualifiers["gene"] = str(idNumber)
                qualifiers["locus_tag"] = headers["phageName"] + '_' + str(idNumber)
                qualifiers["codon_start"] = [1]
                qualifiers["transl_table"] = [11]
                pattern = re.compile("(.*)##(.*)")
                matches = pattern.search(cds.function)
                if matches:
                    qualifiers["product"] = matches.group(1)
                    qualifiers["protein_id"] = matches.group(2)
                qualifiers["translation"] = Seq.translate(helper.get_sequence(genome, cds.strand, cds.start - 1, cds.stop), table=11)[0:-1]
                feature = SeqFeature(FeatureLocation(start=cds.start - 1, end=cds.stop), id=cds.id[0:-1] + str(idNumber), type='CDS', qualifiers=qualifiers)
                record.features.append(feature)
    with open(gb_file, 'w') as genbank:
        SeqIO.write(record, genbank, 'genbank')
    new_lines = []
    with open (gb_file, 'r') as genbank:
        lines = genbank.readlines()
        for index, line in enumerate(lines):
            if index is 0:
                new_lines.append(line[0:-28] + "     linear       " + datetime.now().strftime('%d-%b-%Y').upper() + '\n')
            elif index is 5 or index is 6:
                if headers["source"] != "":
                    new_lines.append(line[0:-2] + headers["source"] + '\n')
            elif index is 7:
                if headers["organism"] != "":
                    new_lines.append(line[0:-2] + headers["organism"] + '\n')
                new_lines.append("REFERENCE   1  (bases 1 to " + str(len(genome)) + ")\n")
                long_line = "  AUTHORS   " + headers["authors"]
                while len(long_line) > 81:
                    new_lines.append(long_line[0:80] + '\n')
                    long_line = long_line[81:]
                new_lines.append(long_line + '\n')
                long_line = "  TITLE     " + headers["title"]
                while len(long_line) > 81:
                    new_lines.append(long_line[0:80] + '\n')
                    long_line = long_line[81:]
                new_lines.append(long_line + '\n')
                long_line = "  JOURNAL   " + headers["journal"]
                while len(long_line) > 81:
                    new_lines.append(long_line[0:80] + '\n')
                    long_line = long_line[81:]
                new_lines.append(long_line + '\n')
            else:
                new_lines.append(line)
    print(new_lines[0])
    with open (gb_file, 'w') as genbank:
        genbank.writelines(new_lines)
    return gb_file