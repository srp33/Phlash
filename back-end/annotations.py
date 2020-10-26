"""
Contains the methods for the Annotations page.
"""
import helper
from Bio import SeqIO, Seq, SeqFeature
from models import *
import re

response_object = {}

# ------------------------------ MAIN FUNCTIONS ------------------------------
def get_dnamaster_data():
    '''
    Queries and returns all CDS data created by DNAMaster.
    '''
    dnamaster = []
    for cds in db.session.query(DNAMaster).order_by(DNAMaster.start):
        dnamaster.append({'id': cds.id,
                            'start': cds.start,
                            'stop': cds.stop,
                            'strand': cds.strand,
                            'function': cds.function,
                            'status': cds.status})
    response_object['dnamaster'] = dnamaster

    return response_object

def get_genbank(UPLOAD_FOLDER):
    '''
    Modifies and returns the genbank file.
    '''
    gb_file = helper.get_file_path("genbank", UPLOAD_FOLDER)
    fasta_file = helper.get_file_path("fasta", UPLOAD_FOLDER)
    out_file = modify_genbank(gb_file, fasta_file)
    f = open(out_file, "r")
    return f.read()

# ---------- GENBANK HELPER FUNCTIONS ----------
def modify_genbank(gb_file, fasta_file):
    '''
    Creates modified GenBank file for input into Sequin
    '''
    gb_filename = re.search(r'(.*/users/.*/uploads/.*).(\w*)', gb_file)
    out_file = str(gb_filename.group(1)) + '_modified.' + str(gb_filename.group(2))

    genome = SeqIO.read(fasta_file, "fasta").seq
    final_annotations = get_final_annotations(genome)
    final_features = []
    for record in SeqIO.parse(open(gb_file, "r"), "genbank"):
        for feature in record.features:
            if feature.type == "gene" or feature.type == "CDS":
                locus_tag = feature.qualifiers["locus_tag"][0]
                if locus_tag in final_annotations.keys():
                    new_start = final_annotations[locus_tag]["start"]
                    feature.location = SeqFeature.FeatureLocation(SeqFeature.ExactPosition(new_start - 1),
                                                                    SeqFeature.ExactPosition(feature.location.end.position),
                                                                    feature.location.strand)
                    if feature.type == "CDS":
                        feature.qualifiers["product"][0] = final_annotations[locus_tag]["function"]
                        if feature.qualifiers["product"][0] == "DELETED":
                            final_features.pop()
                            continue
                        feature.qualifiers["translation"][0] = final_annotations[locus_tag]["translation"]
                else:
                    continue
            final_features.append(feature)  # Append final features
        record.features = final_features
        with open(out_file, "w") as new_gb:
            SeqIO.write(record, new_gb, "genbank")
            
    return out_file

def get_final_annotations(genome):
    '''
    Queries and returns the annotations made in DNAMaster database.
    '''
    final_annotations = {}
    for cds in DNAMaster.query.order_by(DNAMaster.start).all():
        annotation = {}
        annotation['start'] = cds.start
        annotation['strand'] = cds.strand
        annotation['function'] = cds.function
        annotation['translation'] = Seq.translate(sequence=helper.get_sequence(genome, 
                                                                        cds.strand, 
                                                                        cds.start - 1,
                                                                        cds.stop), 
                                                                        table=11)
        final_annotations[cds.id] = annotation

    return final_annotations