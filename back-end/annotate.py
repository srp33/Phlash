"""
Contains functions that are called in app.py.
Functions parse through user's uploaded files, does comparison analysis between files, 
creates files, etc.
"""
from Bio import SeqIO, Seq, SeqFeature
from models import *
import os
import pandas as pd
import re
import json


# HELPER FUNCTIONS -----------------------------------------------------------------------------------------------------

def parse_location(location_stg):
    """Parses through GenBank's gene location string. Used by `parse_genbank`.

    @param location_stg: 
    """
    start = ""
    stop = ""
    frame = ""
    parsed = []
    location_stg = str(location_stg)
    for char in location_stg:
        if char is "[" or char is "]":
            parsed.append(char)
        elif char is ":":
            parsed.append(char)
        elif char is "(" or char is ")":
            parsed.append(char)
        elif char.isdigit() and ":" not in parsed:
            start = start + char
        elif char.isdigit() and ":" in parsed:
            stop = stop + char
        elif "(" in parsed:
            frame = char
    start = int(start) + 1
    stop = int(stop)
    return [start, stop, frame]


def get_keys_by_value(dict, value_to_find):
    """Parses through GeneMark ldata. Used by `parse_genemark_ldata`.

    @param dict:
    @param value_to_find: 
    """
    keys = list()
    items = dict.items()
    for item in items:
        if value_to_find in item[1]:
            keys.append(item[0])
    return keys


# Helper: add probabilities for each frame (used by avg_coding_potential_per_frame)
def add_frame_probs(df, base_positions):
    one = []
    two = []
    three = []
    four = []
    five = []
    six = []
    for base in base_positions:
        one.append(df.loc[base, '1'])
        two.append(df.loc[base, '2'])
        three.append(df.loc[base, '3'])
        four.append(df.loc[base, '4'])
        five.append(df.loc[base, '5'])
        six.append(df.loc[base, '6'])
    return [one, two, three, four, five, six]


# Helper: Find average (used by avg_coding_potential_per_frame)
def calculate_avg_prob(probabilities):
    total = 0
    for probability in probabilities:
        total += probability

    average = total / len(probabilities)
    return round(average, 4)


# Helper: Make dictionary {key: frame #, value: avg probability} (used by failed_genes)
def avg_coding_potential_per_frame(df, start, stop):
    indexes = []
    for index, row in df.iterrows():
        if start <= index <= stop:
            indexes.append(index)

    frames = add_frame_probs(df, indexes)
    avg_probs = {}
    current_frame = 1
    for frame in frames:
        avg_prob = calculate_avg_prob(frame)
        frame_label = "frame_{}".format(current_frame)
        avg_probs[frame_label] = avg_prob
        current_frame = current_frame + 1

    return avg_probs


# Get forward or reverse sequence
def get_sequence(genome, strand, start, stop):
    if strand == '-':
        return genome.reverse_complement()[start : stop]
    else:
        return genome[start : stop]

# MAIN FUNCTIONS ------------------------------------------------------------------------------------------------------

# Get all gene locations from GenBank file and add to sql database
def parse_genbank(genbank_file):
    with open(genbank_file, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            # id_number = 1
            for feature in record.features:
                if feature.type == "CDS":
                    # id = "dnamaster_" + str(id_number)
                    id = feature.qualifiers["locus_tag"][0]
                    gene_info = parse_location(feature.location)
                    cds = DNAMaster(id=id,
                                    start=gene_info[0],
                                    stop=gene_info[1],
                                    strand=gene_info[2],
                                    function="None",
                                    status="None")
                    exists = DNAMaster.query.filter_by(id=id).first()
                    if not exists:
                        db.session.add(cds)
                        db.session.commit()


# Parse through genemark ldata
def parse_genemark_ldata(gm_file):
    genemark_all_genes = dict()
    with open(gm_file, "r") as handle:
        for line in handle:
            if line == '\n':
                break
            if line[0] != '#':
                column = line.strip().split()
                start = int(column[0])
                stop = int(column[1])
                if stop not in genemark_all_genes:
                    genemark_all_genes[stop] = [start]
                else:
                    if start not in genemark_all_genes[stop]:
                        genemark_all_genes[stop].append(start)
                curr_keys = get_keys_by_value(genemark_all_genes, start)
                if len(curr_keys) > 1:
                    max_right = max(curr_keys)
                    for key in curr_keys:
                        if key != max_right:
                            del genemark_all_genes[key]

    genemark_cdss = dict()
    with open(gm_file, "r") as handle:
        for line in handle:
            if line == '\n':
                break
            if line[0] != '#':
                column = line.strip().split()
                curr_start = int(column[0])
                curr_stop = int(column[1])
                frame = int(column[2])
                if 1 <= frame <= 3:
                    curr_frame = "+"
                elif 4 <= frame <= 6:
                    curr_frame = "-"
                id_number = 1
                for stop in genemark_all_genes:
                    min_start = min(genemark_all_genes[stop])
                    if min_start == curr_start and stop == curr_stop:
                        if (min_start, stop) not in genemark_cdss:
                            id = "genemark_" + str(id_number)
                            genemark_cdss[(min_start, stop)] = curr_frame
                            cds = GeneMark(id=id,
                                           start=min_start,
                                           stop=stop,
                                           strand=curr_frame)
                            exists = GeneMark.query.filter_by(id=id).first()
                            if not exists:
                                db.session.add(cds)
                                db.session.commit()
                    id_number += 1


# Compare the gene calls between each tool
def compare():
    for cds in DNAMaster.query.all():
        dnamaster_cds = DNAMaster.query.filter_by(stop=cds.stop).first()
        genemark_cds = GeneMark.query.filter_by(stop=cds.stop).first()
        if dnamaster_cds and genemark_cds:
            if dnamaster_cds.start <= genemark_cds.start:
                dnamaster_cds.status = "Pass"
            else:
                dnamaster_cds.status = "Fail"
        elif not genemark_cds:
            dnamaster_cds.status = "Need more information"
        db.session.commit()


def start_options(cds_id, genome, genemark_gdata_file):
    """
    Finds alternative, possible start positions for a given CDS. 
    Calculates the average coding potential per frame for each start position.

    @return start_options: dictionary of (key) start positions and (value) their
        average coding potential per frame. 
    """
    gdata_df = pd.read_csv(genemark_gdata_file, sep='\t', skiprows=16)
    gdata_df.columns = ['Base', '1', '2', '3', '4', '5', '6']
    gdata_df = gdata_df.set_index('Base')

    bacteria_start_codons = ["ATG", "GTG", "TTG"]
    dnamaster_cds = DNAMaster.query.filter_by(id=cds_id).first()
    genemark_cds = GeneMark.query.filter_by(stop=dnamaster_cds.stop).first()

    original_start = dnamaster_cds.start - 1  # Subtract 1 because python indexing begins with 0.
    prev_start = original_start - 1  # Subtract 4 to account for 3 bases in a codon.

    start_options = []
    start_options.append(dnamaster_cds.start)  # Add original start position info

    min_start = original_start if genemark_cds is None else genemark_cds.start - 1
    num_nucleotides = min_start if min_start < 200 else 200   # Check 200 bp previous to original

    while prev_start >= min_start - num_nucleotides:
        prev_codon = get_sequence(genome, dnamaster_cds.strand, prev_start, prev_start + 3)
        if prev_codon in bacteria_start_codons:
            start_options.append(prev_start + 1)
        prev_start = prev_start - 1

    return start_options


# Create FASTA for input into BLAST
def create_fasta(fasta_file, genemark_gdata_file):
    genome = SeqIO.read(fasta_file, "fasta").seq
    output = ""
    
    for cds in db.session.query(DNAMaster).order_by(DNAMaster.start):
        starts = start_options(cds.id, genome, genemark_gdata_file)
        cds.start_options = ", ".join([str(start) for start in starts])
        db.session.commit()
        for i in range(len(starts[:3])):
            if starts[i] == cds.start:
                output += f">{cds.id}, {starts[i]}-{cds.stop}\n"
                output += f"{Seq.translate(sequence=get_sequence(genome, cds.strand, starts[i]-1, cds.stop), table=11)}\n"
            else:
                output += f">{cds.id}_{i}, {starts[i]}-{cds.stop}\n"
                output += f"{Seq.translate(sequence=get_sequence(genome, cds.strand, starts[i]-1, cds.stop), table=11)}\n"

    return output


def get_final_annotations(genome):
    final_annotations = {}
    for cds in DNAMaster.query.order_by(DNAMaster.start).all():
        annotation = {}
        annotation['start'] = cds.start
        annotation['strand'] = cds.strand
        annotation['function'] = cds.function
        annotation['translation'] = Seq.translate(sequence=get_sequence(genome, 
                                                                        cds.strand, 
                                                                        cds.start - 1,
                                                                        cds.stop), 
                                                                        table=11)
        final_annotations[cds.id] = annotation
    return final_annotations


# Create modified GenBank file for input into Sequin
def modify_genbank(gb_file, fasta_file):
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
                        feature.qualifiers["translation"][0] = final_annotations[locus_tag]["translation"]
                else:
                    continue
            final_features.append(feature)  # Append final features
        record.features = final_features
        with open(out_file, "w") as new_gb:
            SeqIO.write(record, new_gb, "genbank")
    
    return out_file


def parse_blast_multiple(blast_file, cds_id, e_value_thresh):
    blast_results = {}
    with open(blast_file) as f:
        blasts = json.load(f)["BlastOutput2"]
        for blast in blasts:
            search = blast["report"]["results"]["search"]
            title = re.search(
                "([A-Z]+_\d+_*\d*), (\d+)-\d+", search["query_title"])
            if title:
                curr_id = title.group(1)
                curr_start = int(title.group(2))
                cds_id_ = cds_id + "_"
                if cds_id == curr_id or cds_id_ in curr_id:
                    blast_results[curr_start] = []
                    hits = search["hits"]
                    for hit in hits:
                        hsps = hit["hsps"][0]
                        if hsps["evalue"] <= e_value_thresh:
                            alignment = {}
                            description = hit["description"][0]
                            alignment['accession'] = description["accession"]
                            alignment["title"] = description["title"]
                            alignment["sciname"] = description["sciname"]
                            alignment["evalue"] = '{:0.2e}'.format(
                                hsps["evalue"])
                            alignment["query_from"] = hsps["query_from"]
                            alignment["query_to"] = hsps["query_to"]
                            alignment["hit_from"] = hsps["hit_from"]
                            alignment["hit_to"] = hsps["hit_to"]
                            alignment["percent_identity"] = round(
                                hsps["identity"] / hsps["align_len"] * 100, 2)
                            blast_results[curr_start].append(alignment)

    return blast_results
