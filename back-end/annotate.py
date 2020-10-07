"""
Contains functions that are called in main.py, the Flask application.
Functions parse through user's uploaded files, does comparison analysis 
between files, creates files, etc.
"""
from Bio import SeqIO, Seq, SeqFeature
from collections import OrderedDict
from models import *
import json
import os
import pandas as pd
import re
import subprocess
import zipfile
from sys import getsizeof


# HELPER FUNCTIONS ---------------------------------------------------------
def get_keys_by_value(dict, value_to_find):
    """
    Parses through GeneMark ldata. Used by `parse_genemark_ldata`.
    @param dict:
    @param value_to_find: 
    """
    keys = list()
    items = dict.items()
    for item in items:
        if value_to_find in item[1]:
            keys.append(item[0])
    return keys


def get_sequence(genome, strand, start, stop):
    """
    Gets sequence (ranging from start to stop) of the direct or 
    complementary strand of genome.
    * Remember that python indexing starts at 0, so make sure the start 
        value accomadates that when passed in to this function. 
    """
    if strand == '-':
        return genome.reverse_complement()[start : stop]
    else:
        return genome[start : stop]


def zip_files(files, zip_handle):
    for file in files:
        zip_handle.write(file)


# MAIN FUNCTIONS ---------------------------------------------------------
def run_genemark(fasta_file_path):
    """
    Invokes the GeneMarkS utilities and returns .gdata and .ldata files.
    """
    result = subprocess.run(["/genemark_suite_linux_64/gmsuite/gc", fasta_file_path], stdout=subprocess.PIPE)
    gc_percent = result.stdout.decode("utf-8").split(" ")[3]
    gc_percent = "{:d}".format(round(float(gc_percent)))

    subprocess.run(["/genemark_suite_linux_64/gmsuite/gm", "-D", "-g", "0", "-s", "1", "-m", "/genemark_suite_linux_64/gmsuite/heuristic_mat/heu_11_{}.mat".format(gc_percent), "-v", fasta_file_path])

    gdata_file_path = "{}.gdata".format(fasta_file_path)
    ldata_file_path = "{}.ldata".format(fasta_file_path)


def parse_dnamaster_genbank(genbank_file):
    """
    Parses through DNA Master GenBank file to gather DNA Master's CDS calls. 
    Adds each CDS to DNA Master table in user's database. 
    """
    print("in parse_dnamaster_genbank")
    DNAMaster.query.delete()
    with open(genbank_file, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            num = 1
            for feature in record.features:
                if feature.type == "CDS":
                    if "locus_tag" in feature.qualifiers:
                        id = feature.qualifiers["locus_tag"][0]
                    elif "protein_id" in feature.qualifiers:
                        id = feature.qualifiers["protein_id"][0]
                    else:
                        id = f"cds_{num}"
                        num += 1

                    # FIXME: do something for compound locations, e.g. join(1..218,166710..167034)
                    # if isinstance(feature.location, SeqFeature.CompoundLocation):
                    #     print(f"{feature.location} is a compoundlocation")
                    
                    strand = "+" if feature.location.strand == 1 else "-"
                    cds = DNAMaster(id = id,
                                    start = feature.location.start.position + 1,
                                    stop = feature.location.end.position,
                                    strand = strand,
                                    function = "None selected",
                                    status = "None")
                    
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
            dnamaster_cds.status = "Undetermined"
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
    print(start_options)
    return start_options


def create_blast_fasta(current_user, fasta_file, gdata_file):
    """
    Creates fasta file(s) for BLAST input.
    If more than one file is created, then each file should have
    100 sequences max (200 lines).
    @return blast_file_count: number of blast fasta files created.
    """
    filename = re.search('(.*/users/.*)/uploads/.*.\w*', fasta_file)
    genome = SeqIO.read(fasta_file, "fasta").seq
    output = ""
    blast_file_count = 1  # keep track of num blast files created
    out_file = f"{str(filename.group(1))}/{current_user}_blast_{blast_file_count}.fasta"
    files_to_zip = [out_file]
    for cds in db.session.query(DNAMaster).order_by(DNAMaster.start):
        starts = start_options(cds.id, genome, gdata_file)
        cds.start_options = ", ".join([str(start) for start in starts])
        db.session.commit()
        for i in range(len(starts[:5])):
            if starts[i] == cds.start:
                output += f">{cds.id}, {starts[i]}-{cds.stop}\n"
                output += f"{Seq.translate(sequence=get_sequence(genome, cds.strand, starts[i]-1, cds.stop), table=11)}\n"
            else:
                output += f">{cds.id}_{i}, {starts[i]}-{cds.stop}\n"
                output += f"{Seq.translate(sequence=get_sequence(genome, cds.strand, starts[i]-1, cds.stop), table=11)}\n"
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


def parse_blast_results(blast_files, cds_id, e_value_thresh):
    blast_results = {}
    for blast_file in blast_files:
        newLines = []
        with open(blast_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line != "CREATE_VIEW\n":
                    newLines.append(line)
        with open(blast_file, 'w') as f:
            f.writelines(newLines)    
            
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
