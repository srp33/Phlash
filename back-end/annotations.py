"""Contains the functions for the Annotations page.

Returns DNAMaster data.
Modifies and returns the GenBank file.

Attributes:
    response_object:
        The dictionary that is returned by the main functions.
"""
import helper
from Bio import SeqIO, Seq, SeqFeature
from models import *
import re
import json
import pandas as pd
from datetime import datetime

response_object = {}

# ------------------------------ MAIN FUNCTIONS ------------------------------
def get_dnamaster_data():
    """Queries and returns all CDS data created by DNAMaster.

    Returns:
        A dictionary containing the DNAMaster data.
    """
    setting = db.session.query(Settings).order_by(Settings.back_start_range).first()
    response_object['gap'] = setting.gap
    response_object['opposite_gap'] = setting.opposite_gap
    response_object['overlap'] = setting.overlap
    response_object['short'] = setting.short
    dnamaster = []
    for cds in db.session.query(DNAMaster).order_by(DNAMaster.start):
        dnamaster.append({'id': cds.id,
                            'start': cds.start,
                            'stop': cds.stop,
                            'strand': cds.strand,
                            'function': cds.function,
                            'status': cds.status,
                            'frame': cds.frame})
    response_object['dnamaster'] = dnamaster

    return response_object

def get_genbank(UPLOAD_FOLDER):
    """Modifies and returns the GenBank file.

    Args:
        UPLOAD_FOLDER:
            The folder containing all of the uploaded files.
    
    Returns:
        The GenBank file.
    """
    gb_file = helper.get_file_path("genbank", UPLOAD_FOLDER)
    fasta_file = helper.get_file_path("fasta", UPLOAD_FOLDER)
    out_file = modify_genbank(gb_file, fasta_file)
    f = open(out_file, "r")
    return f.read()

# ---------- GENBANK HELPER FUNCTIONS ----------
def modify_genbank(gb_file, fasta_file):
    """Creates modified GenBank file for input into Sequin.

    Args:
        gb_file:
            The GenBank file to be modified.
        fasta_file:
            The fasta file containing the DNA sequence.

    Returns:
        The modified GenBank file.
    """
    gb_filename = re.search(r'(.*/users/.*/uploads/.*).(\w*)', gb_file)
    out_file = str(gb_filename.group(1)) + '_modified.' + str(gb_filename.group(2))

    genome = SeqIO.read(fasta_file, "fasta").seq
    final_annotations = get_final_annotations(genome)
    final_features = []
    for record in SeqIO.parse(open(gb_file, "r"), "genbank"):
        for feature in record.features:
            if feature.type == "gene" or feature.type == "CDS" or feature.type == "tRNA":
                try:
                    locus_tag = feature.qualifiers["locus_tag"][0]
                    if locus_tag in final_annotations.keys():
                        new_start = final_annotations[locus_tag]["start"]
                        feature.location = SeqFeature.FeatureLocation(SeqFeature.ExactPosition(new_start - 1),
                                                                        SeqFeature.ExactPosition(feature.location.end.position),
                                                                        feature.location.strand)
                        if feature.type == "CDS" or feature.type == "tRNA":
                            feature.qualifiers["product"][0] = final_annotations[locus_tag]["function"]
                            if feature.qualifiers["product"][0] == "DELETED":
                                final_features.pop()
                                continue
                            feature.qualifiers["translation"][0] = final_annotations[locus_tag]["translation"]
                except:
                    continue
            else:
                continue
            final_features.append(feature)  # Append final features
        record.features = final_features
        with open(out_file, "w") as new_gb:
            SeqIO.write(record, new_gb, "genbank")
            
    return out_file

def get_final_annotations(genome):
    """Queries and returns the annotations made in DNAMaster database.

    Args:
        genome:
            The Phage genome.
    
    Returns:
        All of the annotations made for each CDS.
    """
    final_annotations = {}
    for cds in DNAMaster.query.order_by(DNAMaster.start).all():
        annotation = {}
        annotation['start'] = cds.start
        annotation['strand'] = cds.strand
        annotation['function'] = cds.function
        if (cds.strand is '+'):
            annotation['translation'] = Seq.translate(sequence=helper.get_sequence(genome, cds.strand, cds.start - 1, cds.stop), table=11)
        else:
            start = len(genome) - cds.stop
            stop = len(genome) - cds.start + 1
            annotation['translation'] = Seq.translate(sequence=helper.get_sequence(genome, cds.strand, start, stop), table=11)
        final_annotations[cds.id] = annotation

    return final_annotations

def parse_blast(UPLOAD_FOLDER):
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
    print(datetime.now())
    blast_files = helper.get_file_path("blast", UPLOAD_FOLDER)
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
                blasts = json.load(f)["BlastOutput2"]
                for blast in blasts:
                    search = blast["report"]["results"]["search"]
                    title = re.search(
                        "(.), (\d+)-(\d+)", search["query_title"])
                    if title:
                        curr_strand = title.group(1)
                        curr_start = title.group(2)
                        curr_stop = title.group(3)
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
                        blast_result = Blast_Results(id = counter,
                                                    start = curr_start,
                                                    stop = curr_stop,
                                                    strand = curr_strand,
                                                    results = str(results))
                        db.session.add(blast_result)
                        db.session.commit()
    except:
        print("leave me alone!")

    print(datetime.now())

def add_cds(request, UPLOAD_FOLDER):
    new_cds_data = request.get_json()
    force = new_cds_data.get('force')
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
    frame, status = helper.get_frame_and_status(new_cds_data.get('start'), new_cds_data.get('stop'), new_cds_data.get('strand'), coding_potential)
    cds = DNAMaster(id = new_cds_data.get('id'),
                    start = new_cds_data.get('start'),
                    stop = new_cds_data.get('stop'),
                    strand = new_cds_data.get('strand'),
                    function = "None selected",
                    status = status,
                    frame = frame)

    exists = DNAMaster.query.filter_by(id=new_cds_data.get('id')).first()
    orf = Blast_Results.query.filter_by(start=new_cds_data.get('start'), stop=new_cds_data.get('stop'), strand=new_cds_data.get('strand')).first()
    if force:
        if exists:
            exists.start = new_cds_data.get('start')
            exists.stop = new_cds_data.get('stop')
            exists.strand = new_cds_data.get('strand')
            db.session.commit()
        else:
            db.session.add(cds)
            db.session.commit()
            add_genbank_cds(UPLOAD_FOLDER, cds)
        response_object['message'] = "Added succesfully."
        return response_object
    if exists:
        response_object['message'] = "ID already exists."
    elif not orf:
        response_object['message'] = "Not orf."
    else:
        db.session.add(cds)
        db.session.commit()
        response_object['message'] = "Added succesfully."
        add_genbank_cds(UPLOAD_FOLDER, cds)
    return response_object

def add_genbank_cds(UPLOAD_FOLDER, cds):
    contents = []
    with open(helper.get_file_path("genbank", UPLOAD_FOLDER), "r") as genbank:
        contents = genbank.readlines()
    index = 0
    for i, line in enumerate(contents):
        if line.startswith("ORIGIN"):
            index = i
            break
    new_contents = []
    prev_gene = int(contents[index - 9][28:31])
    fasta_file = helper.get_file_path("fasta", UPLOAD_FOLDER)
    genome = SeqIO.read(fasta_file, "fasta").seq
    translation = ""
    if cds.strand == '-':
        start = len(genome) - cds.stop
        stop = len(genome) - cds.start + 1
        translation = Seq.translate(sequence=helper.get_sequence(genome, cds.strand, start, stop), table=11)
        new_contents.append("     gene            complement(" + str(cds.start) + ".." + str(cds.stop) + ")\n")
        new_contents.append('                     /gene="' + str((prev_gene + 1)) + '"\n')
        new_contents.append('                     /locus_tag="' + cds.id + '"\n')
        new_contents.append("     CDS             complement(" + str(cds.start) + ".." + str(cds.stop) + ")\n")
    else:
        translation = Seq.translate(sequence=helper.get_sequence(genome, cds.strand, start - 1, stop), table=11)
        new_contents.append("     gene            " + str(cds.start) + ".." + str(cds.stop) + "\n")
        new_contents.append('                     /gene="' + str((prev_gene + 1)) + '"\n')
        new_contents.append('                     /locus_tag="' + cds.id + '"\n')
        new_contents.append("     CDS             " + str(cds.start) + ".." + str(cds.stop) + "\n")
    new_contents.append('                     /gene="' + str((prev_gene + 1)) + '"\n')
    new_contents.append('                     /locus_tag="' + cds.id + '"\n')
    new_contents.append('                     /codon_start=1\n')
    new_contents.append('                     /product="Hypothetical Protein"\n')
    new_contents.append('                     /protein_id="unknown:' + cds.id + '"\n')
    translation = translation[0:-1]
    if len(translation) < 44:
        new_contents.append('                     /translation="' + str(translation) + '"\n')
    else:
        new_contents.append('                     /translation="' + str(translation[0:44]) + '\n')
        translation = translation[44:]
        empty = False
        while not empty:
            if len(translation) < 44:
                new_contents.append('                     ' + str(translation) + '"\n')
                empty = True
            else:
                new_contents.append('                     ' + str(translation[0:44]) + '\n')
                translation = translation[44:]
    print(new_contents)
    for i in range(len(new_contents)):
        contents.insert(i + index, new_contents[i])
    with open(helper.get_file_path("genbank", UPLOAD_FOLDER), "w") as genbank:
        genbank.writelines(contents)