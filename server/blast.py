from Bio import SeqIO, Seq
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import *
import re

# start = 50 - 1
# stop = 364
# record = SeqIO.read("uploads/fern.fasta", "fasta")
# table = 11  # Bacterial code
# genome = record.seq

# sequence = Seq.translate(sequence=genome[start: stop], table=11)
# print(sequence)

def run_blast(fasta_file, start, stop, e_value_thresh):
    start = start - 1
    record = SeqIO.read(fasta_file, "fasta")
    genome = record.seq
    sequence = Seq.translate(sequence=genome[start: stop], table=11)
    print("In `run_blast`")
    result_handle = NCBIWWW.qblast("blastp", "nr", sequence[:-1])
    print("Finished BLAST query")

    with open("blast.xml", "w") as out_handle:
        print("Writing BLAST results to file...")
        blast_results = result_handle.read()
        out_handle.write(blast_results)
    result_handle.close()

    blast = parse_blast(e_value_thresh)
    return blast

# def run_blast(e_value_thresh):
#     print("In `run_blast`")
#     sequences = open("test.fasta").read()
#     result_handle = NCBIWWW.qblast("blastp", "nr", sequences)
#     print("Finished BLAST query")

#     with open("blast.xml", "w") as out_handle:
#         print("Writing BLAST results to file...")
#         blast_results = result_handle.read()
#         out_handle.write(blast_results)
#     result_handle.close()

# run_blast(1e-7)

def parse_blast(e_value_thresh):
    blast = []
    for record in NCBIXML.parse(open("blast.xml")): 
        if record.alignments:
            # blast.append({'title': record.query[:100]})
            for align in record.alignments: 
                alignment = {}
                for hsp in align.hsps: 
                    if hsp.expect <= e_value_thresh: 
                        hit_definition = re.search('([A-Za-z\s]+\[[\.-_A-Za-z\s]+\]).+', align.hit_def)
                        if hit_definition:
                            alignment['hit_def'] = hit_definition.group(1)
                        else:
                            alignment['hit_def'] = align.hit_def
                        # alignment['title'] = align.title
                        alignment['e_value'] = float(hsp.expect)
                        # alignment['query'] = hsp.query
                        alignment['query_start'] = hsp.query_start
                        alignment['query_end'] = hsp.query_end
                        # alignment['match'] = hsp.match
                        # alignment['sbjct'] = hsp.sbjct
                        alignment['sbjct_start'] = hsp.sbjct_start
                        alignment['sbjct_end'] = hsp.sbjct_end
                        blast.append(alignment)
    return blast
