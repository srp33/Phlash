# from Bio import SeqIO, Seq
# from Bio.Blast import NCBIWWW, NCBIXML
# from Bio.Blast.Applications import *
import re
import json

def parse_blast(blast_file, cds_id, e_value_thresh):
   blast_results = []
   with open(blast_file) as f:
      blasts = json.load(f)["BlastOutput2"]
      for blast in blasts:
         search = blast["report"]["results"]["search"]
         title = re.search("(FERN_\d+_*\d*), \d+-\d+", search["query_title"])
         if title:
            curr_id = title.group(1)
            cds_id_ = cds_id + "_"
            if cds_id == curr_id or cds_id_ in curr_id:
               hits = search["hits"]
               for hit in hits:
                  hsps = hit["hsps"][0]
                  if hsps["evalue"] <= e_value_thresh:
                     alignment = {}
                     description = hit["description"][0]
                     alignment['accession'] = description["accession"]
                     alignment["title"] = description["title"]
                     alignment["sciname"] = description["sciname"]
                     alignment["evalue"] = '{:0.2e}'.format(hsps["evalue"])
                     alignment["query_from"] = hsps["query_from"]
                     alignment["query_to"] = hsps["query_to"]
                     alignment["hit_from"] = hsps["hit_from"]
                     alignment["hit_to"] = hsps["hit_to"]
                     alignment["percent_identity"] = round(hsps["identity"] / hsps["align_len"] * 100, 2)
                     blast_results.append(alignment)

   return blast_results

def parse_blast_multiple(blast_file, cds_id, e_value_thresh):
   blast_results = {}
   with open(blast_file) as f:
      blasts = json.load(f)["BlastOutput2"]
      for blast in blasts:
         search = blast["report"]["results"]["search"]
         title = re.search("(FERN_\d+_*\d*), (\d+)-\d+", search["query_title"])
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
                     alignment["evalue"] = '{:0.2e}'.format(hsps["evalue"])
                     alignment["query_from"] = hsps["query_from"]
                     alignment["query_to"] = hsps["query_to"]
                     alignment["hit_from"] = hsps["hit_from"]
                     alignment["hit_to"] = hsps["hit_to"]
                     alignment["percent_identity"] = round(hsps["identity"] / hsps["align_len"] * 100, 2)
                     blast_results[curr_start].append(alignment)

   return blast_results


#-------------------------------------------------------------------
# def run_blast(fasta_file, start, stop, e_value_thresh):
#    start = start - 1
#    record = SeqIO.read(fasta_file, "fasta")
#    genome = record.seq
#    sequence = Seq.translate(sequence=genome[start: stop], table=11)
#    print("In `run_blast`")
#    result_handle = NCBIWWW.qblast("blastp", "nr", sequence[:-1])
#    print("Finished BLAST query")

#    with open("blast.xml", "w") as out_handle:
#       print("Writing BLAST results to file...")
#       blast_results = result_handle.read()
#       out_handle.write(blast_results)
#    result_handle.close()

#    blast = parse_blast(e_value_thresh)
#    return blast

# def parse_blast_xml(blast_file, cds_id, e_value_thresh):
#    blast = []
#    for record in NCBIXML.parse(open(blast_file)): 
#       for iteration in record.descriptions:
#          print(iteration.title + "\n")
#          if cds_id in iteration.title:
#             for align in record.alignments:
#                if cds_id in align.title:
#                   alignment = {}
#                   for hsp in align.hsps: 
#                      if hsp.expect <= e_value_thresh: 
#                         hit_definition = re.search('([A-Za-z\s]+\[[\.-_A-Za-z\s]+\]).+', align.hit_def)
#                         if hit_definition:
#                               alignment['hit_def'] = hit_definition.group(1)
#                         else:
#                               alignment['hit_def'] = align.hit_def
#                         # alignment['title'] = align.title
#                         alignment['e_value'] = float(hsp.expect)
#                         # alignment['query'] = hsp.query
#                         alignment['query_start'] = hsp.query_start
#                         alignment['query_end'] = hsp.query_end
#                         # alignment['match'] = hsp.match
#                         # alignment['sbjct'] = hsp.sbjct
#                         alignment['sbjct_start'] = hsp.sbjct_start
#                         alignment['sbjct_end'] = hsp.sbjct_end
#                         blast.append(alignment)
#    return blast

