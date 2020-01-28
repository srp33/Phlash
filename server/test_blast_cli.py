from Bio import SeqIO, Seq, SeqFeature
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import *
import re


gb_file = "uploads/fern.gb"
out_file = "test/fern_modified.gb"
gene_name = 'FERN_2'
modified_start = 333
strand = '+'

genes_to_modify = {gene_name: {'start': modified_start, 'strand': strand}}

final_features = []
for record in SeqIO.parse(open(gb_file, "r"), "genbank"):
    for feature in record.features:
        if feature.type == "gene" or feature.type == "CDS":
            locus_tag = feature.qualifiers["locus_tag"][0]
            if locus_tag in genes_to_modify:
                new_start = genes_to_modify[locus_tag]['start']
                if feature.location.strand == 1:
                    feature.location = SeqFeature.FeatureLocation(SeqFeature.ExactPosition(new_start - 1),
                                                                    SeqFeature.ExactPosition(
                                                                        feature.location.end.position),
                                                                    feature.location.strand)
                else:
                    feature.location = SeqFeature.FeatureLocation(
                        SeqFeature.ExactPosition(feature.location.start.position),
                        SeqFeature.ExactPosition(new_start), feature.location.strand)
        final_features.append(feature)  # Append final features
    record.features = final_features
    with open(out_file, "w") as new_gb:
        SeqIO.write(record, new_gb, "genbank")



            