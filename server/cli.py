import subprocess

# with open('output.txt', 'w') as f:
p1 = subprocess.run(['./ncbi-blast-2.10.0+/bin/blastp', '-query', 'FASTA'], capture_output=True, text=True, check=True)

# p1.returncode (0 means no error, 1 means error)
# p1.stdout
# p1.stderr
# print(p1.stdout)

./blastp -query first_50.fasta -db nr/nr -out local_db.out -evalue 1e-3 -outfmt "7 delim=, qseqid sseqid qstart qend sstart send evalue pident nident" -num_threads 4 &