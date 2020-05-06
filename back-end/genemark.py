import subprocess

def run_genemark(fasta_file_path):
    """
    This function invokes the GeneMarkS utilities and returns
    .gdata and .ldata files.
    """

    result = subprocess.run(["/genemark_suite_linux_64/gmsuite/gc", fasta_file_path], stdout=subprocess.PIPE)
    gc_percent = result.stdout.decode("utf-8").split(" ")[3]
    gc_percent = "{:d}".format(round(float(gc_percent)))

    subprocess.run(["/genemark_suite_linux_64/gmsuite/gm", "-D", "-g", "0", "-s", "1", "-m", "/genemark_suite_linux_64/gmsuite/heuristic_mat/heu_11_{}.mat".format(gc_percent), "-v", fasta_file_path])

    gdata_file_path = "{}.gdata".format(fasta_file_path)
    ldata_file_path = "{}.ldata".format(fasta_file_path)

    return read_file(gdata_file_path), read_file(ldata_file_path)

def read_file(file_path):
    with open(file_path) as file_handle:
        return file_handle.read()

## Can be used for testing purposes
#import sys
#fasta_file_path = sys.argv[1]
#
#gdata, ldata = run_genemark(fasta_file_path)
#
#print("gdata:")
#print(gdata)
#print("ldata:")
#print(ldata)
