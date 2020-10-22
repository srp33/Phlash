"""
Main back-end script of Flask web application.
"""
from werkzeug.utils import secure_filename
from Bio import SeqIO, Seq
from pathlib import Path
from contextlib import closing
from zipfile import ZipFile
from flask_cors import CORS
from flask import *
from models import *
import os
import shutil
import datetime
import uuid
import pandas as pd
import arrow
import annotate
from datetime import datetime
import services.home as home

# configuration
ROOT = os.path.dirname(os.path.abspath(__file__))
FASTA_EXTENSIONS = set(['.fasta', '.fna'])
GENBANK_EXTENSIONS = set(['.gb', '.gbk'])
GDATA_EXTENSIONS = set(['.gdata'])
LDATA_EXTENSIONS = set(['.ldata'])
BLAST_EXTENSIONS = set(['.json'])
# USERS = []
# CURRENT_USER = ""

# instantiate the app
app = Flask(__name__)
app.config.from_object(__name__)
db.init_app(app)

# enable CORS
CORS(app, resources={r'/*': {'origins': '*'}})

# routers ------------------------------------------------------------------
@app.route('/phlash_api/test', methods=['GET'])
def test():
    return jsonify("Hello, world!")


"""
API endpoint for '/home/:phage_id'.
POST method removes users that have existed for more than 90 days,
            creates a new user if it doesn't exist,
            else gets informations for existing user.
"""
@app.route('/phlash_api/home/<phage_id>', methods=['POST'])
def check_phage_id(phage_id):
    DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', phage_id, f"{phage_id}.db"))
    app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
    return jsonify(home.check_phage_id(phage_id, app))

@app.route('/phlash_api/check_upload/<current_user>', methods=['GET'])
def check_upload(current_user):
    """
    API endpoint for '/upload/:current_user'.
    POST method uploads files accordingly and removes files if necessary.
    """
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')
    DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', current_user, f"{current_user}.db"))
    app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
    response_object = {}

    # check if respective files for fasta and genbank are uploaded
    existing_files = []
    for filename in os.listdir(os.path.join(ROOT, 'users', current_user, 'uploads')):
        ext = os.path.splitext(filename)[1].lower()
        if ext in FASTA_EXTENSIONS:
            existing_files.append("fasta")
        elif ext in GENBANK_EXTENSIONS:
            existing_files.append("genbank")
        elif ext in GDATA_EXTENSIONS:
            existing_files.append("gdata")
        elif ext in LDATA_EXTENSIONS:
            existing_files.append("ldata")

    response_object["fasta"] = True if "fasta" in existing_files and \
                                        "gdata" in existing_files and \
                                        "ldata" in existing_files else False
    response_object["genbank"] = True if "genbank" in existing_files else False

    return jsonify(response_object)

@app.route('/phlash_api/upload/<current_user>/<file_method>/<file_path>', methods=['POST'])
def upload(current_user, file_method, file_path):
    """
    API endpoint for '/upload/:current_user'.
    POST method uploads files accordingly and removes files if necessary.
    """
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')
    DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', current_user, f"{current_user}.db"))
    app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
    response_object = {}
    if file_method == "upload":
        if 'file' not in request.files:
            response_object["status"] = "'file' not in request.files"
        else:
            file = request.files['file']
            fileType = request.form['fileType']
            if file:
                file_name = secure_filename(file.filename)
                if (fileType == "fasta" and allowed_file(file_name, FASTA_EXTENSIONS)) or \
                    (fileType == "genbank" and allowed_file(file_name, GENBANK_EXTENSIONS)):

                    # overwrite existing files with specific extension with newly uploaded files
                    for existing_file in os.listdir(UPLOAD_FOLDER):
                        ext = os.path.splitext(existing_file)[1].lower()
                        if (fileType == "fasta" and ext in ['.fasta', '.fna', '.ldata', '.gdata', '.lst', '.ps']) or \
                            (fileType == "genbank" and ext in GENBANK_EXTENSIONS):
                            os.remove(os.path.join(UPLOAD_FOLDER, existing_file))
                            print(f" * removed {existing_file}")

                    file.save(os.path.join(UPLOAD_FOLDER, file_name))
                    response_object["uploaded"] = file_name
                    print(' * uploaded', file_name)
#Below is create your own fasta file stuff
                    # if fileType == "genbank":
                    #     with open(os.path.join(UPLOAD_FOLDER, file_name), 'r') as gFile:
                    #         name = os.path.join(UPLOAD_FOLDER, current_user+".fasta")
                    #         with open (name, 'w') as fFile:
                    #             fFile.write(">" + current_user)
                    #             startWriting = False
                    #             content = []
                    #             for line in gFile:
                    #                 if line.startswith("ORIGIN"):
                    #                     startWriting = True
                    #                     continue
                    #                 if startWriting:
                    #                     content.append("".join(line[10:-1].strip().split()))
                    #             print(content)
                    #             fFile.writelines(content)
                    #             file.save(name)
                    # parse appropriate files as soon as uploaded
                    # FIXME: check file conent before parsing.
                    file_ext = os.path.splitext(file_name)[1].lower()
                    if file_ext in FASTA_EXTENSIONS:  # run genemark and parse ldata
                        fasta_file = get_file_path("fasta", UPLOAD_FOLDER)
                        annotate.run_genemark(fasta_file)
                        try:
                            ldata_file = get_file_path("ldata", UPLOAD_FOLDER)
                        except FileNotFoundError:
                            print("GeneMark did not save ldata file properly")
                            return("error")
                        annotate.parse_genemark_ldata(ldata_file)
                    elif file_ext in GENBANK_EXTENSIONS:  # parse genbank
                        genbank_file = get_file_path("genbank", UPLOAD_FOLDER)
                        annotate.parse_dnamaster_genbank(genbank_file)
                else:
                    response_object["not_allowed"] = file.filename
            else:
                response_object["status"] = "error"

    elif file_method == "display":
        response_object["fasta_file"] = "Not found"
        response_object["genbank_file"] = "Not found"
        for file in os.listdir(UPLOAD_FOLDER):
            if (file.endswith(".fasta") or file.endswith(".fna") or file.endswith(".txt")):
                response_object["fasta_file"] = file
                print(response_object["fasta_file"])
            elif (file.endswith(".gb") or file.endswith(".gbk")):
                response_object["genbank_file"] = file

    elif file_method == "download":
        try:
            response_object["file_data"] = open(os.path.join(UPLOAD_FOLDER, file_path)).read()
            response_object["status"] = "success"
        except:
            print("error")
            response_object["status"] = "error"

    elif file_method == "delete":
        try:
            os.remove(os.path.join(UPLOAD_FOLDER, file_path))
            DNAMaster.query.delete()
            GeneMark.query.delete()
            Files.query.delete()
            response_object["status"] = "success"
        except:
            print("error")
            response_object["status"] = "error"
    return jsonify(response_object)
@app.route('/phlash_api/dnamaster/<current_user>', methods=['GET', 'POST'])
def dnamaster(current_user):
    """
    API endpoint for '/dnamaster/:current_user'.
    GET method querys database for parsed DNA Master data (from uploaded GenBank file).
    POST method adds a new CDS to the data.
    """
    DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', current_user, f"{current_user}.db"))
    app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
    response_object = {}

    if request.method == "GET":
        dnamaster = []
        for cds in db.session.query(DNAMaster).order_by(DNAMaster.start):
            dnamaster.append({'id': cds.id,
                              'start': cds.start,
                              'stop': cds.stop,
                              'strand': cds.strand})
        response_object['dnamaster'] = dnamaster

    if request.method == "POST":
        post_data = request.get_json()
        cds = DNAMaster(id = post_data.get('id'),
                        start = int(post_data.get('start')),
                        stop = int(post_data.get('stop')),
                        strand = post_data.get('strand'),
                        function = "None",
                        status = "None")
        start_exists = DNAMaster.query.filter_by(start=post_data.get('start')).first()
        stop_exists = DNAMaster.query.filter_by(stop=post_data.get('stop')).first()
        id_exists = DNAMaster.query.filter_by(id=post_data.get('id')).first()
        if cds.start > cds.stop or cds.start == cds.stop:
            response_object['message'] = 'Start is not less than start. CDS not added.'
        elif start_exists and stop_exists:
            if start_exists.id == stop_exists.id:
                response_object['message'] = 'Start and stop already exists. CDS not added.'
        elif id_exists:
            if id_exists.id == cds.id:
                response_object['message'] = 'ID already exists. CDS not added.'
        else:
            db.session.add(cds)
            db.session.commit()
            response_object['message'] = 'CDS added!'

    return jsonify(response_object)


@app.route('/phlash_api/dnamaster/<current_user>/<cds_id>', methods=['PUT', 'DELETE'])
def dnamaster_cds(current_user, cds_id):
    """
    API endpoint for '/dnamaster/:current_user/:cds_id'.
    PUT method updates a CDS.
    DELETE method deletes a CDS.
    """
    DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', current_user, f"{current_user}.db"))
    app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
    response_object = {}

    if request.method == "PUT":
        put_data = request.get_json()
        cds = DNAMaster.query.filter_by(id=put_data.get('id')).first()
        if cds:
            cds.start = put_data.get('start')
            cds.stop = put_data.get('stop')
            cds.strand = put_data.get('strand')
            db.session.commit()
            response_object['message'] = 'CDS updated!'
        else:
            response_object['message'] = 'Error: CDS could not get updated.'

    if request.method == 'DELETE':
        if DNAMaster.query.filter_by(id=cds_id).first():
            DNAMaster.query.filter_by(id=cds_id).delete()
            db.session.commit()
            response_object['message'] = 'CDS removed!'
        else:
            response_object['message'] = 'Error: CDS could not be deleted.'

    return jsonify(response_object)

# TODO: Continue checking code from here.
@app.route('/phlash_api/blast/<current_user>/<file_method>/<file_path>', methods=['GET', 'POST'])
def blast(current_user, file_method, file_path):
    """
    API endpoint for '/blast/:current_user'.
    POST method downloads fasta file for BLAST input or
        uploads json file of BLAST output
    """
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')
    DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', current_user, f"{current_user}.db"))
    app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
    response_object = {}

    if request.method == "GET":
        response_object["blast_uploaded"] = False
        response_object["blast_downloaded"] = False

        # # check if respective file(s) for blast are uploaded
        # for filename in os.listdir(os.path.join(ROOT, 'users', current_user, 'uploads')):
        #     if os.path.splitext(filename)[1].lower() in BLAST_EXTENSIONS:
        #         response_object["blast_uploaded"] = True

        # check if blast input file(s) downloaded
        for filename in os.listdir(os.path.join(ROOT, 'users', current_user)):
            if filename.endswith('.zip'):
                response_object["blast_downloaded"] = True

    if request.method == "POST":
        if file_method == "downloadInput":
            print("Starting comparisons")
            annotate.compare()
            fasta_file = get_file_path("fasta", UPLOAD_FOLDER)
            genemark_gdata_file = get_file_path("gdata", UPLOAD_FOLDER)
            print(datetime.now())
            num_files_downloaded = annotate.create_blast_fasta(current_user, fasta_file, genemark_gdata_file)
            print(datetime.now())
            f = open(os.path.join(ROOT, 'users', current_user, f"{current_user}_blast.zip"), "rb")
            #return {'num_files_downloaded': num_files_downloaded, 'zipfile': f.read()}
            return f.read()
        elif file_method == "displayOutput":
            file_names = []
            for file in os.listdir(UPLOAD_FOLDER):
                if file.endswith(".json"):
                    file_names.append(file)
            print(file_names)
            response_object["file_names"] = file_names
        elif file_method == "downloadOutput":
            try:
                response_object["file_data"] = open(os.path.join(UPLOAD_FOLDER, file_path)).read()
                response_object["status"] = "success"
            except:
                print("error")
                response_object["status"] = "error"
        elif file_method == "deleteOutput":
            try:
                os.remove(os.path.join(UPLOAD_FOLDER, file_path))
                response_object["status"] = "success"
            except:
                print("error")
                response_object["status"] = "error"
        elif file_method == "uploadOutput":
            if 'file' not in request.files:
                print(request.files)
                response_object["status"] = "'file' not in request.files"
                print("in fail")
            else:
                print("in success")
                file = request.files['file']
                if file:
                    file_name = secure_filename(file.filename)
                    if allowed_file(file_name, set(['.json'])):

                        file_ext = file_name.rsplit('.', 1)[1].lower()
                        for existing_file in os.listdir(UPLOAD_FOLDER):
                            if existing_file.endswith(file_name):
                                os.remove(os.path.join(
                                    UPLOAD_FOLDER, existing_file))

                        file.save(os.path.join(UPLOAD_FOLDER, file_name))
                        response_object["uploaded"] = file_name
                        print(' * uploaded', file_name)
                    else:
                        response_object["not_allowed"] = file.filename
                else:
                    response_object["status"] = "error"               
        elif file_method == "numFiles":
            for filename in os.listdir(os.path.join(ROOT, 'users', current_user)):
                if filename.endswith('.zip'):
                    with closing(ZipFile(os.path.join(ROOT, 'users', current_user, filename))) as archive:
                        num_blast_files = len(archive.infolist())
                        return str(num_blast_files)

    return jsonify(response_object)


@app.route('/phlash_api/annotations/<current_user>', methods=['GET', 'POST'])
def annotate_data(current_user):
    """
    Compares DNA Master's predictions against GeneMark's.
    GET method shows all the DNA Master predictions with a status and action item for each.
    """
    DATABASE = "sqlite:///{}".format(os.path.join(ROOT,
                                                  'users', current_user, f"{current_user}.db"))
    app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')
    response_object = {'status': 'success'}

    if request.method == "GET":
        dnamaster = []
        for cds in db.session.query(DNAMaster).order_by(DNAMaster.start):
            dnamaster.append({'id': cds.id,
                              'start': cds.start,
                              'stop': cds.stop,
                              'strand': cds.strand,
                              'function': cds.function,
                              'status': cds.status})
        response_object['dnamaster'] = dnamaster

    if request.method == "POST":
        # -------downloading GENBANK file----------
        gb_file = get_file_path("genbank", UPLOAD_FOLDER)
        fasta_file = get_file_path("fasta", UPLOAD_FOLDER)
        out_file = annotate.modify_genbank(gb_file, fasta_file)
        f = open(out_file, "r")
        return f.read()

    return jsonify(response_object)


@app.route('/phlash_api/annotations/cds/<current_user>/<cds_id>', methods=['GET', 'POST', 'PUT', 'DELETE'])
def cds_annotation(current_user, cds_id):
    """
    Annotation information for each CDS.
    GET method gets cds, its start options, blast results, and graph data.
    PUT method updates the start position and function if the user chooses to do so.
    """
    DATABASE = "sqlite:///{}".format(os.path.join(ROOT,
                                                  'users', current_user, f"{current_user}.db"))
    app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')
    response_object = {'status': 'success'}

    if request.method == "GET":
        cds = DNAMaster.query.filter_by(id=cds_id).first()
        response_object['cds'] = {'id': cds.id,
                                  'start': cds.start,
                                  'stop': cds.stop,
                                  'strand': cds.strand,
                                  'function': cds.function,
                                  'status': cds.status}

        # fasta_file = get_file_path("fasta", UPLOAD_FOLDER)
        print(cds)
        
        starts = [int(start) for start in cds.start_options.split(",")]
        stops = [int(stop) for stop in cds.stop_options.split(",")]
        print(starts)
        print(stops)
        response_object['start_options'] = starts
        response_object['stop_options'] = stops
        blast_files = get_file_path("blast", UPLOAD_FOLDER)
        E_VALUE_THRESH = 1e-7
        print("parsing blast")
        blast_results = annotate.parse_blast_results(blast_files, cds.id, E_VALUE_THRESH)
        response_object['blast'] = blast_results

        genemark_gdata_file = get_file_path("gdata", UPLOAD_FOLDER)
        gdata_df = pd.read_csv(genemark_gdata_file, sep='\t', skiprows=16)
        gdata_df.columns = ['Base', '1', '2', '3', '4', '5', '6']
        gdata_df = gdata_df[gdata_df.Base.isin(
            range(min(starts) - 100, cds.stop + 100))]
        response_object['x_data'] = gdata_df["Base"].to_list()
        response_object['y_data_1'] = gdata_df["1"].to_list()
        response_object['y_data_2'] = gdata_df["2"].to_list()
        response_object['y_data_3'] = gdata_df["3"].to_list()
        response_object['y_data_4'] = gdata_df["4"].to_list()
        response_object['y_data_5'] = gdata_df["5"].to_list()
        response_object['y_data_6'] = gdata_df["6"].to_list()

        dnamaster = []
        reachedCDS = False
        for cds in db.session.query(DNAMaster).order_by(DNAMaster.start):
            if reachedCDS and cds.function == "None selected":
                response_object['nextCDS'] = cds.id
                print(cds.id)
                break
            elif cds.id == cds_id:
                reachedCDS = True
                
    if request.method == "PUT":
        put_data = request.get_json()
        cds = DNAMaster.query.filter_by(id=cds_id).first()
        if cds:
            cds.start = put_data.get('start')
            cds.function = put_data.get('function')
            cds.status = put_data.get('status')
            db.session.commit()
            response_object['message'] = 'CDS updated!'
        else:
            response_object['message'] = 'CDS did not update.'

    if request.method == "DELETE":
        if DNAMaster.query.filter_by(id=cds_id).first():
            DNAMaster.query.filter_by(id=cds_id).delete()
            db.session.commit()
            response_object['message'] = 'CDS removed!'
        else:
            response_object['message'] = 'Error: CDS could not be deleted.'

    return jsonify(response_object)


if __name__ == '__main__':
    app.run()


# ---------- HELPER FUNCTIONS ----------
def allowed_file(filename, allowed_extensions):
    """
    Check if file extension is acceptable.
    """
    return '.' in filename and os.path.splitext(filename)[1].lower() in allowed_extensions

def get_file_path(preference, upload_directory):
    """
    Gets path of required file.
    """
    blast_files = []
    for filename in os.listdir(upload_directory):
        file_ext = os.path.splitext(filename)[1].lower()
        if preference == "fasta":
            if file_ext in FASTA_EXTENSIONS:
                return os.path.join(upload_directory, filename)
        elif preference == "genbank":
            if file_ext in GENBANK_EXTENSIONS:
                return os.path.join(upload_directory, filename)
        elif preference == "ldata":
            if file_ext in LDATA_EXTENSIONS:
                return os.path.join(upload_directory, filename)
        elif preference == "gdata":
            if file_ext in GDATA_EXTENSIONS:
                return os.path.join(upload_directory, filename)
        elif preference == "blast":
            if file_ext in BLAST_EXTENSIONS:
                blast_files.append(os.path.join(upload_directory, filename))
        else:
            return("Couldn't find file.")
    return blast_files
