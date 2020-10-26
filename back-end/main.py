"""
Main back-end script of Flask web application.
"""
from flask_cors import CORS
from flask import *
import os
from models import *
import models
from home import *
import home
from upload import *
import upload
from dna_master import *
import dna_master
from blast import *
import blast
from annotations import *
import annotations
from annotations_cds import *
import annotations_cds

# Configuration
ROOT = os.path.dirname(os.path.abspath(__file__))
FASTA_EXTENSIONS = set(['.fasta', '.fna'])
GENBANK_EXTENSIONS = set(['.gb', '.gbk'])
GDATA_EXTENSIONS = set(['.gdata'])
LDATA_EXTENSIONS = set(['.ldata'])
BLAST_EXTENSIONS = set(['.json'])

# Instantiates the app
app = Flask(__name__)
app.config.from_object(__name__)
db.init_app(app)

# Enables CORS
CORS(app, resources={r'/*': {'origins': '*'}})

# routers ------------------------------------------------------------------
@app.route('/phlash_api/test', methods=['GET'])
def test():
    return jsonify("Hello, world!")

@app.route('/phlash_api/home/<phage_id>', methods=['POST'])
def check_phage_id(phage_id):
    """
    API endpoint for '/home/:phage_id'.
    POST method removes users that have existed for more than 90 days,
                creates a new user if it doesn't exist,
                else gets informations for existing user.
    """
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
    return jsonify(check_uploaded_files(UPLOAD_FOLDER))

@app.route('/phlash_api/upload/<current_user>/<file_method>/<file_path>', methods=['POST'])
def upload(current_user, file_method, file_path):
    """
    API endpoint for '/upload/:current_user'.
    POST method uploads files accordingly and removes files if necessary.
    """
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')
    DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', current_user, f"{current_user}.db"))
    app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE

    if file_method == "upload":
        return jsonify(upload_file(request, UPLOAD_FOLDER))

    elif file_method == "display":
        return jsonify(display_files(UPLOAD_FOLDER))
       
    elif file_method == "download":
        return jsonify(download_file(file_path, UPLOAD_FOLDER))

    elif file_method == "delete":
        return jsonify(delete_file(file_path, UPLOAD_FOLDER))

@app.route('/phlash_api/dnamaster/<current_user>', methods=['GET', 'POST'])
def dnamaster(current_user):
    """
    API endpoint for '/dnamaster/:current_user'.
    GET method querys database for parsed DNA Master data (from uploaded GenBank file).
    """
    DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', current_user, f"{current_user}.db"))
    app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE

    if request.method == "GET":
        return jsonify(get_all_cds())

@app.route('/phlash_api/dnamaster/<current_user>/<cds_id>', methods=['PUT', 'DELETE'])
def dnamaster_cds(current_user, cds_id):
    """
    API endpoint for '/dnamaster/:current_user/:cds_id'.
    PUT method updates a CDS.
    DELETE method deletes a CDS.
    """
    DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', current_user, f"{current_user}.db"))
    app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')

    if request.method == "PUT":
        return jsonify(update_cds(request, UPLOAD_FOLDER))

    if request.method == 'DELETE':
        return jsonify(delete_cds(cds_id, UPLOAD_FOLDER))

# TODO: Continue checking code from here.
@app.route('/phlash_api/blast/<current_user>/<file_method>/<file_path>', methods=['GET', 'POST'])
def blast(current_user, file_method, file_path):
    """
    API endpoint for '/blast/:current_user/:file_method/:file_path'.
    GET method determines if the Blast zip has already been downloaded
    POST method downloads fasta file for BLAST input or
        uploads json file of BLAST output or
        deletes json file of BLAST output or
        downloads json file of BLAST output or
        returns json file names of BLAST output
    """
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')
    DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', current_user, f"{current_user}.db"))
    app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE

    if request.method == "GET":
        return jsonify(find_blast_zip(current_user))
        
    if request.method == "POST":
        if file_method == "downloadInput":
            return download_blast_input(UPLOAD_FOLDER, current_user)

        elif file_method == "uploadOutput":
            return jsonify(upload_blast_output(UPLOAD_FOLDER, request))

        elif file_method == "displayOutput":
            return jsonify(get_blast_output_names(UPLOAD_FOLDER))

        elif file_method == "downloadOutput":
            return jsonify(download_blast_output(UPLOAD_FOLDER, file_path))

        elif file_method == "deleteOutput":
            return jsonify(delete_blast_output(UPLOAD_FOLDER, file_path))

        elif file_method == "numFiles":
            return get_num_blast_files(current_user)

@app.route('/phlash_api/annotations/<current_user>', methods=['GET', 'POST'])
def annotate_data(current_user):
    """
    Compares DNA Master's predictions against GeneMark's.
    GET method shows all the DNA Master predictions with a status and action item for each.
    POST method returns the genbank file
    """
    DATABASE = "sqlite:///{}".format(os.path.join(ROOT,
                                                  'users', current_user, f"{current_user}.db"))
    app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')
    response_object = {'status': 'success'}

    if request.method == "GET":
        return jsonify(get_dnamaster_data())

    if request.method == "POST":
        return get_genbank(UPLOAD_FOLDER)

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
        return jsonify(get_cds_data(UPLOAD_FOLDER, cds_id))
                
    if request.method == "PUT":
        return jsonify(annotate_cds(request, cds_id))

if __name__ == '__main__':
    app.run()