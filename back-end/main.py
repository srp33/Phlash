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
from blast import *
import blast
from annotations import *
import annotations
from annotations_cds import *
import annotations_cds
from gene_map import *
import gene_map
from genbank import *
import genbank
from settings import *
import settings
import threading
import time
import subprocess
import asyncio

# Configuration
ROOT = os.path.dirname(os.path.abspath(__file__))

# Instantiates the app
app = Flask(__name__)
app.config.from_object(__name__)
db.init_app(app)

# Enables CORS
CORS(app)

# Database instantiation
DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', "Phlash.db"))
app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
with app.app_context():
    # db.drop_all()
    db.create_all()

# routers ------------------------------------------------------------------
@app.route('/phlash_api/annotations/geneMap/<phage_id>/', methods=['GET'])
def gene_map(phage_id):
    """
    Builds and returns the gene map.
    """
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', phage_id, 'uploads')
    if request.method == "GET":
        return jsonify(get_map(phage_id, UPLOAD_FOLDER))

@app.route('/phlash_api/settings/<phage_id>/<payload>/', methods=['GET'])
def settings(phage_id, payload):
    """
    Updates default settings.
    """
    if request.method == "GET":
        if payload == "none":
            return jsonify(get_settings(phage_id))
        else:
            return jsonify(update_settings(phage_id, payload))

@app.route('/phlash_api/share/<phage_id>/<email>/', methods=['POST'])
def share(phage_id, email):
    """
    Adds another user to the phage_id for view access only.
    """
    if request.method == "POST":
        original_user = db.session.query(Users).filter_by(id=phage_id).first()
        if original_user != None and db.session.query(Users).filter_by(user=email).first() and not db.session.query(Users).filter_by(user=email).filter_by(id=phage_id).first():
            new_user = Users(user = email,
                        phage_id = original_user.phage_id,
                        creation_date = "view",
                        deletion_date = original_user.user,
                        id = phage_id)
            db.session.add(new_user)
            db.session.commit()
            message = "The user " + email + " was successfully given view permissions."
            return jsonify(message)
        elif db.session.query(Users).filter_by(user=email).filter_by(id=phage_id).first():
            message = "The user " + email + " already has permission to view this phage."
            return jsonify(message)
        else:
            message = "The user " + email + " does not have an account."
            return jsonify(message)

@app.route('/phlash_api/test', methods=['GET'])
def test():
    return jsonify("Hello, world!")

@app.route('/phlash_api/home/<current_user>/<phage_id>', methods=['POST', 'DELETE'])
def check_phage(phage_id, current_user):
    """
    API endpoint for '/home/:current_user/:phage_id'.
    POST method removes users that have existed for more than 90 days,
                creates a new user if it doesn't exist,
                else gets informations for existing user.
    DELETE method removes all data associated with a phage ID.
    """
    exists = False
    users = db.session.query(Users).filter_by(phage_id=phage_id).filter_by(user=current_user)
    if users:
        for user in users:
            if user.creation_date != 'view':
                exists = True
                break
    if request.method == "POST":
        if not exists:
            return jsonify(home.check_phage_id(current_user, phage_id, app))
        else:
            response_object['id_status'] = "ID already exists. Please continue below."
            return jsonify(response_object)
    else:
        user = db.session.query(Users).filter_by(id=phage_id).filter_by(user=current_user).first()
        if user is not None:
            return jsonify(remove_user(user.id, user))

@app.route('/phlash_api/check_upload/<phage_id>', methods=['GET'])
def check_upload(phage_id):
    """
    API endpoint for '/upload/:phage_id'.
    POST method uploads files accordingly and removes files if necessary.
    """
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', phage_id, 'uploads')

    return jsonify(check_uploaded_files(UPLOAD_FOLDER, phage_id))

@app.route('/phlash_api/check_user/<current_user>/<phage_id>', methods=['GET'])
def check_user(current_user, phage_id):
    """
    API endpoint for '/check_user/:current_user/:phage_id'.
    GET method checks to see if a user ID matches a given phage ID.
    """
    if db.session.query(Users).filter_by(id=phage_id).filter_by(user=current_user).first() is None:
        return jsonify("fail")
    elif db.session.query(Users).filter_by(id=phage_id).filter_by(user=current_user).first().creation_date == "view":
        response_object = {}
        response_object['view'] = True
        response_object['phage_id'] = db.session.query(Users).filter_by(id=phage_id).filter_by(user=current_user).first().phage_id
        return jsonify(response_object)
    return jsonify(db.session.query(Users).filter_by(id=phage_id).filter_by(user=current_user).first().phage_id)

@app.route('/phlash_api/get_user_data/<email>', methods=['GET'])
def get_user_data(email):
    """
    API endpoint for '/get_user_data/:email'.
    GET method finds all of the phages and their data associated with a user.
    """
    return jsonify(get_phage_data(email))

@app.route('/phlash_api/upload/<phage_id>/<file_method>/<file_path>', methods=['POST'])
def upload(phage_id, file_method, file_path):
    """
    API endpoint for '/upload/:phage_id'.
    POST method uploads files accordingly and removes files if necessary.
    """
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', phage_id, 'uploads')

    if file_method == "display":
        return jsonify(display_files(UPLOAD_FOLDER))

    elif file_method == "delete":
        return jsonify(delete_file(phage_id, file_path, UPLOAD_FOLDER))

    elif file_method == "uploadFasta":
        dropzone_fasta(UPLOAD_FOLDER, request, phage_id)
        return jsonify("success")

@app.route('/phlash_api/blast/<phage_id>/<file_method>/<file_path>', methods=['GET', 'POST'])
def blast(phage_id, file_method, file_path):
    """
    API endpoint for '/blast/:phage_id/:file_method/:file_path'.
    GET method determines if all files have been downloaded and uploaded or
        auto annotates the genome.
    POST method downloads fasta file for BLAST input or
        creates the input files or
        uploads json file of BLAST output or
        deletes json file of BLAST output or
        returns json file names of BLAST output or
        returns the number of needed json files or
        checks for partially uploaded files.

    """
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', phage_id, 'uploads')

    if request.method == "GET":
        if file_method == "checkFiles":
            return jsonify(find_blast_zip(phage_id))

        elif file_method == "autoAnnotate":
            return jsonify(add_annotation_task(phage_id, UPLOAD_FOLDER))
        
    if request.method == "POST":
        if file_method == "downloadInput":
            return download_blast_input(phage_id)

        elif file_method == "createInput":
            return jsonify(add_blast_input_task(UPLOAD_FOLDER, phage_id))

        elif file_method == "displayOutput":
            return jsonify(get_blast_output_names(phage_id, UPLOAD_FOLDER, file_path))

        elif file_method == "deleteOutput":
            return jsonify(delete_blast_output(phage_id, UPLOAD_FOLDER, file_path))

        elif file_method == "numFiles":
            return jsonify(get_num_blast_files(phage_id))

        elif file_method == "drop":
            dropzone(phage_id, UPLOAD_FOLDER, request)
            return jsonify("success")

        elif file_method == "deleteBlastResults":
            return jsonify(delete_all_blast(phage_id, UPLOAD_FOLDER))

        else:
            return jsonify(new_file(phage_id, file_path, file_method))

@app.route('/phlash_api/annotations/<phage_id>/<file_method>', methods=['GET', 'POST', 'PUT'])
def annotate_data(phage_id, file_method):
    """
    Parses blast data and returns annotation data.
    GET method shows all the annotation predictions with a status and action item for each.
    PUT method adds a CDS
    """
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', phage_id, 'uploads')

    if request.method == "GET":
        if file_method == "check":
            return jsonify(check_blast_task(phage_id))
        elif file_method == "blast":
            return jsonify(add_blast_task(phage_id, UPLOAD_FOLDER))
        else:
            return jsonify(get_annotations_data(phage_id))

    if request.method == "PUT":
        return jsonify(add_cds(request, UPLOAD_FOLDER, phage_id))

@app.route('/phlash_api/genbank/<phage_id>/<file_method>', methods=['GET', 'POST', 'PUT'])
def create_genbank(phage_id, file_method):
    """
    Creates the Genbank file from data stored in database.
    POST method creates and returns the genbank file.
    """
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', phage_id, 'uploads')
    response_object = {'status': 'success'}

    if request.method == "POST":
        return get_genbank(UPLOAD_FOLDER, phage_id, request)

@app.route('/phlash_api/annotations/cds/<phage_id>/<cds_id>', methods=['GET', 'PUT'])
def cds_annotation(phage_id, cds_id):
    """
    Annotation information for each CDS.
    GET method gets cds, its left options, blast results, and graph data.
    PUT method updates the left position and function if the user chooses to do so.
    """
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', phage_id, 'uploads')
    response_object = {'status': 'success'}

    if request.method == "GET":
        return jsonify(get_cds_data(phage_id, UPLOAD_FOLDER, cds_id))
                
    if request.method == "PUT":
        return jsonify(annotate_cds(phage_id, request, cds_id, UPLOAD_FOLDER))

if __name__ == '__main__':
    app.run(debug=False)