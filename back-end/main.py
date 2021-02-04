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
from annotations_gene_map import *
import annotations_gene_map
from settings import *
import settings
import threading
import time

# Configuration
ROOT = os.path.dirname(os.path.abspath(__file__))

# Instantiates the app
app = Flask(__name__)
app.config.from_object(__name__)
db.init_app(app)

# Enables CORS
CORS(app)

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

    if file_method == "display":
        return jsonify(display_files(UPLOAD_FOLDER))

    elif file_method == "delete":
        return jsonify(delete_file(file_path, UPLOAD_FOLDER))

    elif file_method == "uploadFasta":
        dropzone_fasta(UPLOAD_FOLDER, request, current_user)
        return jsonify("success")

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
        if file_method == "checkFiles":
            return jsonify(find_blast_zip(current_user))
        elif file_method == "autoAnnotate":
            auto_annotate(UPLOAD_FOLDER, current_user)
            return jsonify("success")
        
    if request.method == "POST":
        if file_method == "downloadInput":
            return download_blast_input(current_user)

        elif file_method == "createInput":
            return jsonify(create_blast_input(UPLOAD_FOLDER, current_user))

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

        elif file_method == "drop":
            dropzone(UPLOAD_FOLDER, request)
            return jsonify("success")

        elif file_method == "deleteBlastResults":
            size1 = db.session.query(Blast_Results).count()
            print(size1)
            time.sleep(.2)
            size2 = db.session.query(Blast_Results).count()
            print(size2)
            if (size1 == size2):
                db.session.query(Blast_Results).delete()
                db.session.commit()
                return jsonify("success")
            else:
                return jsonify("fail")
        else:
            index = file_path.find(".json")
            name = file_path[0:index + 5]
            size = file_path[index + 5:]
            # file_data = db.session.query(Files).filter_by(name=name).first()
            # if file_data != None:
            #     db.session.delete(file_data)
            #     db.session.commit()
            file_data = Files(name=name,
                            date=file_method,
                            size=size,
                            complete=False)
            try:
                db.session.add(file_data)
                db.session.commit()
            except:
                return jsonify("already added")
            return jsonify("success")

@app.route('/phlash_api/annotations/<current_user>/<file_method>', methods=['GET', 'POST', 'PUT'])
def annotate_data(current_user, file_method):
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
        print("get")
        if file_method == "delete":
            db.session.query(Blast_Results).delete()
            db.session.commit()
            print("delete")
            return jsonify("success")
        if file_method == "blast":
            if (db.session.query(Blast_Results).first() is None):
                print("empty")
                parse_blast(UPLOAD_FOLDER)
                return jsonify("empty")
            else:
                return jsonify("not empty")
        else:
            print("CDS")
            return jsonify(get_dnamaster_data())
        print("done")

    if request.method == "POST":
        return get_genbank(UPLOAD_FOLDER, current_user, request)

    if request.method == "PUT":
        return jsonify(add_cds(request, UPLOAD_FOLDER, current_user))

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
        return jsonify(annotate_cds(request, cds_id, UPLOAD_FOLDER))

@app.route('/phlash_api/annotations/geneMap/<current_user>/', methods=['GET'])
def gene_map(current_user):
    """
    Builds and returns the gene map.
    """
    DATABASE = "sqlite:///{}".format(os.path.join(ROOT,
                                                  'users', current_user, f"{current_user}.db"))
    app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')
    if request.method == "GET":
        return jsonify(get_map(UPLOAD_FOLDER))

@app.route('/phlash_api/settings/<current_user>/<payload>/', methods=['PUT', 'GET'])
def settings(current_user, payload):
    """
    Updates default settings.
    """
    DATABASE = "sqlite:///{}".format(os.path.join(ROOT,
                                                  'users', current_user, f"{current_user}.db"))
    app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
    if request.method == "GET":
        if payload == "none":
            return jsonify(get_settings())
        else:
            return jsonify(update_settings(payload))

if __name__ == '__main__':
    app.run()