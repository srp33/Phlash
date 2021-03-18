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
    db.create_all()

def run_tasks():
    app.app_context().push()
    task = db.session.query(Tasks).filter_by(complete=False).order_by(Tasks.time).first()
    if task is not None:
        if task.function == "parse_blast":
            args = list(task.arguments.split(" "))
            task.result = parse_blast(args[0], args[1])
            task.complete = True
            db.session.commit()
        elif task.function == "auto_annotate":
            args = list(task.arguments.split(" "))
            auto_annotate(args[0], args[1])
            task.complete = True
            db.session.commit()
    t = threading.Timer(5.0, run_tasks)
    t.daemon = True
    t.start()
run_tasks()

# routers ------------------------------------------------------------------
@app.route('/phlash_api/test', methods=['GET'])
def test():
    return jsonify("Hello, world!")

@app.route('/phlash_api/home/<phage_id>', methods=['POST'])
def check_phage(phage_id):
    """
    API endpoint for '/home/:phage_id'.
    POST method removes users that have existed for more than 90 days,
                creates a new user if it doesn't exist,
                else gets informations for existing user.
    """
    return jsonify(home.check_phage_id(phage_id, app))

@app.route('/phlash_api/check_upload/<current_user>', methods=['GET'])
def check_upload(current_user):
    """
    API endpoint for '/upload/:current_user'.
    POST method uploads files accordingly and removes files if necessary.
    """
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')

    return jsonify(check_uploaded_files(UPLOAD_FOLDER))

@app.route('/phlash_api/upload/<current_user>/<file_method>/<file_path>', methods=['POST'])
def upload(current_user, file_method, file_path):
    """
    API endpoint for '/upload/:current_user'.
    POST method uploads files accordingly and removes files if necessary.
    """
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')

    if file_method == "display":
        return jsonify(display_files(UPLOAD_FOLDER))

    elif file_method == "delete":
        return jsonify(delete_file(current_user, file_path, UPLOAD_FOLDER))

    elif file_method == "uploadFasta":
        dropzone_fasta(UPLOAD_FOLDER, request, current_user)
        return jsonify("success")

@app.route('/phlash_api/blast/<current_user>/<file_method>/<file_path>', methods=['GET', 'POST'])
def blast(current_user, file_method, file_path):
    """
    API endpoint for '/blast/:current_user/:file_method/:file_path'.
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
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')

    if request.method == "GET":
        if file_method == "checkFiles":
            return jsonify(find_blast_zip(current_user))

        elif file_method == "autoAnnotate":
            args = UPLOAD_FOLDER + " " + current_user
            task = Tasks(phage_id=current_user,
                        function="auto_annotate",
                        arguments=args,
                        complete=False,
                        time=datetime.now())
            try:
                db.session.add(task)
                db.session.commit()
            except:
                return jsonify("Error in adding task to queue")
            # auto_annotate(UPLOAD_FOLDER, current_user)
            return jsonify("success")
        
    if request.method == "POST":
        if file_method == "downloadInput":
            return download_blast_input(current_user)

        elif file_method == "createInput":
            return jsonify(create_blast_input(UPLOAD_FOLDER, current_user))

        elif file_method == "displayOutput":
            return jsonify(get_blast_output_names(current_user, UPLOAD_FOLDER, file_path))

        elif file_method == "deleteOutput":
            return jsonify(delete_blast_output(current_user, UPLOAD_FOLDER, file_path))

        elif file_method == "numFiles":
            return jsonify(get_num_blast_files(current_user))

        elif file_method == "drop":
            dropzone(current_user, UPLOAD_FOLDER, request)
            return jsonify("success")

        elif file_method == "deleteBlastResults":
            if db.session.query(Tasks).filter_by(phage_id=current_user).filter_by(function="parse_blast").first() is None:
                db.session.query(Blast_Results).filter_by(phage_id=current_user).delete()
                db.session.query(Files).filter_by(phage_id=current_user).delete()
                db.session.commit()
                for filename in os.listdir(UPLOAD_FOLDER):
                    if filename.endswith('.json'):
                        os.remove(os.path.join(UPLOAD_FOLDER, filename))
                return jsonify("success")
            else:
                return jsonify("fail")

        else:
            index = file_path.find(".json")
            name = secure_filename(file_path[0:index + 5])
            size = file_path[index + 5:]
            # file_data = db.session.query(Files).filter_by(name=name).first()
            # if file_data != None:
            #     db.session.delete(file_data)
            #     db.session.commit()
            file_data = Files(phage_id=current_user,
                            name=name,
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
    Parses blast data and returns annotation data.
    GET method shows all the annotation predictions with a status and action item for each.
    PUT method adds a CDS
    """
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')
    response_object = {'status': 'success'}

    if request.method == "GET":
        if file_method == "check":
            curr_tasks = db.session.query(Tasks).filter_by(complete=False).order_by(Tasks.time)
            counter = 0
            for curr_task in curr_tasks:
                if curr_task.phage_id == current_user:
                    break
                counter += 1
            task = db.session.query(Tasks).filter_by(phage_id=current_user).filter_by(function="parse_blast").first()
            if task is not None and task.complete:
                result = task.result
                db.session.delete(task)
                db.session.commit()
                return jsonify(result)
            elif task is None:
                return jsonify("complete")
            else:
                return jsonify(str(counter))
        elif file_method == "delete":
            db.session.query(Blast_Results).filter_by(phage_id=current_user).delete()
            db.session.commit()
            return jsonify("success")
        elif file_method == "blast":
            print(db.session.query(Blast_Results).filter_by(phage_id=current_user).first())
            print(db.session.query(Tasks).filter_by(phage_id=current_user).filter_by(function="parse_blast").first())
            if db.session.query(Blast_Results).filter_by(phage_id=current_user).first() is None and db.session.query(Tasks).filter_by(phage_id=current_user).filter_by(function="parse_blast").filter_by(complete=False).first() is None:
                args = current_user + " " + UPLOAD_FOLDER
                task = Tasks(phage_id=current_user,
                                function="parse_blast",
                                arguments=args,
                                complete=False,
                                time=datetime.now())
                try:
                    db.session.add(task)
                    db.session.commit()
                except:
                    return jsonify("Error in adding task to queue")
                return jsonify("empty")
            if db.session.query(Tasks).filter_by(phage_id=current_user).filter_by(function="parse_blast").first() is not None:
                return jsonify("empty")
            else:
                return jsonify("not empty")
        else:
            return jsonify(get_annotations_data(current_user))

    if request.method == "PUT":
        return jsonify(add_cds(request, UPLOAD_FOLDER, current_user))

@app.route('/phlash_api/genbank/<current_user>/<file_method>', methods=['GET', 'POST', 'PUT'])
def create_genbank(current_user, file_method):
    """
    Creates the Genbank file from data stored in database.
    POST method creates and returns the genbank file.
    """
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')
    response_object = {'status': 'success'}

    if request.method == "POST":
        return get_genbank(UPLOAD_FOLDER, current_user, request)

@app.route('/phlash_api/annotations/cds/<current_user>/<cds_id>', methods=['GET', 'POST', 'PUT', 'DELETE'])
def cds_annotation(current_user, cds_id):
    """
    Annotation information for each CDS.
    GET method gets cds, its left options, blast results, and graph data.
    PUT method updates the left position and function if the user chooses to do so.
    """
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')
    response_object = {'status': 'success'}

    if request.method == "GET":
        return jsonify(get_cds_data(current_user, UPLOAD_FOLDER, cds_id))
                
    if request.method == "PUT":
        return jsonify(annotate_cds(current_user, request, cds_id, UPLOAD_FOLDER))

@app.route('/phlash_api/annotations/geneMap/<current_user>/', methods=['GET'])
def gene_map(current_user):
    """
    Builds and returns the gene map.
    """
    UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')
    if request.method == "GET":
        return jsonify(get_map(current_user, UPLOAD_FOLDER))

@app.route('/phlash_api/settings/<current_user>/<payload>/', methods=['PUT', 'GET'])
def settings(current_user, payload):
    """
    Updates default settings.
    """
    if request.method == "GET":
        if payload == "none":
            return jsonify(get_settings(current_user))
        else:
            return jsonify(update_settings(current_user, payload))

if __name__ == '__main__':
    app.run()