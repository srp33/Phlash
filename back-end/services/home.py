import main
from pathlib import Path
from flask import *
from models import *
import os
import shutil
import arrow
from builtins import FileExistsError

response_object = {}
USERS = []
ROOT = os.path.dirname(os.path.abspath(__file__))
FASTA_EXTENSIONS = set(['.fasta', '.fna'])
GENBANK_EXTENSIONS = set(['.gb', '.gbk'])
GDATA_EXTENSIONS = set(['.gdata'])
LDATA_EXTENSIONS = set(['.ldata'])

def check_phage_id(phage_id, app):
    # remove 90 day old users
    remove_old_users()
    # get list of existing users
    get_users()
    if phage_id in USERS:
        handle_existing_users(phage_id)
    else:
        handle_new_users(phage_id, app)
    # get date that ID will be deleted
    get_deletion_date(phage_id)
    return response_object

def remove_old_users():
    critical_time = arrow.now().shift(days=-90)
    for user in Path(os.path.join(ROOT, 'users')).glob('*'):
        user_time = arrow.get(user.stat().st_mtime)
        if user_time < critical_time:
            shutil.rmtree(user)

def get_users():
    for dir in os.listdir(os.path.join(ROOT, 'users')):
        USERS.append(dir)

def handle_existing_users(phage_id):
    response_object['id_status'] = "ID already exists. If this is your ID, please continue. If not, enter a new one."

    # check if all required files for phage_id are uploaded
    required_files = ["fasta", "genbank", "gdata", "ldata"]
    for filename in os.listdir(os.path.join(ROOT, 'users', phage_id, 'uploads')):
        ext = os.path.splitext(filename)[1].lower()
        if ext in FASTA_EXTENSIONS:
            required_files.remove("fasta")
        elif ext in GENBANK_EXTENSIONS:
            required_files.remove("genbank")
        elif ext in GDATA_EXTENSIONS:
            required_files.remove("gdata")
        elif ext in LDATA_EXTENSIONS:
            required_files.remove("ldata")

    response_object['uploaded_all_files'] = True if len(required_files) == 0 else False

def handle_new_users(phage_id, app):
    create_directory(os.path.join(ROOT, 'users', phage_id))
    create_directory(os.path.join(ROOT, 'users', phage_id, 'uploads'))
    with app.app_context():
        db.drop_all()
        db.create_all()
    response_object["id_status"] = "ID created. Please continue."
    response_object["uploaded_all_files"] = False

def get_deletion_date(phage_id):
        critical_time = arrow.now().shift(days=-90)
        for user in Path(os.path.join(ROOT, 'users')).glob('*'):
            if phage_id in str(user):
                user_time = arrow.get(user.stat().st_mtime).shift(days=+90)
                response_object["delete_time"] = str(user_time)

def create_directory(directory):
    """
    Create specified directory.
    """
    try:
        os.mkdir(directory)
        print("Directory \'" + directory + "\' created.")
    except FileExistsError:
        print("Directory \'" + directory + "\' already exists.")