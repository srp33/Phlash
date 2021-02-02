"""Contains the functions for the Home page.

Attributes:
    response_object:
        The dictionary that the main functions return.
    USERS:
        A list of users.
    ROOT:
        The root directory.
    FASTA_EXTENSIONS:
        The allowed fasta extensions.
    GENBANK_EXTENSIONS:
        The allowed genbank extensions.
    GDATA_EXTENSIONS:
        The allowed gdata extensions.
    LDATA_EXTENSIONS:
        The allowed ldata extensions.
"""
import main
from pathlib import Path
from flask import *
from models import *
import os
import shutil
import arrow
from builtins import FileExistsError
import models

response_object = {}
USERS = []
ROOT = os.path.dirname(os.path.abspath(__file__))
FASTA_EXTENSIONS = set(['.fasta', '.fna'])
GENBANK_EXTENSIONS = set(['.gb', '.gbk'])
GDATA_EXTENSIONS = set(['.gdata'])
LDATA_EXTENSIONS = set(['.ldata'])

# ------------------------------ MAIN FUNCTIONS ------------------------------
def check_phage_id(phage_id, app):
    """Checks to see if phage ID exists.

    Args:
        phage_id:
            The ID being logged in or registered.
        app:
            The application.

    Returns:
        A dictionary containing success message.
    """
    remove_old_users() # remove 90 day old users
    get_users()
    if phage_id in USERS:
        handle_existing_users(phage_id)
    else:
        handle_new_users(phage_id, app)
    get_deletion_date(phage_id)
    return response_object

def remove_old_users():
    """Removes 90 day old users.
    """
    critical_time = arrow.now().shift(days=-1)
    for user in Path(os.path.join(ROOT, 'users')).glob('*'):
        user_time = arrow.get(user.stat().st_mtime)
        if user_time < critical_time:
            shutil.rmtree(user)

def get_users():
    """Gets existing users.
    """
    for dir in os.listdir(os.path.join(ROOT, 'users')):
        USERS.append(dir)

def handle_existing_users(phage_id):
    """Checks to see if required files are uploaded.

    Args:
        phage_id:
            The ID the existing user.
    """
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
    response_object['blast_complete'] = False if db.session.query(Blast_Results).first() is None else True

def handle_new_users(phage_id, app):
    """Creates a new user and the associated directories.

    Args:
        phage_id:
            The new ID to be created.
        app:
            The application.
    """
    create_directory(os.path.join(ROOT, 'users', phage_id))
    create_directory(os.path.join(ROOT, 'users', phage_id, 'uploads'))
    with app.app_context():
        db.drop_all()
        db.create_all()
    setting = Settings(back_start_range = 300,
                        forward_start_range = 100,
                        gap = 10,
                        overlap = 10,
                        opposite_gap = 50,
                        short = 200)
    db.session.add(setting)
    db.session.commit()
    response_object["id_status"] = "ID created. Please continue."
    response_object["uploaded_all_files"] = False

def get_deletion_date(phage_id):
    """Gets the date that the phage_id will be removed.

    Args:
        phage_id:
            The ID of the new phage.
    """
    for user in Path(os.path.join(ROOT, 'users')).glob('*'):
        if phage_id in str(user):
            user_time = arrow.get(user.stat().st_mtime).shift(days=+90)
            response_object["delete_time"] = str(user_time)

def create_directory(directory):
    """Creates specified directory.

    Args:
        directory:
            The directory path to be created.
    """
    try:
        os.mkdir(directory)
        print("Directory \'" + directory + "\' created.")
    except FileExistsError:
        print("Directory \'" + directory + "\' already exists.")