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
import random
import string
import datetime

response_object = {}
USERS = []
ROOT = os.path.dirname(os.path.abspath(__file__))
FASTA_EXTENSIONS = set(['.fasta', '.fna'])
GENBANK_EXTENSIONS = set(['.gb', '.gbk'])
GDATA_EXTENSIONS = set(['.gdata'])
LDATA_EXTENSIONS = set(['.ldata'])

# ------------------------------ MAIN FUNCTIONS ------------------------------
def check_phage_id(current_user, phage_id, app):
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
    stored_successfully = False
    while not stored_successfully:
        try:
            future_date = datetime.datetime.today() + datetime.timedelta(days=90)
            creation_date = str(datetime.datetime.today().strftime("%B %d, %Y"))
            deletion_date = str(future_date.strftime("%B %d, %Y"))
            random_id = ''.join(random.choice(string.ascii_letters) for _ in range(15))
            user = Users(user = current_user,
                        phage_id = phage_id,
                        creation_date = creation_date,
                        deletion_date = deletion_date,
                        id = random_id)
            db.session.add(user)
            db.session.commit()
            response_object["phage_id"] = user.phage_id
            response_object["id"] = user.id
            response_object["creation_date"] = user.creation_date
            response_object["deletion_date"] = user.deletion_date
            response_object["phage_page"] = "upload"
            stored_successfully = True
        except:
            stored_successfully = False
        
    handle_new_users(str(user.id), app)
    return response_object

def get_phage_data(email):
    """Gets all of the data associated with a user.

    Args:
        email:
            The email of a user.

    Returns:
        A dictionary containing all associated data or a message indicating if no data is found.
    """
    response_object = {}
    phage_id_list = []
    phage_creation_date_list = []
    phage_deletion_date_list = []
    id_list = []
    phage_view_id_list = []
    phage_view_email_list = []
    id_view_list = []
    phage_view_pages = []
    phages = db.session.query(Users).filter_by(user=email)
    if phages is not None:
        for phage in phages:
            if phage.creation_date == 'view':
                phage_view_id_list.append(phage.phage_id)
                phage_view_email_list.append(phage.deletion_date)
                id_view_list.append(phage.id)
                phage_view_pages.append('annotations')
            else:    
                phage_id_list.append(phage.phage_id)
                phage_creation_date_list.append(phage.creation_date)
                phage_deletion_date_list.append(phage.deletion_date)
                id_list.append(phage.id)
        response_object["phage_id_list"] = phage_id_list
        response_object["phage_creation_date_list"] = phage_creation_date_list
        response_object["phage_deletion_date_list"] = phage_deletion_date_list
        response_object["id_list"] = id_list
        response_object["phage_pages"] = find_current_page(id_list)
        response_object["phage_view_id_list"] = phage_view_id_list
        response_object["phage_view_email_list"] = phage_view_email_list
        response_object["id_view_list"] = id_view_list
        response_object["phage_view_pages"] = phage_view_pages
        return response_object
    return "empty"

def remove_old_users():
    """Removes 90 day old users.
    """
    critical_time = arrow.now().shift(days=-90)
    for user in Path(os.path.join(ROOT, 'users')).glob('*'):
        phage_id = str(user)[str(user).rfind('/') + 1:]
        user_time = arrow.get(user.stat().st_mtime)
        if phage_id != "Phlash.db" and user_time < critical_time:
            shutil.rmtree(user)
            db.session.query(Files).filter_by(phage_id=phage_id).delete()
            db.session.query(Annotations).filter_by(phage_id=phage_id).delete()
            db.session.query(Blast_Results).filter_by(phage_id=phage_id).delete()
            db.session.query(Gene_Calls).filter_by(phage_id=phage_id).delete()
            db.session.query(Tasks).filter_by(phage_id=phage_id).delete()
            db.session.query(Users).filter_by(id=phage_id).delete()
            db.session.commit()

def remove_user(phage_id, user):
    """Removes all data associated with a phage ID.

    Args:
        phage_id:
            The ID of the to be deleted phage.

    Returns:
        A success message.
    """
    if user.creation_date == 'view':
        db.session.delete(user)
        db.session.commit()
    else:
        for user in Path(os.path.join(ROOT, 'users')).glob('*'):
            if phage_id == str(user)[str(user).rfind('/') + 1:]:
                shutil.rmtree(user)
                db.session.query(Files).filter_by(phage_id=phage_id).delete()
                db.session.query(Annotations).filter_by(phage_id=phage_id).delete()
                db.session.query(Blast_Results).filter_by(phage_id=phage_id).delete()
                db.session.query(Gene_Calls).filter_by(phage_id=phage_id).delete()
                db.session.query(Tasks).filter_by(phage_id=phage_id).delete()
                db.session.query(Users).filter_by(id=phage_id).delete()
                db.session.commit()
                break
    return("success")

def get_users():
    """Gets existing users.
    """
    for dir in os.listdir(os.path.join(ROOT, 'users')):
        USERS.append(dir)

def find_current_page(phages):
    """Finds the page that a user is currently on per a phage ID.

    Args:
        phages:
            A list of phage IDs.

    Returns:
        A list of the current pages in the same order of the phage IDs in phages.
    """
    phage_pages = []
    for phage_id in phages:
        fasta = False
        for filename in os.listdir(os.path.join(ROOT, 'users', phage_id, 'uploads')):
            ext = os.path.splitext(filename)[1].lower()
            if ext in FASTA_EXTENSIONS:
                fasta = True
                if db.session.query(Blast_Results).filter_by(phage_id=phage_id).first() is None:
                    phage_pages.append('blast')
                else:
                    phage_pages.append('annotations')
                break
        if not fasta:
            phage_pages.append('upload')
    return phage_pages 

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
    setting = Settings(phage_id = phage_id,
                        back_left_range = 300,
                        forward_left_range = 100,
                        gap = 10,
                        overlap = 10,
                        opposite_gap = 50,
                        short = 200)
    db.session.add(setting)
    db.session.commit()
    response_object["id_status"] = "ID created. Please continue below."
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
    except FileExistsError:
        print("Directory \'" + directory + "\' already exists.")