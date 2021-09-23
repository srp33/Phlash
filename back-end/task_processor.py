import time
from flask import *
import os
from models import *
from sqlalchemy import *
from sqlalchemy.orm import sessionmaker
from background_tasks import *

# Configuration
ROOT = os.path.dirname(os.path.abspath(__file__))
DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', "Phlash.db"))
engine = create_engine(DATABASE)
Session = sessionmaker(bind=engine)
session = Session()

def get_task():
    return session.query(Tasks).filter_by(complete=False).filter_by(result="waiting").order_by(Tasks.time).first()

def run_task(task):
    task.result = "executing"
    session.commit()
    args = list(task.arguments.split(" "))
    if task.function == "parse_blast":
        try:
            task.result = parse_blast(args[0], args[1], session)
        except:
            task.result = "error"
            return None   
    elif task.function == "blast_input":
        try:
            task.result = create_blast_input(args[0], args[1])
        except:
            task.result = "error"
            return None
    elif task.function == "auto_annotate":
        setting = session.query(Settings).filter_by(phage_id=args[1]).order_by(Settings.back_left_range).first()
        try:
            run_genemark(args[0], args[1])
        except:
            task.result = "error"
            return None
        if setting.glimmer:
            try:
                run_glimmer(args[0], args[1])
            except:
                task.result = "error"
                return None
        if setting.aragorn:
            try:
                run_aragorn(args[0], args[1])
            except:
                task.result = "error"
                return None
        if setting.phanotate:
            try:
                run_phanotate(args[0], args[1])
            except:
                task.result = "error"
                return None
        if setting.prodigal:
            try:
                run_prodigal(args[0], args[1])
            except:
                task.result = "error"
                return None
        try:
            task.result = auto_annotate(args[0], args[1], session)
        except:
            task.result = "error"
            return None

def handle_unfinished_tasks():
    old_tasks = session.query(Tasks).filter_by(complete=False).filter_by(result="executing").order_by(Tasks.time)
    if old_tasks:
        for old_task in old_tasks:
            args = list(old_task.arguments.split(" "))
            if old_task.function == "auto_annotate":
                session.query(Annotations).filter_by(phage_id=args[1]).delete()
                session.query(Gene_Calls).filter_by(phage_id=args[1]).delete()
                for filename in os.listdir(args[0]):
                    if not filename.endswith(".fasta") and not filename.endswith(".json"):
                        os.remove(os.path.join(args[0], filename))
                old_task.result = "waiting"
                session.commit()
            elif old_task.function == "blast_input":
                try:
                    USER_FOLDER = args[0][:-8]
                    for file in os.listdir(USER_FOLDER):
                        if file.endswith(".zip") or file.endswith(".fasta"):
                            os.remove(os.path.join(USER_FOLDER, file))
                except FileNotFoundError:
                    print("zip file not found")
                old_task.result = "waiting"
                session.commit()
            else:
                session.query(Blast_Results).filter_by(phage_id=args[0]).delete()
                session.query(Files).filter_by(phage_id=args[0]).delete()
                for filename in os.listdir(args[1]):
                    if filename.endswith('.json'):
                        os.remove(os.path.join(args[1], filename))
                old_task.result = "waiting"
                session.commit()

handle_unfinished_tasks()

while True:
    task = get_task()
    if task:
        run_task(task)
        task.complete = True
        session.commit()
    else:
        time.sleep(2)