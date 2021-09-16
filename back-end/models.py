"""
Declaring models (tables) using Flask-SQLAlchemy for use in each phage's database. 
"""
from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()
class Users(db.Model):
    __tablename__="users"
    user = db.Column(db.Text,
                             nullable=False,
                             primary_key=True)
    phage_id = db.Column(db.Text,
                         nullable=False,
                         primary_key=True)
    creation_date = db.Column(db.Text,
                              nullable=False,
                              primary_key=True)
    deletion_date = db.Column(db.Text,
                              nullable=True)
    id = db.Column(db.Text,
                   primary_key=True,
                   nullable=False)
class Annotations(db.Model):
    __tablename__ = "annotations"
    phage_id = db.Column(db.Text,
                         nullable=False,
                         primary_key=True)
    id = db.Column(db.Text,
                   nullable=False,
                   primary_key=True)
    left = db.Column(db.Integer,
                      nullable=False)
    right = db.Column(db.Integer,
                     nullable=False)
    strand = db.Column(db.Text,
                       nullable=False)
    function = db.Column(db.Text,
                         nullable=False)
    status = db.Column(db.Text,
                       nullable=False)
    frame = db.Column(db.Integer,
                        nullable=False)
    notes = db.Column(db.Text,
                        nullable=True)

class Gene_Calls(db.Model):
    __tablename__ = "gene_calls"
    phage_id = db.Column(db.Text,
                         nullable=False,
                         primary_key=True)
    id = db.Column(db.Text,
                   nullable=False,
                   primary_key=True)
    calls = db.Column(db.Text,
                      nullable=True)

class Blast_Results(db.Model):
    __tablename__ = "blast_results"
    phage_id = db.Column(db.Text,
                         nullable=False,
                         primary_key=True)
    id = db.Column(db.Text,
                   nullable=False,
                   primary_key=True)
    left = db.Column(db.Integer,
                      nullable=False)
    right = db.Column(db.Integer,
                     nullable=False)
    strand = db.Column(db.Text,
                       nullable=False)
    results = db.Column(db.Text,
                       nullable=True)

class Settings(db.Model):
    __tablename__ = "settings"
    phage_id = db.Column(db.Text,
                         nullable=False,
                         primary_key=True)
    back_left_range = db.Column(db.Integer,
                                nullable=False,
                                primary_key=True)
    forward_left_range = db.Column(db.Integer,
                                    nullable=False)
    opposite_gap = db.Column(db.Integer,
                            nullable=False)
    gap = db.Column(db.Integer,
                    nullable=False)
    overlap = db.Column(db.Integer,
                        nullable=False)
    short = db.Column(db.Integer,
                      nullable=False)
    prodigal = db.Column(db.Boolean,
                         nullable=False)
    glimmer = db.Column(db.Boolean,
                         nullable=False)
    genemark = db.Column(db.Boolean,
                         nullable=False)
    aragorn = db.Column(db.Boolean,
                         nullable=False)
    phanotate = db.Column(db.Boolean,
                         nullable=False)

class Files(db.Model):
    __tablename__ = "files"
    phage_id = db.Column(db.Text,
                         nullable=False,
                         primary_key=True)
    name = db.Column(db.Text,
                     nullable=False,
                     primary_key=True)
    date = db.Column(db.Text,
                     nullable=False)
    size = db.Column(db.Text,
                     nullable=False)
    complete = db.Column(db.Boolean,
                         nullable=True)

class Tasks(db.Model):
    __tablename__ = "tasks"
    phage_id = db.Column(db.Text,
                         nullable=False,
                         primary_key=True)
    function = db.Column(db.Text,
                         nullable=False,
                         primary_key=True)
    arguments = db.Column(db.Text,
                          nullable=True)
    complete = db.Column(db.Boolean,
                         nullable=False)
    result = db.Column(db.Text,
                       nullable=True)
    time = db.Column(db.Text,
                     nullable=False)

class MetaData(db.Model):
    __tablename__= "metadata"
    version = db.Column(db.Integer,
                        primary_key=True)