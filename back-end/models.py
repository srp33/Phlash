"""
Declaring models (tables) using Flask-SQLAlchemy for use in each phage's database. 
"""
from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()

class DNAMaster(db.Model):
    __tablename__ = "dnamaster"
    id = db.Column(db.Text,
                   nullable=False,
                   primary_key=True)
    start = db.Column(db.Integer,
                      nullable=False)
    stop = db.Column(db.Integer,
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


class GeneMark(db.Model):
    __tablename__ = "genemark"
    id = db.Column(db.Text,
                   nullable=False,
                   primary_key=True)
    start = db.Column(db.Integer,
                      nullable=False)
    stop = db.Column(db.Integer,
                     nullable=False)
    strand = db.Column(db.Text,
                       nullable=False)

class Gene_Calls(db.Model):
    __tablename__ = "gene_calls"
    id = db.Column(db.Text,
                   nullable=False,
                   primary_key=True)
    calls = db.Column(db.Text,
                      nullable=True)

class Blast_Results(db.Model):
    __tablename__ = "blast_results"
    id = db.Column(db.Text,
                   nullable=False,
                   primary_key=True)
    start = db.Column(db.Integer,
                      nullable=False)
    stop = db.Column(db.Integer,
                     nullable=False)
    strand = db.Column(db.Text,
                       nullable=False)
    results = db.Column(db.Text,
                       nullable=True)

class Settings(db.Model):
    __tablename__ = "settings"
    back_start_range = db.Column(db.Integer,
                                nullable=False,
                                primary_key=True)
    forward_start_range = db.Column(db.Integer,
                                    nullable=False)
    opposite_gap = db.Column(db.Integer,
                            nullable=False)
    gap = db.Column(db.Integer,
                    nullable=False)
    overlap = db.Column(db.Integer,
                        nullable=False)
    short = db.Column(db.Integer,
                        nullable=False)
