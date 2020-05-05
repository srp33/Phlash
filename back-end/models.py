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
    start_options = db.Column(db.Text)


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


# Model class of files uploaded
class Files(db.Model):
    __tablename__ = "files"
    id = db.Column(db.Text,
                   primary_key=True)
    name = db.Column(db.Text)
    date = db.Column(db.DateTime)

    def __init__(self, id, name, date):
        self.id = id
        self.name = name
        self.date = date
