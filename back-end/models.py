# from flask import Flask
from flask_sqlalchemy import SQLAlchemy
# from sqlalchemy.orm.attributes import QueryableAttribute

db = SQLAlchemy()

# class Users(db.Model):
#    __tablename__ = "users"
#    username = db.Column(db.Text,
#                         nullable=False,
#                         primary_key=True)
   
#    def __init__(self, username):
#       self.username = username

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
