from flask import *
import os
from models import *
from sqlalchemy import *
from sqlalchemy.orm import sessionmaker
from database_migration import *

ROOT = os.path.dirname(os.path.abspath(__file__))
DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', "Phlash.db"))
engine = create_engine(DATABASE)
Session = sessionmaker(bind=engine)
session = Session()
print_tables(os.path.join(ROOT, 'users', "Phlash.db"))
database_version = session.query(Database_Version).first()
with open(os.path.join(ROOT, "VERSION"), 'r') as version_num:
    curr_version = int(version_num.read())
    if not database_version:
        database_version = Database_Version(version = curr_version)
        session.add(database_version)       
    elif database_version.version != curr_version:
        try:
            print(eval('merge' + str(database_version.version) + '_' + str(curr_version) + '("' + os.path.join(ROOT, 'users', "Phlash.db") + '")'))
            database_version.version = curr_version
        except Exception as e:
            print(e)
            print("Write a method in database_migration called " + 'merge' + str(database_version.version) + '_' + str(curr_version) + '()' " to merge database changes.")
        print(database_version.version)
    session.commit()