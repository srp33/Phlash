from flask import *
import os

from sqlalchemy.sql.ddl import DropConstraint
from models import *
from sqlalchemy import *
from sqlalchemy.orm import sessionmaker
from sqlalchemy_utils import database_exists, create_database, drop_database
from database_migration import *

ROOT = os.path.dirname(os.path.abspath(__file__))
DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', "Phlash.db"))
engine = create_engine(DATABASE)
Session = sessionmaker(bind=engine)
session = Session()
print_tables(os.path.join(ROOT, 'users', "Phlash.db"))

version_file_path = os.path.join(ROOT, "VERSION")

if not os.path.exists(version_file_path):
    raise Exception("No VERSION file exists.")

with open(version_file_path, 'r') as version_num:
    version = version_num.read()

    if version == "":
        raise EOFError("The file VERSION needs an integer value indicating the version number to continue.")

    if not database_exists(engine.url):
        app = Flask(__name__)
        app.config.from_object(__name__)
        db.init_app(app)
        app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE

        with app.app_context():
            db.create_all()

    curr_version = int(version)
    database_version = session.query(Database_Version).first()

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
    session.commit()
