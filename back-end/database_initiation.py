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
with open(os.path.join(ROOT, "VERSION"), 'r+') as version_num:
    version = version_num.read()
    if version == "":
        version_num.write("0")
        version = "0"
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
        print(database_version.version)
    session.commit()