from werkzeug.utils import secure_filename
from Bio import SeqIO, Seq
from flask_cors import CORS
from flask import *
from models import *
import annotate
import os, datetime, uuid
import pandas as pd

# configuration
ROOT = os.path.dirname(os.path.abspath(__file__))
FASTA_EXTENSIONS = set(['fasta', 'fna'])
GENBANK_EXTENSIONS = set(['gb', 'gbk'])
GDATA_EXTENSION = set(['gdata'])
LDATA_EXTENSION = set(['ldata'])
USERS = []
CURRENT_USER = ""

# instantiate the app
app = Flask(__name__)
app.config.from_object(__name__)
db.init_app(app)

# enable CORS
CORS(app, resources={r'/*': {'origins': '*'}})

# routers ------------------------------------------------------------------
@app.route('/api/home/<username>', methods=['POST'])
def check_username(username):
   """
   """
   response_object = {'status': 'success'}

   if request.method == "POST":
      USERS = []
      for dir in os.listdir(os.path.join(ROOT, 'users')):
         USERS.append(dir)
      print(USERS)

      if username in USERS:
         DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', username, f"{username}.db"))
         app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
         response_object['message'] = "ID already exists. If this is your ID, please continue. If not, enter a new one."
      else:
         create_directory(os.path.join(ROOT, 'users', username))
         create_directory(os.path.join(ROOT, 'users', username, 'uploads'))
         DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', username, f"{username}.db"))
         app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
         with app.app_context():
            db.drop_all()
            db.create_all()
         response_object["message"] = "ID created. Please continue."

   return jsonify(response_object)


@app.route('/api/upload/<current_user>', methods=['POST'])
def upload_file(current_user):
   """
   User uploads all required files.
   POST method uploads file.
   """
   UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')
   DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', current_user, f"{current_user}.db"))
   app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE

   response_object = {'status': 'success'}

   if request.method == "POST":
      if 'file' not in request.files:
         response_object["status"]: "'file' not in request.files"
         print("in fail")
      else:
         print("in success")
         file = request.files['file']
         fileType = request.form['fileType']
         if file:
            file_name = secure_filename(file.filename)
            if (fileType == "fasta" and allowed_file(file_name, FASTA_EXTENSIONS)) or \
               (fileType == "genbank" and allowed_file(file_name, GENBANK_EXTENSIONS)) or \
               (fileType == "gdata" and allowed_file(file_name, GDATA_EXTENSION)) or \
               (fileType == "ldata" and allowed_file(file_name, LDATA_EXTENSION)):

               file_ext = file_name.rsplit('.', 1)[1].lower()
               for existing_file in os.listdir(UPLOAD_FOLDER):
                  if existing_file.endswith(f".{file_ext}"):
                     os.remove(os.path.join(UPLOAD_FOLDER, existing_file))

               file.save(os.path.join(UPLOAD_FOLDER, file_name))
               response_object["uploaded"] = file_name
               print(' * uploaded', file_name)

               if file_ext in ['gb', 'gbk']:
                  print(UPLOAD_FOLDER)
                  genbank_file = get_file("GenBank", UPLOAD_FOLDER)
                  print(genbank_file)
                  annotate.parse_genbank(genbank_file)
               elif file_ext == "ldata":
                  genemark_ldata_file = get_file("GeneMark_ldata", UPLOAD_FOLDER)
                  annotate.parse_genemark_ldata(genemark_ldata_file)
            else:
               response_object["not_allowed"] = file.filename
         else:
            response_object["status"] = "error"

   return jsonify(response_object)


@app.route('/api/dnamaster/<current_user>', methods=['GET', 'POST'])
def view_dnamaster(current_user):
   """
   User can view DNA Master data in database.
   GET method displays all data in table format.
   POST method allows user to add a new CDS to the data.
   """
   DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', current_user, f"{current_user}.db"))
   app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE

   response_object = {'status': 'success'}

   if request.method == "GET":
      dnamaster = []
      for cds in db.session.query(DNAMaster).order_by(DNAMaster.start):
         dnamaster.append({'id': cds.id,
                        'start': cds.start,
                        'stop': cds.stop,
                        'strand': cds.strand})
      response_object['dnamaster'] = dnamaster

   if request.method == "POST":
      post_data = request.get_json()
      cds = DNAMaster(id=post_data.get('id'),
                     start=int(post_data.get('start')),
                     stop=int(post_data.get('stop')),
                     strand=post_data.get('strand'),
                     function="None",
                     status="None")
      start_exists = DNAMaster.query.filter_by(start=post_data.get('start')).first()
      stop_exists = DNAMaster.query.filter_by(stop=post_data.get('stop')).first()
      id_exists = DNAMaster.query.filter_by(id=post_data.get('id')).first()
      if cds.start > cds.stop or cds.start == cds.stop:
         response_object['message'] = 'Start is not less than start. CDS not added.'
      elif start_exists and stop_exists:
         if start_exists.id == stop_exists.id:
            response_object['message'] = 'Start and stop already exists. CDS not added.'
      elif id_exists:
         if id_exists.id == cds.id:
            response_object['message'] = 'ID already exists. CDS not added.'
      else:
         db.session.add(cds)
         db.session.commit()
         response_object['message'] = 'CDS added!'

   return jsonify(response_object)


@app.route('/api/dnamaster/<current_user>/<cds_id>', methods=['PUT', 'DELETE'])
def single_cds(current_user, cds_id):
   """
   User can update or delete a CDS from the DNA Master table in database.
   PUT method updates a CDS.
   DELETE method deletes a CDS.
   """
   DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', current_user, f"{current_user}.db"))
   app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
   response_object = {'status': 'success'}

   if request.method == "PUT":
      put_data = request.get_json()
      cds = DNAMaster.query.filter_by(id=put_data.get('id')).first()
      if cds:
         cds.start = put_data.get('start')
         cds.stop = put_data.get('stop')
         cds.strand = put_data.get('strand')
         db.session.commit()
         response_object['message'] = 'CDS updated!'
      else:
         response_object['message'] = 'Error: CDS could not get updated.'

   if request.method == 'DELETE':
      if DNAMaster.query.filter_by(id=cds_id).first():
         DNAMaster.query.filter_by(id=cds_id).delete()
         db.session.commit()
         response_object['message'] = 'CDS removed!'
      else:
         response_object['message'] = 'Error: CDS could not be deleted.'

   return jsonify(response_object)


@app.route('/api/blast/<current_user>/<file_method>', methods=['POST'])
def upload_blast_file(current_user, file_method):
   """
   POST method uploads file.
   GET method adds data in file to database.
   """
   DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', current_user, f"{current_user}.db"))
   app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
   UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')
   response_object = {'status': 'success'}

   if request.method == "POST":
      if file_method == "download":
         print("Starting comparisons")
         annotate.compare()
         fasta_file = get_file("Fasta", UPLOAD_FOLDER)
         genemark_gdata_file = get_file("GeneMark_gdata", UPLOAD_FOLDER)
         output = annotate.create_fasta(fasta_file, genemark_gdata_file)

         with open(os.path.join(ROOT, 'users', current_user, f"{current_user}_blast.fasta"), "w") as out_handle:
            out_handle.write(output)

         f = open(os.path.join(ROOT, 'users', current_user, f"{current_user}_blast.fasta"), "r")
         return f.read()
      elif file_method == "upload":
         if 'file' not in request.files:
            response_object["status"]: "'file' not in request.files"
            print("in fail")
         else:
            print("in success")
            file = request.files['file']
            if file:
               file_name = secure_filename(file.filename)
               if allowed_file(file_name, set(['json'])):

                  file_ext = file_name.rsplit('.', 1)[1].lower()
                  for existing_file in os.listdir(UPLOAD_FOLDER):
                     if existing_file.endswith(f".{file_ext}"):
                        os.remove(os.path.join(UPLOAD_FOLDER, existing_file))

                  file.save(os.path.join(UPLOAD_FOLDER, file_name))
                  response_object["uploaded"] = file_name
                  print(' * uploaded', file_name)
               else:
                  response_object["not_allowed"] = file.filename
            else:
               response_object["status"] = "error"

   return jsonify(response_object)


@app.route('/api/annotations/<current_user>', methods=['GET', 'POST'])
def annotate_data(current_user):
   """
   Compares DNA Master's predictions against GeneMark's.
   GET method shows all the DNA Master predictions with a status and action item for each.
   """
   DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', current_user, f"{current_user}.db"))
   app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
   UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')
   response_object = {'status': 'success'}

   if request.method == "GET":
      dnamaster = []
      for cds in db.session.query(DNAMaster).order_by(DNAMaster.start):
         dnamaster.append({'id': cds.id,
                           'start': cds.start,
                           'stop': cds.stop,
                           'strand': cds.strand,
                           'function': cds.function,
                           'status': cds.status})
      response_object['dnamaster'] = dnamaster

   if request.method == "POST":
      # -------downloading GENBANK file----------
      gb_file = get_file("GenBank", UPLOAD_FOLDER)
      fasta_file = get_file("Fasta", UPLOAD_FOLDER)
      out_file = annotate.modify_genbank(gb_file, fasta_file)
      f = open(out_file, "r")
      return f.read()

   return jsonify(response_object)




@app.route('/api/annotations/cds/<current_user>/<cds_id>', methods=['GET', 'POST', 'PUT', 'DELETE'])
def cds_annotation(current_user, cds_id):
   """
   Annotation information for each CDS. 
   GET method gets cds, its start options, blast results, and graph data.
   PUT method updates the start position and function if the user chooses to do so.
   """
   DATABASE = "sqlite:///{}".format(os.path.join(ROOT, 'users', current_user, f"{current_user}.db"))
   app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
   UPLOAD_FOLDER = os.path.join(ROOT, 'users', current_user, 'uploads')
   response_object = {'status': 'success'}

   if request.method == "GET":
      cds = DNAMaster.query.filter_by(id=cds_id).first()
      response_object['cds'] = {'id': cds.id,
                              'start': cds.start,
                              'stop': cds.stop,
                              'strand': cds.strand,
                              'function': cds.function,
                              'status': cds.status}

      # fasta_file = get_file("Fasta", UPLOAD_FOLDER)
      starts = [int(start) for start in cds.start_options.split(",")]
      print(starts)
      response_object['start_options'] = starts


      blast_file = get_file("Blast", UPLOAD_FOLDER)
      E_VALUE_THRESH = 1e-7
      print("parsing blast")
      blast_results = annotate.parse_blast_multiple(blast_file, cds.id, E_VALUE_THRESH)
      print("writing results to post response")
      response_object['blast'] = blast_results

      genemark_gdata_file = get_file("GeneMark_gdata", UPLOAD_FOLDER)
      gdata_df = pd.read_csv(genemark_gdata_file, sep='\t', skiprows=16)
      gdata_df.columns = ['Base', '1', '2', '3', '4', '5', '6']
      gdata_df = gdata_df[gdata_df.Base.isin(range(min(starts) - 100, cds.stop + 100))]
      response_object['x_data'] = gdata_df["Base"].to_list()
      response_object['y_data_1'] = gdata_df["1"].to_list()
      response_object['y_data_2'] = gdata_df["2"].to_list()
      response_object['y_data_3'] = gdata_df["3"].to_list()
      response_object['y_data_4'] = gdata_df["4"].to_list()
      response_object['y_data_5'] = gdata_df["5"].to_list()
      response_object['y_data_6'] = gdata_df["6"].to_list()

   if request.method == "PUT":
      put_data = request.get_json()
      cds = DNAMaster.query.filter_by(id=cds_id).first()
      if cds:
         cds.start = put_data.get('start')
         cds.function = put_data.get('function')
         cds.status = put_data.get('status')
         db.session.commit()
         response_object['message'] = 'CDS updated!'
      else:
         response_object['message'] = 'CDS did not update.'

   if request.method == "DELETE":
      if DNAMaster.query.filter_by(id=cds_id).first():
         DNAMaster.query.filter_by(id=cds_id).delete()
         db.session.commit()
         response_object['message'] = 'CDS removed!'
      else:
         response_object['message'] = 'Error: CDS could not be deleted.'

   return jsonify(response_object)


if __name__ == '__main__':
    app.run()


# ---------- HELPER FUNCTIONS ----------

# check file extension. only allow specific ones
def allowed_file(filename, allowed_extensions):
   return '.' in filename and filename.rsplit('.', 1)[1].lower() in allowed_extensions

# create directory
def create_directory(directory):
   try:
      os.mkdir(directory)
      print("Directory \'" + directory + "\' created.")
   except FileExistsError:
      print("Directory \'" + directory + "\' already exists.")

# Get necessary file from uploads dir depending on preference
def get_file(preference, upload_directory):
   for filename in os.listdir(upload_directory):
      if preference == "GenBank":
         if filename.endswith(".gb") or filename.endswith(".genbank"):
            return os.path.join(upload_directory, filename)
      elif preference == "GeneMark_ldata":
         if filename.endswith(".ldata"):
            return os.path.join(upload_directory, filename)
      elif preference == "GeneMark_gdata":
         if filename.endswith(".gdata"):
            return os.path.join(upload_directory, filename)
      elif preference == "Fasta":
         if filename.endswith(".fasta"):
            return os.path.join(upload_directory, filename)
      elif preference == "Blast":
         if filename.endswith(".json"):
            return os.path.join(upload_directory, filename)
      else:
         return("Couldn't find file.")

