from werkzeug.utils import secure_filename
from Bio import SeqIO
from flask_cors import CORS
from flask import *
from models import *
import annotate
import graph
import os, datetime, uuid
import pandas as pd
import blast

# configuration
DEBUG = True
ROOT = os.path.dirname(os.path.abspath(__file__))
UPLOAD_FOLDER = os.path.join(ROOT, 'uploads')
DATABASE = "sqlite:///{}".format(os.path.join(ROOT, "database.db"))
ALLOWED_EXTENSIONS = set(['gb', 'gbk', 'fasta', 'fna', 'faa', 'ffn', 'frn', 'gdata', 'ldata'])

# instantiate the app
app = Flask(__name__)
app.config.from_object(__name__)
app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
db.init_app(app)
with app.app_context():
    db.drop_all()
    db.create_all()
# db = SQLAlchemy(app)

# enable CORS
CORS(app, resources={r'/*': {'origins': '*'}})


# foo router check --------------------------------------------------------------
@app.route('/ping', methods=['GET'])
def ping_pong():
    return jsonify('pong!')

# real routers ------------------------------------------------------------------
@app.route('/api/upload', methods=['GET', 'POST'])
def upload_file():
   """
   User uploads DNA Master file (genbank). Add data from file to the database.
   POST method uploads file.
   GET method adds data in file to database.
   """
   create_directory(UPLOAD_FOLDER)
   response_object = {'status': 'success'}
   response_object["uploaded"] = []
   response_object["not_allowed"] = []
   response_object["required"] = ["fasta", "genbank", "gdata", "ldata"]

   # # generate GUUID
   # file_id = str(uuid.uuid4())
   # file.save(os.path.join(app.config['UPLOAD_FOLDER'], file_name))
   # save file info into DB
   # file = Files(id=file_id, name=file_name, \
   #              date=datetime.datetime.now(datetime.timezone.utc))
   # db.session.add(file)
   #    db.session.commit()
   # response_object['message'] = f"{file_name} uploaded successfully."

   if request.method == "POST":
      if request.files.getlist('files') is None: 
         response_object["status"]: "request.files.getlist('files') is None"
      else:
         for file in request.files.getlist('files'):
            if file and allowed_file(file.filename):
               print(allowed_file(file.filename))
               file_name = secure_filename(file.filename)
               file.save(os.path.join(app.config['UPLOAD_FOLDER'], file_name))
               print(' * file uploaded', file_name)
               response_object["uploaded"].append(f"{file_name}")
               file_ext = file_name.rsplit('.', 1)[1].lower()
               print(file_ext)
               if file_ext in ['fasta', 'fna', 'faa', 'ffn', 'frn']:
                  response_object["required"].remove("fasta")
               elif file_ext in ['gb', 'gbk']:
                  response_object["required"].remove("genbank")
               elif file_ext == "gdata":
                  response_object["required"].remove("gdata")
               elif file_ext == "ldata":
                  response_object["required"].remove("ldata")
            elif file and allowed_file(file.filename) == False:
               response_object["not_allowed"].append(file.filename)
            else:
               response_object["status"] = "error"

   # if request.method == 'GET':
   #    genbank_file = get_file("GenBank")
   #    annotate.parse_genbank(genbank_file)
   #    response_object['message'] = "Data added to database."

   return jsonify(response_object)


@app.route('/database', methods=['GET', 'POST'])
def view_database():
    """
    User can view DNA Master data in database.
    GET method displays all data in table format.
    POST method allows user to add a new CDS to the data.
    """
    response_object = {'status': 'success'}

    if request.method == "GET":
        dnamaster = []
        for cds in db.session.query(DNAMaster).order_by(DNAMaster.start):
            dnamaster.append({'id': cds.id,
                             'start': cds.start,
                             'stop': cds.stop,
                             'strand': cds.strand})
                            #  'status': cds.status})
        response_object['dnamaster'] = dnamaster

    if request.method == "POST":
        post_data = request.get_json()
        cds = DNAMaster(id=post_data.get('id'),
                        start=int(post_data.get('start')),
                        stop=int(post_data.get('stop')),
                        strand=post_data.get('strand'),
                        function="None",
                        status="None")
        exists = DNAMaster.query.filter_by(id=post_data.get('id')).first()
        if not exists:
            db.session.add(cds)
            db.session.commit()
            response_object['message'] = 'CDS added!'
        else:
            response_object['message'] = 'CDS already exists'
    
    return jsonify(response_object)


@app.route('/database/<cds_id>', methods=['PUT', 'DELETE'])
def single_cds(cds_id):
    """
    User can update or delete a CDS from the DNA Master table in database.
    PUT method updates a CDS.
    DELETE method deletes a CDS.
    """
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
            response_object['message'] = 'CDS did not updated.'

    if request.method == 'DELETE':
        if DNAMaster.query.filter_by(id=cds_id).first():
            DNAMaster.query.filter_by(id=cds_id).delete()
            db.session.commit()
            response_object['message'] = 'CDS removed!'
    
    return jsonify(response_object)


@app.route('/api/upload_genemark', methods=['GET', 'POST'])
def upload_gm_file():
    """
    User uploads GeneMark file (ldata). Add data from file to the database.
    POST method uploads file.
    GET method adds data in file to database.
    """
    create_directory(UPLOAD_FOLDER)
    response_object = {'status': 'success'}

    if request.method == 'POST':
        file = request.files['file']
        if file and allowed_file(file.filename):
            file_name = secure_filename(file.filename)
            # generate GUUID
            file_id = str(uuid.uuid4())
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], file_name))
            # save file info into DB
            file = Files(id=file_id, name=file_name, \
                         date=datetime.datetime.now(datetime.timezone.utc))
            db.session.add(file)
            db.session.commit()
        response_object['message'] = f"{file_name} uploaded successfully."

    if request.method == 'GET':
        genemark_ldata_file = get_file("GeneMark_ldata")
        annotate.parse_genemark_ldata(genemark_ldata_file)
        annotate.compare()
        response_object['message'] = "Data added to database."

    return jsonify(response_object)


@app.route('/annotate_data', methods=['GET', 'POST'])
def annotate_data():
    """
    Compares DNA Master's predictions against GeneMark's. 
    GET method shows all the DNA Master predictions with a status and action item for each.
    """
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
        # f = open("example.gb", "r")
        # return f.read()
        # ----------------------------
        gb_file = get_file("GenBank")
        out_file = annotate.modify_gb(gb_file)
        f = open(out_file, "r")
        return f.read()
        


    return jsonify(response_object)


@app.route('/annotate_data/failed/<cds_id>', methods=['GET', 'POST', 'PUT'])
def failed_annotation(cds_id):
    """
    Gives user new start options for CDSs with 'failed' status.
    GET method parses through failed genes and finds new start options.
    PUT method updates the start position if the user chooses to do so.
    """
    response_object = {'status': 'success'}

    if request.method == "GET":
        cds = DNAMaster.query.filter_by(id=cds_id).first()
        response_object['cds'] = {'id': cds.id,
                                'start': cds.start,
                                'stop': cds.stop,
                                'strand': cds.strand,
                                'function': cds.function,
                                'status': cds.status}
        
        fasta_file = get_file("Fasta")
        genemark_gdata_file = get_file("GeneMark_gdata")
        starts = annotate.failed_gene(cds_id, fasta_file, genemark_gdata_file)

        start_options = []
        # modal_options = []
        lowest_start = 2**100
        for start, frames in starts.items():
            frames['start_position'] = start
            if start < lowest_start:
                lowest_start = start
            start_options.append(frames)
        response_object['lowest_start'] = lowest_start
        response_object['start_options'] = start_options

        gdata_df = pd.read_csv(genemark_gdata_file, sep='\t', skiprows=16)
        gdata_df.columns = ['Base', '1', '2', '3', '4', '5', '6']
        gdata_df = gdata_df[gdata_df.Base.isin(range(cds.start-100, cds.stop+100))]
        # response_object['base'] = list(gdata_df['Base'])
        # response_object['one'] = list(gdata_df['1'])
        direct_graph = graph.make_graph_direct(gdata_df, cds.start, cds.stop, start_options)
        comp_graph = graph.make_graph_complementary(gdata_df, cds.start, cds.stop, start_options)
        response_object['direct_graph'] = direct_graph
        response_object['comp_graph'] = comp_graph

    if request.method == "POST":
        print("In fail post")
        cds = DNAMaster.query.filter_by(id=cds_id).first()
        fasta_file = get_file("Fasta")
        E_VALUE_THRESH = 1e-35
        print("starting blast")
        blast_results = blast.run_blast(fasta_file, cds.start, cds.stop, E_VALUE_THRESH)
        print("writing results to post response")
        response_object['blast'] = blast_results

    if request.method == "PUT":
        put_data = request.get_json()
        cds = DNAMaster.query.filter_by(id=put_data.get('id')).first()
        if cds:
            cds.start = put_data.get('start')
            cds.function = put_data.get('function')
            cds.status = put_data.get('status')
            db.session.commit()
            response_object['message'] = 'CDS updated!'
        else:
            response_object['message'] = 'CDS did not update.'

    return jsonify(response_object)


@app.route('/annotate_data/more/<cds_id>', methods=['GET', 'PUT', 'POST'])
def do_more_annotation(cds_id):
    """

    """
    response_object = {'status': 'success'}

    if request.method == "GET":
        cds = DNAMaster.query.filter_by(id=cds_id).first()
        response_object['cds'] = {'id': cds.id,
                                'start': cds.start,
                                'stop': cds.stop,
                                'strand': cds.strand,
                                'function': cds.function,
                                'status': cds.status}
        genemark_gdata_file = get_file("GeneMark_gdata")
        probs = annotate.need_more_info_genes(cds_id, genemark_gdata_file)
        response_object['probs'] = probs
        #---
        fasta_file = get_file("Fasta")

        gdata_df = pd.read_csv(genemark_gdata_file, sep='\t', skiprows=16)
        gdata_df.columns = ['Base', '1', '2', '3', '4', '5', '6']
        gdata_df = gdata_df[gdata_df.Base.isin(range(cds.start-100, cds.stop+100))]
        direct_graph = graph.make_graph_direct(gdata_df, cds.start, cds.stop, [])
        comp_graph = graph.make_graph_complementary(gdata_df, cds.start, cds.stop, [])
        response_object['direct_graph'] = direct_graph
        response_object['comp_graph'] = comp_graph

    if request.method == "POST":
        print("In more post")
        cds = DNAMaster.query.filter_by(id=cds_id).first()
        fasta_file = get_file("Fasta")
        E_VALUE_THRESH = 1e-35
        print("starting blast")
        blast_results = blast.run_blast(fasta_file, cds.start, cds.stop, E_VALUE_THRESH)
        print("writing results to post response")
        response_object['blast'] = blast_results
        print(blast_results)

    if request.method == "PUT":
        put_data = request.get_json()
        cds = DNAMaster.query.filter_by(id=put_data.get('id')).first()
        if cds:
            # cds.start = put_data.get('start')
            cds.function = put_data.get('function')
            cds.status = put_data.get('status')
            db.session.commit()
            response_object['message'] = 'CDS updated!'
        else:
            response_object['message'] = 'CDS did not update.'

    return jsonify(response_object)


@app.route('/blast/<cds_id>', methods=['GET', 'POST', 'PUT'])
def blast_for_function(cds_id):
    """

    """
    response_object = {'status': 'success'}

    if request.method == "GET":
        cds = DNAMaster.query.filter_by(id=cds_id).first()
        response_object['cds'] = {'id': cds.id,
                                'start': cds.start,
                                'stop': cds.stop,
                                'strand': cds.strand,
                                'function': cds.function,
                                'status': cds.status}
        
    if request.method == "POST":
        print("In pass post")
        cds = DNAMaster.query.filter_by(id=cds_id).first()
        fasta_file = get_file("Fasta")
        E_VALUE_THRESH = 1e-35
        print("starting blast")
        blast_results = blast.run_blast(fasta_file, cds.start, cds.stop, E_VALUE_THRESH)
        print("writing results to post response")
        response_object['blast'] = blast_results

    if request.method == "PUT":
        put_data = request.get_json()
        cds = DNAMaster.query.filter_by(id=put_data.get('id')).first()
        if cds:
            cds.function = put_data.get('function')
            db.session.commit()
            response_object['message'] = 'CDS updated!'
        else:
            response_object['message'] = 'CDS did not update.'
    
    return jsonify(response_object)


if __name__ == '__main__':
    app.run(debug=True)


# ---------- HELPER FUNCTIONS ----------

# check file extension. only allow specific ones
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

# create directory
def create_directory(directory):
    try:
        os.mkdir(directory)
        print("Directory \'" + directory + "\' created.")
    except FileExistsError:
        print("Directory \'" + directory + "\' already exists.")

# Get necessary file from uploads dir depending on preference
def get_file(preference):
    for filename in os.listdir("uploads"):
        if preference == "GenBank":
            if filename.endswith(".gb") or filename.endswith(".genbank"):
                return os.path.join("uploads", filename)
        elif preference == "GeneMark_ldata":
            if filename.endswith(".ldata"):
                return os.path.join("uploads", filename)
        elif preference == "GeneMark_gdata":
            if filename.endswith(".gdata"):
                return os.path.join("uploads", filename)
        elif preference == "Fasta":
            if filename.endswith(".fasta"):
                return os.path.join("uploads", filename)
        else:
            return("Couldn't find file.")

