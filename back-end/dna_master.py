"""
Contains the methods for the DNA Master page.
"""
import os
import models
from models import *
import helper
from helper import *

# ------------------------------ MAIN FUNCTIONS ------------------------------
def get_all_cds():
    '''
    Queries all cds data from the dnamaster database and returns it.
    '''
    response_object = {}
    dnamaster = []
    for cds in db.session.query(DNAMaster).order_by(DNAMaster.start):
        dnamaster.append({'id': cds.id,
                            'start': cds.start,
                            'stop': cds.stop,
                            'strand': cds.strand})
    response_object['dnamaster'] = dnamaster

    return response_object

def update_cds(request, UPLOAD_FOLDER):
    '''
    Updates a cds given the data.
    '''
    response_object = {}
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

    response_object['status'] = delete_blast_zip(UPLOAD_FOLDER)

    return response_object

def delete_cds(cds_id, UPLOAD_FOLDER):
    '''
    Deletes a cds given the ID.
    '''
    response_object = {}
    if DNAMaster.query.filter_by(id=cds_id).first():
        DNAMaster.query.filter_by(id=cds_id).delete()
        db.session.commit()
        response_object['message'] = 'CDS removed!'
    else:
        response_object['message'] = 'Error: CDS could not be deleted.'

    response_object['status'] = delete_blast_zip(UPLOAD_FOLDER)

    return response_object

# ---------- DELETE AND UPDATE HELPER FUNCTIONS ----------
