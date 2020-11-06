"""Contains the functions for the DNA Master page.

Gets all of the DNAMaster data.
Updates a CDS.
Deletes a CDS.
"""
import os
import models
from models import *
import helper
from helper import *

# ------------------------------ MAIN FUNCTIONS ------------------------------
def get_all_cds():
    """Queries all cds data from the dnamaster database and returns it.

    Returns:
        A dictionary containing all of DNAMaster's data.
    """
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
    """Updates a cds given the data.
    
    Args:
        request:
            A dictionary containing the new data.
        UPLOAD_FOLDER:
            The directory containing all of the uploaded files.

    Returns:
        A dictionary containing a success message.
    """
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
    """Deletes a cds given the ID.

    Args:
        cds_id:
            The ID of the CDS to be deleted.
        UPLOAD_FOLDER:
            The directory containing all of the uploaded files.

    Returns:
        A dictionary containing a success message.
    """
    response_object = {}
    if DNAMaster.query.filter_by(id=cds_id).first():
        DNAMaster.query.filter_by(id=cds_id).delete()
        db.session.commit()
        response_object['message'] = 'CDS removed!'
    else:
        response_object['message'] = 'Error: CDS could not be deleted.'

    response_object['status'] = delete_blast_zip(UPLOAD_FOLDER)

    return response_object
