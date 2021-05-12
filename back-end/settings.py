"""Contains the functions for the settings modal
    Gets the current settings.
    Updates settings.
"""

from werkzeug.utils import secure_filename
from Bio import SeqIO, Seq
from models import *
import os
from builtins import FileNotFoundError
import subprocess
import models
import shutil
import pandas as pd

def update_settings(phage_id, payload):
    """Updates settings given new setting data.
    Args:
        payload:
            The new settings.
    Returns:
        A dictionary containing a success message.
    """
    response_object = {}
    new_settings_data = payload.split(',')
    setting = db.session.query(Settings).filter_by(phage_id=phage_id).order_by(Settings.back_left_range).first()
    setting.back_left_range = new_settings_data[3]
    setting.forward_left_range = new_settings_data[4]
    setting.gap = new_settings_data[0]
    setting.overlap = new_settings_data[1]
    setting.opposite_gap = new_settings_data[2]
    setting.short = new_settings_data[5]
    db.session.commit()
    response_object['message'] = 'Settings updated!'
    return response_object

def get_settings(phage_id):
    """Returns the current settings.

        Returns:
            A dictionary containing the current settings.
    """
    response_object = {}
    setting = db.session.query(Settings).filter_by(phage_id=phage_id).order_by(Settings.back_left_range).first()
    response_object['back_left_range'] = setting.back_left_range
    response_object['forward_left_range'] = setting.forward_left_range
    response_object['gap'] = setting.gap
    response_object['opposite_gap'] = setting.opposite_gap
    response_object['overlap'] = setting.overlap
    response_object['short'] = setting.short
    return response_object