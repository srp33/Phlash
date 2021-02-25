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
import helper
import shutil
import pandas as pd

def update_settings(payload):
    """Updates settings given new setting data.
    Args:
        payload:
            The new settings.
    Returns:
        A dictionary containing a success message.
    """
    response_object = {}
    print(payload)
    new_settings_data = payload.split(',')
    setting = db.session.query(Settings).order_by(Settings.back_start_range).first()
    setting.back_start_range = new_settings_data[3]
    setting.forward_start_range = new_settings_data[4]
    setting.gap = new_settings_data[0]
    setting.overlap = new_settings_data[1]
    setting.opposite_gap = new_settings_data[2]
    setting.short = new_settings_data[5]
    db.session.commit()
    response_object['message'] = 'Settings updated!'
    return response_object

def get_settings():
    """Returns the current settings.

        Returns:
            A dictionary containing the current settings.
    """
    response_object = {}
    setting = db.session.query(Settings).order_by(Settings.back_start_range).first()
    response_object['back_start_range'] = setting.back_start_range
    response_object['forward_start_range'] = setting.forward_start_range
    response_object['gap'] = setting.gap
    response_object['opposite_gap'] = setting.opposite_gap
    response_object['overlap'] = setting.overlap
    response_object['short'] = setting.short
    return response_object