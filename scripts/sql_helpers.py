from __future__ import print_function, division

# %matplotlib inline
# from matplotlib import pyplot as plt
# import seaborn as sns; sns.set(context="poster")
# import ipywidgets
# import yt
import glob
import os
import shutil
# import warnings
# import h5py

import numpy as np
import pandas as pd

from units import M_solar, m_proton, pc, yr, Myr, km, s, gamma

# from injection_helpers import get_SNe

from visualize_helpers import \
    get_snapshot_filenames, \
    snapshot_filename_to_number, \
    get_snapshot_time

from snapshot_helpers import \
    total_mass_of_snapshot, \
    total_radial_momentum_of_snapshot, \
    total_kinetic_energy_of_snapshot, \
    total_internal_energy_of_snapshot

import sqlite3

########################## CLASSES


class Snapshot(object):
    """basic struct representing a record in the a `snapshots` table
    
        Inputs / Attributes
        ------
        number   : integer 
        time     : double   - Myr
        e_kin    : double   - ergs
        e_int    : double   - ergs
        e_mag    : double   - ergs
        mass     : double   - M_solar
        momentum : double   - g cm s^-1
        run_name : str
        """
    def __init__(self, number, time, e_kin, e_int, e_mag, mass, momentum, run_name):
        super(Snapshot, self).__init__()
        
        self.number   = number
        self.time     = time
        self.e_kin    = e_kin
        self.e_int    = e_int
        self.e_mag    = e_mag
        self.mass     = mass
        self.momentum = momentum
        self.run_name = run_name


########################## FUNCTIONS
def get_db_dirname_tmp():
    tmp_drive = "/dev/shm"
    user = "egentry"

    db_dirname = os.path.join(tmp_drive, user)
    return db_dirname

def get_db_filename_tmp(run_name):
    db_filename = os.path.join(get_db_dirname_tmp(), "{}.db".format(run_name))
    return db_filename

def get_db_filename_permanent(run_name):
    runs_dir = "../runs"
    db_filename = os.path.join(runs_dir, "{0}/outputs/{0}.db".format(run_name))
    return db_filename

def copy_database_from_permanent_to_tmp(run_name):
    db_filename_tmp       = get_db_filename_tmp(run_name)
    db_filename_permanent = get_db_filename_permanent(run_name)
    
    shutil.copy2(db_filename_permanent, db_filename_tmp)

def copy_database_from_tmp_to_permanent(run_name):
    db_filename_tmp       = get_db_filename_tmp(run_name)
    db_filename_permanent = get_db_filename_permanent(run_name)
    
    shutil.copy2(db_filename_tmp, db_filename_permanent)

def open_connection(run_name):
    db_filename = get_db_filename_tmp(run_name)
    conn = sqlite3.connect(db_filename)
    return conn

def create_new_table(run_name):
    """run_name should be a string (e.g. `cluster_cooling_100`)"""
    conn = open_connection(run_name)
    c = conn.cursor()
    

    try:
        c.execute("""CREATE TABLE snapshots 
                (number INTEGER PRIMARY KEY, 
                time DOUBLE PRECISION,
                e_kin DOUBLE PRECISION,
                e_int DOUBLE PRECISION,
                e_mag DOUBLE PRECISION,
                mass DOUBLE PRECISION,
                momentum DOUBLE PRECISION,
                run_name date text) """)
    except:
        conn.close()
        raise
    
    conn.commit()
    conn.close()


def add_entry(snapshot):
    conn = open_connection(snapshot.run_name)
    c = conn.cursor()

    add_command = "INSERT INTO snapshots VALUES (?,?,?,?,?,?,?,?)"
    
    entry = [(snapshot.number,
              snapshot.time,
              snapshot.e_kin,
              snapshot.e_int,
              snapshot.e_mag,
              snapshot.mass,
              snapshot.momentum,
              snapshot.run_name,
             ),
             ]
    
    try:
        c.executemany(add_command, entry)
    except:
        conn.close()
        raise
    
    conn.commit()
    conn.close()
    
def is_snapshot_missing(run_name, snapsnot_number):
    """snapshot_number should be an int"""
    conn = open_connection(run_name)
    c = conn.cursor()
    get_existing_primary_keys_command = "SELECT number from snapshots"
    
    existing_primary_keys = c.execute(get_existing_primary_keys_command).fetchall()

    conn.close()

    if (snapsnot_number,) in existing_primary_keys:
        return False
    else:
        return True

def add_entry_if_missing(snapshot):
    if is_snapshot_missing(snapshot.run_name, snapshot.number):
        add_entry(snapshot)


def open_as_DataFrame(run_name, tmp_location=True, copy_first=False):
    if tmp_location:
        db_filename = get_db_filename_tmp(run_name)
        if copy_first:
            copy_database_from_permanent_to_tmp(run_name)
    else:
        print("Warning: sqlite might not work on a pfs drive")
        db_filename = get_db_filename_permanent(run_name)
    df = pd.read_sql_table("snapshots", "sqlite:///{}".format(db_filename))
    return df

def add_snapshot_if_missing(snapshot_filename, run_name, verbose=False):
    snapshot_number = snapshot_filename_to_number(snapshot_filename)
    if not is_snapshot_missing(run_name, snapshot_number):
        if verbose:
            print("Snapshot #{} already found for run: {}".format(snapshot_number, run_name) )
        return

    mass = total_mass_of_snapshot(snapshot_filename)
    e_int = total_internal_energy_of_snapshot(snapshot_filename)
    e_kin = total_kinetic_energy_of_snapshot(snapshot_filename)
    e_mag = 0
    momentum = total_radial_momentum_of_snapshot(snapshot_filename)
    time = get_snapshot_time(snapshot_filename)
    
    snapshot = Snapshot(snapshot_number, time, e_kin, e_int, e_mag, mass, momentum, run_name)
    
    add_entry_if_missing(snapshot)

def add_simulation(run_name, verbose=True):
    try:
        tmp_dir = get_db_dirname_tmp()
        os.makedirs(tmp_dir, exist_ok=True)
        copy_database_from_permanent_to_tmp(run_name)
    except FileNotFoundError:
        if not os.path.exists(get_db_filename_tmp(run_name)):
            create_new_table(run_name)
        else:
            raise


    run_dir = "../runs/{}/".format(run_name)
    outputs_dir = os.path.join(run_dir, "outputs")

    snapshot_filenames = get_snapshot_filenames(outputs_dir)

    for snapshot_filename in snapshot_filenames:
        if verbose:
            print("adding: ", snapshot_filename)
        add_snapshot_if_missing(snapshot_filename, run_name, verbose=verbose)

    copy_database_from_tmp_to_permanent(run_name)
    
    
