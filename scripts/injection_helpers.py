import h5py
import shutil
import os
import glob
import numpy as np
from scipy import stats

from units import M_solar, m_proton, pc, yr, Myr, km, s, gamma

default_SNe_datafile = "1D_data/25451948-485f-46fe-b87b-f4329d03b203_SNe.dat"


##### CLASSES ###########


class Supernova(object):
    def __init__(self, time, ejecta_mass, ejecta_mass_Z):
        self.time          = time
        self.ejecta_mass   = ejecta_mass
        self.ejecta_mass_Z = ejecta_mass_Z



##### FUNCTIONS ###########
SNe_filename_format = "*SNe.dat"
def get_SNe(SNe_dir, SNe_filename_format=SNe_filename_format):
    """Load datafile from SNe_dir, parse into a list of Supernova objects"""
    possible_SN_files = glob.glob(os.path.join(SNe_dir, SNe_filename_format))

    if len(possible_SN_files) == 0: 
        raise FileNotFoundError("No SN data files found in {}".format(SNe_dir))
    elif len(possible_SN_files) > 1:
        raise RuntimeError("Too many SN data files found in {}".format(SNe_dir))
    SN_file = possible_SN_files[0]


    SN_data = np.loadtxt(SN_file, ndmin=2)
    SN_data = SN_data[np.argsort(SN_data[:,0])]

    # change into the code units I'm using in GIZMO
    # Base units: Gyr, M_solar, pc
    SN_times         = SN_data[:,0] / (Myr)
    SN_ejecta_mass   = SN_data[:,2] / M_solar
    SN_ejecta_mass_Z = SN_data[:,3] / M_solar

    SN_times -= SN_times[0]

    SNe = [Supernova(SN_times[i], SN_ejecta_mass[i], SN_ejecta_mass_Z[i])
           for i in range(SN_times.size)]

    return SNe


def snapshot_number_from_basename(basename):
    return int(basename.replace("snapshot_", "").replace(".hdf5", ""))

def snapshot_basename_from_number(snapshot_number):
    return "snapshot_{:03d}.hdf5".format(snapshot_number)


def get_last_snapshot_file_in_dir(outputs_dir):
    snapshot_files = glob.glob(os.path.join(outputs_dir, "snapshot_*.hdf5"))
    snapshot_numbers = [snapshot_number_from_basename(os.path.basename(f)) 
                        for f in snapshot_files]
    
    i_max = np.argmax(snapshot_numbers)

    last_snapshot_file   = snapshot_files[i_max]
    
    return last_snapshot_file
    

def which_SN_is_about_to_explode(current_time, SNe, tol=1e-4):
    """Returns the *zero-indexed* number of the SNe about to explode"""

    SN_times = np.array([SN.time for SN in SNe])
        
    possible_SNe_about_to_occur = np.where(   (SN_times <= (current_time + 1e-3 + tol)) \
                                            & (SN_times > current_time-tol))[0]
    
    if possible_SNe_about_to_occur.size == 0:
        raise RuntimeError("No SNe are going to explode within about 1 kyr")
    elif possible_SNe_about_to_occur.size > 1:
        raise RuntimeError("More than 1 SN are going to explode within 1 kyr")
    
    return possible_SNe_about_to_occur[0]


def find_params_file_base(inputs_dir):
    
    possible_base_files = glob.glob(os.path.join(inputs_dir, "*params.base"))
    if len(possible_base_files) == 0:
        raise RuntimeError("No base params files found")
    elif len(possible_base_files) > 1:
        raise RuntimeError("More than one base params file found")
        
    params_file_base = possible_base_files[0]
    
    return params_file_base


def get_default_snapshot_spacing(params_file):
    # try one iteration ignoring commented lines
    for line in open(params_file, mode="r"):
        if line[0] == "%":
            continue
        if "TimeBetSnapshot" in line:
            return float(line.split()[-1])
    
    # now only look at commented lines if nothing was found in commented lines
    for line in open(params_file, mode="r"):
        if "TimeBetSnapshot" in line:
            return float(line.split()[-1])
        
    raise RuntimeError("No lines found with `TimeBetSnapshot`")


def create_restart_params(snapshot_file_after_SN, inputs_dir, SNe):
    """There must be a SNe between `snapshot_file_after_SN` and the previous snapshot"""


    snapshot_number_after_SN = snapshot_number_from_basename(
                                    os.path.basename(snapshot_file_after_SN)
                                )

    snapshot_file_before_SN = snapshot_file_after_SN.replace(
                                    "{:03}".format(snapshot_number_after_SN),
                                    "{:03}".format(snapshot_number_after_SN-1),
                                )

    with h5py.File(snapshot_file_before_SN) as f:
        time_before = f["Header"].attrs["Time"]

    with h5py.File(snapshot_file_after_SN) as f:
        time_after = f["Header"].attrs["Time"]    


    SN_times = np.array([SN.time for SN in SNe])
    possible_SNe = np.where((SN_times >= time_before) & (SN_times < time_after))[0]
    if possible_SNe.size == 0:
        raise RuntimeError("No SNe exploded between 2 most recent snapshots")
    elif possible_SNe.size > 1:
        raise RuntimeError("More than 1 SN exploded between 2 most recent snapshots")
    i_SN = possible_SNe[0]

    params_file_base = find_params_file_base(inputs_dir)
    params_new_file = params_file_base.replace("base", "restart")
    shutil.copy2(params_file_base, params_new_file)
    
    with open(params_new_file, mode="a") as f:
        
        if i_SN+1 < len(SNe):
            time_of_previous_SN = SNe[i_SN].time
            time_of_next_SN     = SNe[i_SN+1].time
        else:
            time_of_previous_SN = SNe[i_SN].time
            time_of_next_SN     = 40

        
        time_min = time_of_previous_SN + 1e-3
        time_max = time_of_next_SN     - 1e-3

        
        default_snapshot_spacing = get_default_snapshot_spacing(params_file_base)
        
        num_snapshots = max(5, int((time_max - time_min)/default_snapshot_spacing))

        snapshot_spacing = (time_max - time_min) / num_snapshots
        
        print("InitCondFile                       {}".format(
            os.path.join(*snapshot_file_after_SN.replace(".hdf5", "").split(os.sep)[1:])), file=f
        )

        
        print("TimeBegin                          {}".format(time_min), file=f)
        print("TimeMax                            {}".format(time_max), file=f)

        
        print("TimeBetSnapshot                    {}".format(snapshot_spacing), file=f)
        print("TimeOfFirstSnapshot                {}".format(time_min + snapshot_spacing), file=f)
        
