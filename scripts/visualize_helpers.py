import glob
import os
import h5py
import numpy as np
import yt

from units import M_solar, m_proton, pc, yr, Myr, gamma


##########################
snapshot_filename_format = "snapshot_???.hdf5"

def get_snapshot_filenames(outputs_dir, 
                           snapshot_filename_format=snapshot_filename_format):

    return sorted(glob.glob(os.path.join(outputs_dir, snapshot_filename_format)))


def load_snapshots(outputs_dir):
    unit_base = {
        "UnitLength_in_cm" : (pc),
        "UnitVelocity_in_cm_per_s"   : (pc / Myr),
        "UnitMass_in_g"   : (M_solar)
    }

    snapshot_filenames = get_snapshot_filenames(outputs_dir)

    n_files_ready = len(snapshot_filenames)
    if n_files_ready == 0:
        raise FileNotFoundError("No snapshots found in {}".format(outputs_dir))

    ts = yt.load(os.path.join(outputs_dir, snapshot_filename_format),
                 unit_base=unit_base)
    return ts


def get_snapshot_times(ts):
    """`ts` should be a yt timeseries"""
    snapshot_times = np.array([ds.current_time.convert_to_cgs() 
                               for ds in ts]) / Myr
    return snapshot_times



def total_mass_of_snapshot(snapshot_filename):
    f = h5py.File(snapshot_filename, mode="r")

    total_mass = np.sum(f["PartType0"]["Masses"], dtype=float)
    
    f.close()
    
    return total_mass


def total_radial_momentum_of_snapshot(snapshot_filename):
    f = h5py.File(snapshot_filename, mode="r")
    masses_shape = f["PartType0"]["Masses"].shape
    new_shape = masses_shape + (1,)

    r_hat = f["PartType0"]["Coordinates"] - (f["Header"].attrs["BoxSize"]/2)
    r_hat = r_hat / (np.sum(r_hat**2, axis=1, dtype=float).reshape(new_shape))**.5
    
    mom = np.sum(r_hat * f["PartType0"]["Velocities"] \
    * np.reshape(f["PartType0"]["Masses"], new_shape), dtype=float)
    
    f.close()
    return mom * M_solar * pc / (Myr)


def total_kinetic_energy_of_snapshot(snapshot_filename):
    f = h5py.File(snapshot_filename, mode="r")
    masses_shape = f["PartType0"]["Masses"].shape
    new_shape = masses_shape + (1,)

    E_kin = 0.5 * np.sum(np.array(f["PartType0"]["Velocities"], dtype=float)**2 \
    * np.reshape(f["PartType0"]["Masses"], new_shape), dtype=float)
    
    f.close()
    return E_kin * M_solar * (pc / Myr)**2


def total_internal_energy_of_snapshot(snapshot_filename):
    f = h5py.File(snapshot_filename, mode="r")
    masses_shape = f["PartType0"]["Masses"].shape
    new_shape = masses_shape + (1,)

    E_int = np.sum(f["PartType0"]["InternalEnergy"] \
    * f["PartType0"]["Masses"], dtype=float)
    
    f.close()
    return E_int * M_solar * (pc / Myr)**2

def total_internal_energy_of_snapshot(snapshot_filename):
    E_tot =   total_internal_energy_of_snapshot(snapshot_filename) \
            + total_kinetic_energy_of_snapshot(snapshot_filename)
    return E_tot



def map_to_all_snapshots(outputs_dir, mapped_function):
    """To do: option for multiprocessing map?"""
    snapshot_filenames = get_snapshot_filenames(outputs_dir)

    results = np.array(list(map(mapped_function, snapshot_filenames)), ndmin=1)

    return results


