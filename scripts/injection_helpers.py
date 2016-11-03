import h5py
import shutil
import os
import glob
import numpy as np
from scipy import stats
from enum import Enum, unique

from collections import OrderedDict


from units import M_solar, m_proton, pc, yr, Myr, km, s, gamma
from astropy import units as u

default_SNe_datafile = "1D_data/25451948-485f-46fe-b87b-f4329d03b203_SNe.dat"




##### CLASSES ###########


class Supernova(object):
    def __init__(self, time, ejecta_mass, ejecta_mass_Z):
        self.time          = time
        self.ejecta_mass   = ejecta_mass
        self.ejecta_mass_Z = ejecta_mass_Z

@unique
class RestartType(Enum):
    RESTARTFILE = 1
    SNAPSHOTS = 2


class Params(object):
    types = OrderedDict([
    ("OutputDir",                          str),
    ("RestartFile",                        str),
    ("SnapshotFileBase",                   str),
    ("OutputListFilename",                 str),
    ("ICFormat",                           int),
    ("SnapFormat",                         int),
    ("TimeLimitCPU",                       int),
    ("CpuTimeBetRestartFile",              int),
    ("ResubmitOn",                         int),
    ("ResubmitCommand",                    str),
    ("MaxMemSize",                         int),
    ("PartAllocFactor",                    float),
    ("BufferSize",                         float),
    ("ComovingIntegrationOn",              float),
    ("Omega0",                             float),
    ("OmegaLambda",                        float),
    ("OmegaBaryon",                        float),
    ("HubbleParam",                        float),
    ("BoxSize",                            float),
    ("OutputListOn",                       int),
    ("TimeBetStatistics",                  float),
    ("NumFilesPerSnapshot",                int),
    ("NumFilesWrittenInParallel",          int),
    ("ErrTolIntAccuracy",                  float),
    ("CourantFac",                         float),
    ("MaxRMSDisplacementFac",              float),
    ("MaxSizeTimestep",                    float),
    ("MinSizeTimestep",                    float),
    ("ErrTolTheta",                        float),
    ("ErrTolForceAcc",                     float),
    ("TreeDomainUpdateFrequency",          float),
    ("DesNumNgb",                          int),
    ("MaxNumNgbDeviation",                 float),
    ("ArtBulkViscConst",                   float),
    ("UnitLength_in_cm",                   float),
    ("UnitMass_in_g",                      float),
    ("UnitVelocity_in_cm_per_s",           float),
    ("GravityConstantInternal",            float),
    ("MaxHsml",                            float),
    ("MinGasHsmlFractional",               float),
    ("SofteningGas",                       float),
    ("SofteningHalo",                      float),
    ("SofteningDisk",                      float),
    ("SofteningBulge",                     float),
    ("SofteningStars",                     float),
    ("SofteningBndry",                     float),
    ("SofteningGasMaxPhys",                float),
    ("SofteningHaloMaxPhys",               float),
    ("SofteningDiskMaxPhys",               float),
    ("SofteningBulgeMaxPhys",              float),
    ("SofteningStarsMaxPhys",              float),
    ("SofteningBndryMaxPhys",              float),
    ("InitGasTemp",                        float),
    ("MinGasTemp",                         float),
    ("InitMetallicity",                    int),
    ("MetalCoolingOn",                     int),
    ("UVBackgroundOn",                     int),
    ("GrackleDataFile",                    str),
    ("SNeDataFile",                        str),
    ("InitCondFile",                       str),
    ("TimeBegin",                          float),
    ("TimeMax",                            float),
    ("TimeBetSnapshot",                    float),
    ("TimeOfFirstSnapshot",                float),
    ])
    
    def __init__(self, **kwargs):
        for key in kwargs:
            if key not in self.types.keys():
                raise RuntimeError("'{}' identifier not recognized".format(key))
        for key in self.types.keys():
            if key not in kwargs:
                raise RuntimeError("missing entry for: {}".format(key))
        self.__dict__ = kwargs
        

    @classmethod
    def from_filename(cls, filename):
        param_dict = {}
        with open(filename, mode="r") as f:
            for line in f:
                line = line.rstrip()
                if len(line) == 0:
                    continue
                if line[0] is "%":
                    continue
                key, value = line.split()
                param_dict[key] = cls.types[key](value)
        p = Params(**param_dict)
        return p
                        

    def to_file(self, file):
        """file must be a file stream, not a filename!"""
        for key in self.__class__.types:
            value = self.__dict__[key]
            print("{0:35s}{1:<45}".format(key, value), file=file)
        
    def copy(self):
        return Params(self.__dict__)



##### FUNCTIONS ###########
SNe_filename_format = "*SNe.dat"
def get_SNe(inputs_dir, SNe_filename_format=SNe_filename_format):
    """Load datafile from inputs_dir, parse into a list of Supernova objects"""
    possible_SN_files = glob.glob(os.path.join(inputs_dir, SNe_filename_format))

    if len(possible_SN_files) == 0: 
        raise FileNotFoundError("No SN data files found in {}".format(inputs_dir))
    elif len(possible_SN_files) > 1:
        raise RuntimeError("Too many SN data files found in {}".format(inputs_dir))
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
    return get_ith_snapshot_file_in_dir(outputs_dir, -1)

def get_ith_snapshot_file_in_dir(outputs_dir, i):
    snapshot_files = np.array(glob.glob(os.path.join(outputs_dir, "snapshot_*.hdf5")))
    snapshot_numbers = np.array([snapshot_number_from_basename(os.path.basename(f)) 
                        for f in snapshot_files])

    sorted_indices = np.argsort(snapshot_numbers)
    snapshot_files = snapshot_files[sorted_indices]
    
    return snapshot_files[i]
    

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
        

# clean this up (copied from a notebook)
def create_snapshot_with_new_SN(run_dir):
    inputs_dir  = os.path.join(run_dir, "inputs")
    outputs_dir = os.path.join(run_dir, "outputs")
    
    SNe = get_SNe(inputs_dir)
    SN_times           = np.array([SN.time          for SN in SNe])
    SN_ejecta_masses   = np.array([SN.ejecta_mass   for SN in SNe])
    SN_ejecta_masses_Z = np.array([SN.ejecta_mass_Z for SN in SNe])

    snapshot_file_old = get_last_snapshot_file_in_dir(outputs_dir)
    snapshot_number_old = snapshot_number_from_basename(os.path.basename(snapshot_file_old))

    snapshot_file_new = os.path.join(outputs_dir, snapshot_basename_from_number(snapshot_number_old+1))


    f_old = h5py.File(snapshot_file_old, mode="r")
    # this next line should throw a RuntimeError if you shouldn't add a new SN now
    i_SN = which_SN_is_about_to_explode(f_old["Header"].attrs["Time"], SNe)

    snapshot_file_new_tmp = os.path.join(outputs_dir, "tmp.hdf5")
    f_new = h5py.File(snapshot_file_new_tmp, mode="w")

    f_new.create_group("Header")
    f_new.create_group("PartType0")

    SN = SNe[i_SN]



    average_particle_mass = np.sum(f_old["PartType0"]["Masses"]) / f_old["Header"].attrs["NumPart_Total"][0]
    # print(average_particle_mass)

    num_new_particles_needed = int(np.ceil(SN.ejecta_mass / average_particle_mass))
    # print(num_new_particles_needed)

    mass_per_new_particle = SN.ejecta_mass / num_new_particles_needed
    # print(mass_per_new_particle)


    for key in f_old["Header"].attrs:
        f_new["Header"].attrs[key] = f_old["Header"].attrs[key]

        
    # This is necessary since h5py attrs don't allow slicing/element-wise changes
    new_NumPart_ThisFile = f_new["Header"].attrs["NumPart_ThisFile"]
    new_NumPart_Total    = f_new["Header"].attrs["NumPart_Total"]

    new_NumPart_ThisFile[0] += num_new_particles_needed
    new_NumPart_Total[0]    += num_new_particles_needed

    f_new["Header"].attrs["NumPart_ThisFile"] = new_NumPart_ThisFile
    f_new["Header"].attrs["NumPart_Total"]    = new_NumPart_Total

    old_total_particles = f_old["Header"].attrs["NumPart_Total"][0]
    new_total_particles = f_new["Header"].attrs["NumPart_Total"][0]


    # Make other changes
    f_new["Header"].attrs["Time"] = SN.time + 1e-3 # 1 kyr after the SN


    # check data
    for key in f_new["Header"].attrs:
        # print(key, ":", f_new["Header"].attrs[key])
        pass

    for key in f_old["Header"].attrs:
        # print(key, ":", f_old["Header"].attrs[key]) 
        pass

    for key in f_old["PartType0"].keys():    
        new_shape = (new_total_particles,) + f_old["PartType0"][key].shape[1:]
        dtype = f_old["PartType0"][key].dtype
        
        f_new["PartType0"].require_dataset(key, new_shape, dtype=dtype)
        
        f_new["PartType0"][key][:old_total_particles] = f_old["PartType0"][key]
        
        # print(key, ":", f_new["PartType0"][key])
        
    r_SN = 2 # pc    
    injection_kernel = stats.norm(loc=f_new["Header"].attrs["BoxSize"]/2, scale=r_SN)
    f_new["PartType0"]["Coordinates"][-num_new_particles_needed:,0] = injection_kernel.rvs(size=num_new_particles_needed)
    f_new["PartType0"]["Coordinates"][-num_new_particles_needed:,1] = injection_kernel.rvs(size=num_new_particles_needed)
    f_new["PartType0"]["Coordinates"][-num_new_particles_needed:,2] = injection_kernel.rvs(size=num_new_particles_needed)

    f_new["PartType0"]["Density"][-num_new_particles_needed:] = np.mean(f_old["PartType0"]["Density"]) # EDIT: will this mess with my yt plots of this snapshot?

    f_new["PartType0"]["ElectronAbundance"][-num_new_particles_needed:] = 1.0 # EDIT: unused?

    energy_in_correct_units = 1e51 * u.erg.to(u.pc**2 * u.M_sun / u.Myr**2)
    f_new["PartType0"]["InternalEnergy"][-num_new_particles_needed:] = energy_in_correct_units/SN.ejecta_mass

    f_new["PartType0"]["Masses"][-num_new_particles_needed:] = mass_per_new_particle

    f_new["PartType0"]["Metallicity"][-num_new_particles_needed:] = SN.ejecta_mass / SN.ejecta_mass_Z

    f_new["PartType0"]["NeutralHydrogenAbundance"][-num_new_particles_needed:] = 0.0 # EDIT: unused?

    f_new["PartType0"]["ParticleChildIDsNumber"][-num_new_particles_needed:] = 0 # is this right?

    f_new["PartType0"]["ParticleIDGenerationNumber"][-num_new_particles_needed:] = 0 # is this right?

    new_particle_ids = np.max(f_old["PartType0"]["ParticleIDs"]) + 1 + np.arange(num_new_particles_needed)
    f_new["PartType0"]["ParticleIDs"][-num_new_particles_needed:] = new_particle_ids

    f_new["PartType0"]["SmoothingLength"][-num_new_particles_needed:] = np.mean(f_old["PartType0"]["SmoothingLength"])

    f_new["PartType0"]["Velocities"][-num_new_particles_needed:,0] = 0
    f_new["PartType0"]["Velocities"][-num_new_particles_needed:,1] = 0
    f_new["PartType0"]["Velocities"][-num_new_particles_needed:,2] = 0


    f_old.close()
    f_new.close()
    os.rename(snapshot_file_new_tmp, snapshot_file_new)


def create_restart_type_file(inputs_dir, restart_type, basename="RESTART"):
    if restart_type not in {entry.value for entry in RestartType}:
        raise ValueError("restart_type must be a valid value of RestartType enum")

    with open(os.path.join(inputs_dir, basename), mode="w") as f:
        print("RESTART_TYPE={}".format(restart_type), file=f)

################### Wrapped functions (clean these up)



def create_restart_params_wrapped(run_dir):
    inputs_dir  = os.path.join(run_dir, "inputs")
    outputs_dir = os.path.join(run_dir, "outputs")

    SNe = get_SNe(inputs_dir)

    snapshot_file_after_SN = get_last_snapshot_file_in_dir(outputs_dir)
    snapshot_number_after_SN = snapshot_number_from_basename(os.path.basename(snapshot_file_after_SN))

    snapshot_file_before_SN = snapshot_file_after_SN.replace("{:03}".format(snapshot_number_after_SN),
                                                             "{:03}".format(snapshot_number_after_SN-1),)

    ## check that we should have added a SN between the most recent 2 snapshots
    ## throws a RuntimeError if that isn't true
    with h5py.File(snapshot_file_before_SN) as f:
        _ = which_SN_is_about_to_explode(f["Header"].attrs["Time"], SNe)


    create_restart_params(snapshot_file_after_SN, inputs_dir, SNe)


