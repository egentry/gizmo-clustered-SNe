import h5py
import shutil
import os
import glob
import numpy as np
from scipy import stats, interpolate
from enum import IntEnum, unique

from collections import OrderedDict

from units import M_solar, m_proton, pc, yr, Myr, km, s, gamma
from astropy import units as u

from visualize_helpers import total_mass_of_snapshot, \
                              total_internal_energy_of_snapshot, \
                              total_kinetic_energy_of_snapshot

default_SNe_datafile = "1D_data/25451948-485f-46fe-b87b-f4329d03b203_SNe.dat"


# #### CLASSES ###########


class Supernova(object):
    def __init__(self, time, ejecta_mass, ejecta_mass_Z):
        self.time = time
        self.ejecta_mass = ejecta_mass
        self.ejecta_mass_Z = ejecta_mass_Z


@unique
class RestartType(IntEnum):
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

    types_MHD = OrderedDict([
        ("UnitMagneticField_in_gauss",           float),
    ])

    def __init__(self, with_MHD=False, **kwargs):
        self.with_MHD = with_MHD
        for key in kwargs:
            if key not in (set(self.types.keys()) | set(self.types_MHD.keys())):
                raise RuntimeError("'{}' identifier not recognized".format(key))
        for key in self.types.keys():
            if key not in kwargs:
                raise RuntimeError("missing entry for: {}".format(key))
        if self.with_MHD:
            for key in self.types_MHD.keys():
                if key not in kwargs:
                    raise RuntimeError("missing entry for: {}".format(key))
        self.__dict__.update(kwargs)

    @classmethod
    def from_filename(cls, filename):
        if "mhd" in filename:
            with_MHD = True
        else:
            with_MHD = False

        param_dict = {}
        with open(filename, mode="r") as f:
            for line in f:
                line = line.rstrip()
                if len(line) == 0:
                    continue
                if line[0] is "%":
                    continue
                key, value = line.split()
                if key in cls.types:
                    param_dict[key] = cls.types[key](value)
                else:
                    param_dict[key] = cls.types_MHD[key](value)
        p = Params(with_MHD=with_MHD, **param_dict)
        return p

    def to_file(self, file):
        """file must be a file stream, not a filename!"""
        types = self.__class__.types.copy()
        if self.with_MHD:
            types.update(self.__class__.types_MHD)

        for key in types:
            value = self.__dict__[key]
            print("{0:35s}{1:<45}".format(key, value), file=file)

    def copy(self):
        return Params(with_MHD=self.with_MHD, **self.__dict__)


# #### FUNCTIONS ###########


SNe_filename_format = "*SNe.dat"


def get_SNe(inputs_dir, SNe_filename_format=SNe_filename_format):
    """Load datafile from inputs_dir, parse into a list of Supernova objects"""
    possible_SN_files = glob.glob(os.path.join(inputs_dir,
                                               SNe_filename_format))

    if len(possible_SN_files) == 0:
        raise FileNotFoundError("No SN data files found in {}".format(
            inputs_dir))
    elif len(possible_SN_files) > 1:
        raise RuntimeError("Too many SN data files found in {}".format(
            inputs_dir))
    SN_file = possible_SN_files[0]

    SN_data = np.loadtxt(SN_file, ndmin=2)
    SN_data = SN_data[np.argsort(SN_data[:, 0])]

    # change into the code units I'm using in GIZMO
    # Base units: Gyr, M_solar, pc
    SN_times         = SN_data[:,0] / Myr
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
    snapshot_files = np.array(glob.glob(os.path.join(outputs_dir,
                                                     "snapshot_*.hdf5")))
    snapshot_numbers = np.array([
        snapshot_number_from_basename(os.path.basename(f))
        for f in snapshot_files])

    sorted_indices = np.argsort(snapshot_numbers)
    snapshot_files = snapshot_files[sorted_indices]

    return snapshot_files[i]


def which_SN_is_about_to_explode(current_time, SNe, tol=1e-4):
    """Returns the *zero-indexed* number of the SNe about to explode"""

    SN_times = np.array([SN.time for SN in SNe])

    possible_SNe_about_to_occur = np.where(   (SN_times <= (current_time + 1e-3 + tol))
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


def create_restart_params(inputs_dir, outputs_dir, SNe):
    """There must be a SN between the two most recent snapshots"""

    snapshot_file_after_SN  = get_ith_snapshot_file_in_dir(outputs_dir, -1)
    snapshot_file_before_SN = get_ith_snapshot_file_in_dir(outputs_dir, -2)

    with h5py.File(snapshot_file_before_SN) as f:
        time_before = f["Header"].attrs["Time"]

    with h5py.File(snapshot_file_after_SN) as f:
        time_after = f["Header"].attrs["Time"]

    SN_times = np.array([SN.time for SN in SNe])
    possible_SNe = np.argwhere(  (SN_times >= time_before)
                               & (SN_times < time_after))[0]
    if possible_SNe.size == 0:
        raise RuntimeError("No SNe exploded between 2 most recent snapshots")
    elif possible_SNe.size > 1:
        raise RuntimeError("More than 1 SN exploded between 2 most recent snapshots")
    i_SN = possible_SNe[0]

    params_file_base = find_params_file_base(inputs_dir)
    params_new_file = params_file_base.replace("base", "restart")

    params = Params.from_filename(params_file_base)

    if i_SN+1 < len(SNe):
        time_of_previous_SN = SNe[i_SN].time
        time_of_next_SN     = SNe[i_SN+1].time
    else:
        time_of_previous_SN = SNe[i_SN].time
        time_of_next_SN     = 40

    time_min = time_of_previous_SN + 1e-3
    time_max = time_of_next_SN     - 1e-3

    # determine number of snapshots
    # (minimum of 5, maximum set by TimeBetSnapshot in `*.params.base`)
    num_snapshots = max(5, int((time_max - time_min)/params.TimeBetSnapshot))

    snapshot_spacing = (time_max - time_min) / num_snapshots

    params.TimeBetSnapshot = snapshot_spacing
    params.InitCondFile = os.path.join(*snapshot_file_after_SN.replace(".hdf5", "").split(os.sep)[1:])
    params.TimeBegin = time_min
    params.TimeMax   = time_max
    params.TimeOfFirstSnapshot = params.TimeBegin + params.TimeBetSnapshot

    with open(params_new_file, mode="w") as f:
        params.to_file(f)


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

    num_new_particles_needed = int(np.ceil(SN.ejecta_mass / average_particle_mass))

    mass_per_new_particle = SN.ejecta_mass / num_new_particles_needed

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

    print("old_total_particles: ", old_total_particles)

    # Make other changes
    f_new["Header"].attrs["Time"] = SN.time + 1e-3  # 1 kyr after the SN

    for key in f_old["PartType0"].keys():
        new_shape = (new_total_particles,) + f_old["PartType0"][key].shape[1:]
        dtype = f_old["PartType0"][key].dtype

        f_new["PartType0"].require_dataset(key, new_shape, dtype=dtype)

        f_new["PartType0"][key][:old_total_particles] = f_old["PartType0"][key]

    r_SN = 2  # pc
    injection_kernel = stats.norm(loc=f_new["Header"].attrs["BoxSize"]/2, scale=r_SN)
    xs = injection_kernel.rvs(size=num_new_particles_needed)
    ys = injection_kernel.rvs(size=num_new_particles_needed)
    zs = injection_kernel.rvs(size=num_new_particles_needed)

    f_new["PartType0"]["Coordinates"][-num_new_particles_needed:, 0] = xs
    f_new["PartType0"]["Coordinates"][-num_new_particles_needed:, 1] = ys
    f_new["PartType0"]["Coordinates"][-num_new_particles_needed:, 2] = zs

    old_density_interp = get_central_interpolator(f_old, "Density")
    old_density = old_density_interp([0, 0, 0])

    new_density = old_density * SN.ejecta_mass \
        * injection_kernel.pdf(xs) \
        * injection_kernel.pdf(ys) \
        * injection_kernel.pdf(zs)

    f_new["PartType0"]["Density"][-num_new_particles_needed:] = new_density

    f_new["PartType0"]["ElectronAbundance"][-num_new_particles_needed:] = 1.0  # EDIT: unused?

    energy_in_correct_units = 1e51 * u.erg.to(u.pc**2 * u.M_sun / u.Myr**2)
    f_new["PartType0"]["InternalEnergy"][-num_new_particles_needed:] = energy_in_correct_units/SN.ejecta_mass

    f_new["PartType0"]["Masses"][-num_new_particles_needed:] = mass_per_new_particle

    f_new["PartType0"]["Metallicity"][-num_new_particles_needed:] = SN.ejecta_mass / SN.ejecta_mass_Z

    f_new["PartType0"]["NeutralHydrogenAbundance"][-num_new_particles_needed:] = 0.0  # EDIT: unused?

    f_new["PartType0"]["ParticleChildIDsNumber"][-num_new_particles_needed:] = 0  # is this right?

    f_new["PartType0"]["ParticleIDGenerationNumber"][-num_new_particles_needed:] = 0  # is this right?

    new_particle_ids = np.max(f_old["PartType0"]["ParticleIDs"]) + 1 + np.arange(num_new_particles_needed)
    f_new["PartType0"]["ParticleIDs"][-num_new_particles_needed:] = new_particle_ids

    f_new["PartType0"]["SmoothingLength"][-num_new_particles_needed:] = np.mean(f_old["PartType0"]["SmoothingLength"])

    f_new["PartType0"]["Velocities"][-num_new_particles_needed:, 0] = 0
    f_new["PartType0"]["Velocities"][-num_new_particles_needed:, 1] = 0
    f_new["PartType0"]["Velocities"][-num_new_particles_needed:, 2] = 0

    if "MagneticField" in set(f_old["PartType0"]):
        with_MHD = True

        interp = get_central_interpolator(f_old, "MagneticField")

        # For now, interpolate all the the origin
        # In the future, maybe interpolate to their real locations?
        #   - I think that would still be divergence free, if you use a linear interpolator?
        B_new = interp([0, 0, 0])

        f_new["PartType0"]["MagneticField"][-num_new_particles_needed:] = B_new

        f_new["PartType0"]["DivergenceOfMagneticField"  ][-num_new_particles_needed:] = 0.0  # unused, I think
        f_new["PartType0"]["DivBcleaningFunctionPhi"    ][-num_new_particles_needed:] = 0.0  # unused, I think
        f_new["PartType0"]["DivBcleaningFunctionGradPhi"][-num_new_particles_needed:] = 0.0  # unused, I think

    f_old.close()
    f_new.close()
    os.rename(snapshot_file_new_tmp, snapshot_file_new)


def get_central_interpolator(hdf_file, field, k_nearest_neighbors=100):
    """ Create a linear interpolator for interpolation `field` back to the origin

    Inputs
    ------
        hdf_file : h5py.File object
            - should conform to GIZMO snapshot file standards
        field : string
            - should be a field that is defined for `PartType0` of `hdf_file`
              (e.g. "Density")

    Returns
    -------
        interp : scipy.interpolate.LinearNDInterpolator object
            - will be a 3D linear interpolator
                - inputs: coordinates *relative to box center*
                - outputs: interpolated field value (could be 1D or multi-D)


    Notes
    -----
        - May have trouble if the origin is not contained within the envelop
          of the k nearest neighbors (e.g. if the nearest neighbors are
          preferentially lopsided). In practice this hasn't been a big deal,
          since the problem is mostly symmetrical.
           - If this does become a problem, you could try to increase
             `k_nearest_neighbors` or re-write this to ensure that you are
             choosing a subset of nearest neighbors *from all octants*.

    """

    coords = hdf_file["PartType0"]["Coordinates"].value - hdf_file["Header"].attrs["BoxSize"]/2
    field_values = hdf_file["PartType0"][field].value
    dist = (coords**2).sum(axis=1)
    arg_part = np.argpartition(dist, k_nearest_neighbors)
    interpolate_from = arg_part[:k_nearest_neighbors]

    # now interpolate from those `k_nearest_neighbors` particles
    interp = interpolate.LinearNDInterpolator(coords[interpolate_from],
                                              field_values[interpolate_from])

    return interp


def create_restart_type_file(inputs_dir, restart_type, basename="RESTART"):
    if restart_type not in {entry.value for entry in RestartType}:
        raise ValueError("restart_type must be a valid value of RestartType enum; instead was {}".format(restart_type))

    with open(os.path.join(inputs_dir, basename), mode="w") as f:
        print("RESTART_TYPE={}".format(restart_type), file=f)


def snapshot_to_energy_file(snapshot_filename, energy_filename):
    """Takes a snapshot file, calculates the relevant entries for the GIZMO
    energy files, and then adds the appropriate row to the energy file.

    Note: Appends the energy information; doesn't gaurentee correct
    placement in time"""

    with h5py.File(snapshot_filename, "r") as snapshot_file:
        time = snapshot_file["Header"].attrs["Time"]

        particle_types = sorted([key
                                 for key in snapshot_file.keys()
                                 if "PartType" in key])

    N_particle_types = 6
    int_energy = np.zeros(N_particle_types, dtype=float)
    pot_energy = np.zeros(N_particle_types, dtype=float)
    kin_energy = np.zeros(N_particle_types, dtype=float)
    masses     = np.zeros(N_particle_types, dtype=float)

    for i, particle_type in enumerate(particle_types):
        masses[i] = total_mass_of_snapshot(snapshot_filename,
                                               particle_type=particle_type)

        int_energy[i] = total_internal_energy_of_snapshot(snapshot_filename,
                                                              particle_type=particle_type)

        kin_energy[i] = total_kinetic_energy_of_snapshot(snapshot_filename,
                                                              particle_type=particle_type)

        # to do:
        #  - create a way to determine the grav. potential energy
        #  - create a flag for whether a simulation is/isn't using gravity
        #  - output potential energy, conditional on that flag
        pot_energy[i] = 0

    # convert to code units
    int_energy /= M_solar * (pc / Myr)**2
    kin_energy /= M_solar * (pc / Myr)**2
    pot_energy /= M_solar * (pc / Myr)**2

    formatter = "{:<20.14f} " + " ".join(["{:18.10e}"]*27)
    with open(energy_filename, mode="a") as energy_file:
        print(formatter.format(time,
                               int_energy.sum(), pot_energy.sum(), kin_energy.sum(),
                               int_energy[0], pot_energy[0], kin_energy[0],
                               int_energy[1], pot_energy[1], kin_energy[1],
                               int_energy[2], pot_energy[2], kin_energy[2],
                               int_energy[3], pot_energy[3], kin_energy[3],
                               int_energy[4], pot_energy[4], kin_energy[4],
                               int_energy[5], pot_energy[5], kin_energy[5],
                               *masses,
                               ), file=energy_file)
