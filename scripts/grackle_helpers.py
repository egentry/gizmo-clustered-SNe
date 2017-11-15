from __future__ import print_function, division

import os

from pygrackle import chemistry_data, FluidContainer

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    sec_per_Myr, \
    cm_per_mpc

from units import M_solar, pc, Myr

import yt


def wrapped_initializer(mass_unit_cgs, length_unit_cgs, time_unit_cgs,
                        verbose=False):
    """ Create a `grackle.chemistry_data` object

    Inputs
    ------
    mass_units_cgs : float
    length_unit_cgs : float
    time_unit_cgs : float
    verbose : Optional(bool)

    Returns
    -------
    rval : int
        return value from Grackle. Returns 1 for success, 0 for failure.
    my_chemistry : `grackle.chemistry_data` object

    """
    if verbose:
        print("""`grackle_helpers.wrapped_initializer` is setting up grackle assuming code units:
    mass   : {:.3} M_solar
    length : {:.3} pc
    time   : {:.3} Myr
""".format(
          mass_unit_cgs / M_solar,
          length_unit_cgs / pc,
          time_unit_cgs / Myr,
        ))

    current_redshift = 0.

    # Set solver parameters
    my_chemistry = chemistry_data()
    my_chemistry.use_grackle = 1
    my_chemistry.with_radiative_cooling = 1
    my_chemistry.primordial_chemistry = 0
    my_chemistry.metal_cooling = 1
    my_chemistry.UVbackground = 1
    my_chemistry.self_shielding_method = 0
    my_chemistry.H2_self_shielding = 0

    grackle_dir = "/pfs/home/egentry/local/grackle"
    my_chemistry.grackle_data_file = os.path.join(
        grackle_dir, "input", "CloudyData_UVB=HM2012.h5")

    print("grackle cooling file: ", my_chemistry.grackle_data_file)

    # Set units
    my_chemistry.comoving_coordinates = 0  # proper units
    my_chemistry.a_units = 1.0
    my_chemistry.a_value = 1. / (1. + current_redshift) / \
        my_chemistry.a_units
    my_chemistry.density_units = mass_unit_cgs / length_unit_cgs**3
    my_chemistry.length_units = length_unit_cgs
    my_chemistry.time_units = time_unit_cgs
    my_chemistry.velocity_units = my_chemistry.a_units * \
        (my_chemistry.length_units / my_chemistry.a_value) / \
        my_chemistry.time_units

    rval = my_chemistry.initialize()

    return rval, my_chemistry


def temperature(field, data, my_chemistry):

    fc = FluidContainer(my_chemistry,
                        data["PartType0", "Masses"].shape[0])

    fc["density"][:] = data["PartType0", "Density"]
    fc["metal"][:] = data["PartType0", "Density"] * data["PartType0",
                                                         "Metallicity"]
    fc["x-velocity"][:] = data["PartType0", "Velocities"][:, 0]
    fc["y-velocity"][:] = data["PartType0", "Velocities"][:, 1]
    fc["z-velocity"][:] = data["PartType0", "Velocities"][:, 2]
    fc["energy"][:] = data["PartType0", "InternalEnergy"]

    fc.calculate_temperature()

    return fc["temperature"] * yt.units.Kelvin


def cooling_rate(field, data, my_chemistry):

    fc = FluidContainer(my_chemistry, data["PartType0", "Masses"].shape[0])

    fc["density"][:] = data["PartType0", "Density"]
    fc["metal"][:] = data["PartType0", "Density"] * data["PartType0", "Metallicity"]
    fc["x-velocity"][:] = data["PartType0", "Velocities"][:,0]
    fc["y-velocity"][:] = data["PartType0", "Velocities"][:,1]
    fc["z-velocity"][:] = data["PartType0", "Velocities"][:,2]
    fc["energy"][:] = data["PartType0", "InternalEnergy"]

    fc.calculate_temperature()
    fc.calculate_cooling_time()

    dt = yt.units.s
    fc.solve_chemistry(dt.value / my_chemistry.time_units)

    de = (data["PartType0", "InternalEnergy"].value - fc["energy"])
    de *= yt.units.Msun * yt.units.pc**2 / yt.units.Myr**2

    de_dt = de / dt

    return de_dt
