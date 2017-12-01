import h5py
import numpy as np

from units import M_solar, m_proton, pc, yr, Myr, km, s, gamma


def total_mass_of_snapshot(snapshot_filename, particle_type="PartType0"):
    with h5py.File(snapshot_filename, mode="r") as f:
        total_mass = np.sum(f[particle_type]["Masses"], dtype=float)

    return total_mass


def total_radial_momentum_of_snapshot(snapshot_filename,
                                      particle_type="PartType0"):
    with h5py.File(snapshot_filename, mode="r") as f:
        masses_shape = f[particle_type]["Masses"].shape
        new_shape = masses_shape + (1,)

        r_hat = f[particle_type]["Coordinates"] \
            - (f["Header"].attrs["BoxSize"] / 2)
        r_hat = r_hat \
            / (np.sum(r_hat**2, axis=1, dtype=float).reshape(new_shape))**.5

        mom = np.sum(r_hat * f[particle_type]["Velocities"] *
                     np.reshape(f[particle_type]["Masses"], new_shape),
                     dtype=float)

    return mom * M_solar * pc / (Myr)


def total_kinetic_energy_of_snapshot(snapshot_filename,
                                     particle_type="PartType0"):
    with h5py.File(snapshot_filename, mode="r") as f:
        masses_shape = f[particle_type]["Masses"].shape
        new_shape = masses_shape + (1,)

        E_kin = 0.5 * np.sum(
            np.array(f[particle_type]["Velocities"], dtype=float)**2 *
            np.reshape(f[particle_type]["Masses"], new_shape),
            dtype=float)

    return E_kin * M_solar * (pc / Myr)**2


def total_internal_energy_of_snapshot(snapshot_filename,
                                      particle_type="PartType0"):
    with h5py.File(snapshot_filename, mode="r") as f:

        E_int = np.sum(
            np.array(f[particle_type]["InternalEnergy"], dtype=float) *
            np.array(f[particle_type]["Masses"], dtype=float),
            dtype=float)

    return E_int * M_solar * (pc / Myr)**2


def total_magnetic_energy_of_snapshot(snapshot_filename,
                                      particle_type="PartType0"):
    with h5py.File(snapshot_filename, mode="r") as f:
        if "MagneticField" not in f[particle_type].keys():
            return 0

        masses_shape = f[particle_type]["Masses"].shape
        new_shape = masses_shape + (1,)

        E_mag = np.sum(
            np.array(f[particle_type]["MagneticField"], dtype=float)**2 *
            np.reshape(f[particle_type]["Masses"], new_shape) /
            np.reshape(f[particle_type]["Density"], new_shape),
            dtype=float) / (8 * np.pi)

    return E_mag * (pc)**3 * (1e-6)**2


def total_energy_of_snapshot(snapshot_filename):
    E_tot = total_internal_energy_of_snapshot(snapshot_filename) \
            + total_magnetic_energy_of_snapshot(snapshot_filename) \
            + total_kinetic_energy_of_snapshot(snapshot_filename)
    return E_tot
