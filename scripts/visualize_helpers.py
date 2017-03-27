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


def get_snapshot_time(snapshot_filename):
    """excepts code units to be Myr for time variables"""
    with h5py.File(snapshot_filename, mode="r") as f:
        time = f["Header"].attrs["Time"]
        
    return time


def get_snapshot_times(ts):
    """`ts` should be a yt timeseries"""
    snapshot_times = np.array([ds.current_time.convert_to_cgs() 
                               for ds in ts]) / Myr
    return snapshot_times


def total_mass_of_snapshot(snapshot_filename, particle_type="PartType0"):
    with h5py.File(snapshot_filename, mode="r") as f:
        total_mass = np.sum(f[particle_type]["Masses"], dtype=float)
        
    return total_mass


def total_radial_momentum_of_snapshot(snapshot_filename, particle_type="PartType0"):
    with h5py.File(snapshot_filename, mode="r") as f:
        masses_shape = f[particle_type]["Masses"].shape
        new_shape = masses_shape + (1,)

        r_hat = f[particle_type]["Coordinates"] - (f["Header"].attrs["BoxSize"]/2)
        r_hat = r_hat / (np.sum(r_hat**2, axis=1, dtype=float).reshape(new_shape))**.5
        
        mom = np.sum(r_hat * f[particle_type]["Velocities"] \
        * np.reshape(f[particle_type]["Masses"], new_shape), dtype=float)
    
    return mom * M_solar * pc / (Myr)


def total_kinetic_energy_of_snapshot(snapshot_filename, particle_type="PartType0"):
    with h5py.File(snapshot_filename, mode="r") as f:
        masses_shape = f[particle_type]["Masses"].shape
        new_shape = masses_shape + (1,)

        E_kin = 0.5 * np.sum(np.array(f[particle_type]["Velocities"], dtype=float)**2 \
        * np.reshape(f[particle_type]["Masses"], new_shape), dtype=float)
    
    return E_kin * M_solar * (pc / Myr)**2


def total_internal_energy_of_snapshot(snapshot_filename, particle_type="PartType0"):
    with h5py.File(snapshot_filename, mode="r") as f:
        masses_shape = f[particle_type]["Masses"].shape
        new_shape = masses_shape + (1,)

        E_int = np.sum(np.array(f[particle_type]["InternalEnergy"], dtype=float) \
        * np.array(f[particle_type]["Masses"], dtype=float), dtype=float)
    
    return E_int * M_solar * (pc / Myr)**2

def total_energy_of_snapshot(snapshot_filename):
    E_tot =   total_internal_energy_of_snapshot(snapshot_filename) \
            + total_kinetic_energy_of_snapshot(snapshot_filename)
    return E_tot



def map_to_all_snapshots(outputs_dir, mapped_function):
    snapshot_filenames = get_snapshot_filenames(outputs_dir)

    results = np.array(list(map(mapped_function, snapshot_filenames)), ndmin=1)

    return results


def yt_plot_saver(plot, plot_name, plots_dir="./",
    # with_tight_layout=True,
    ):
    """ Takes a *yt plot* object, and saves it as png, eps, pdf

    This can't be done simply by passing a matplotlib figure object,
    because yt does its plotting in a way that doesn't play nice with matplotlib
    (Even if you do plot.show(), the current figure is empty) 

    
    Inputs
    ------
        plot : yt plot object, not a matplotlib figure
        plot_name : string
        plots_dir : Optional(string)

    Notes
    -----
        - assumes that any subdirectories in `plot_name` or `plots_dir` are 
          correctly specified for the host OS (e.g. using the right directory
          sepearating character such as `/` or `\`)

        - `plots_dir` could be precombined with `plot_name`. They ultimately
          get joined using an `os.path.join` call

    """

    print("Plotting: ", plot_name, "in", plots_dir)
    
    # if with_tight_layout:
        # fig.tight_layout()
    plot_filename = os.path.join(plots_dir, plot_name)
    plot.save(plot_filename + ".eps")
    plot.save(plot_filename + ".pdf")
    plot.save(plot_filename + ".png")


def mpl_plot_saver(fig, plot_name, plots_dir="./",
    with_tight_layout=True,
    ):
    """ Takes a *matplotlib figure* object, and saves it as png, eps, pdf

    This can't be used for plots created in yt.
    See `visualize_helpers.yt_plot_saver` for details
    
    Inputs
    ------
        fig : matplotlib figure object
        plot_name : string
        plots_dir : Optional(string)

    Notes
    -----
        - assumes that any subdirectories in `plot_name` or `plots_dir` are 
          correctly specified for the host OS (e.g. using the right directory
          sepearating character such as `/` or `\`)

        - `plots_dir` could be precombined with `plot_name`. They ultimately
          get joined using an `os.path.join` call

    """

    print("Plotting: ", plot_name, "in", plots_dir)
    
    if with_tight_layout:
        fig.tight_layout()
    plot_filename = os.path.join(plots_dir, plot_name)
    fig.savefig(plot_filename + ".eps")
    fig.savefig(plot_filename + ".pdf")
    fig.savefig(plot_filename + ".png")

