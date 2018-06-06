import glob
import os
import h5py
import numpy as np
import pandas as pd
import yt
import seaborn as sns
import pickle
from matplotlib import pyplot as plt
import warnings

from units import M_solar, m_proton, pc, yr, Myr, gamma, km, s
import MHD
from grackle_helpers import \
    wrapped_initializer, \
    temperature, \
    cooling_rate, \
    FluidContainer

from snapshot_helpers import \
    total_mass_of_snapshot, \
    total_radial_momentum_of_snapshot, \
    total_kinetic_energy_of_snapshot, \
    total_internal_energy_of_snapshot, \
    total_energy_of_snapshot

#  # Boilerplate path hack to give access to full clustered_SNe package
import sys
import os
if os.pardir not in sys.path[0]:
    file_dir = os.getcwd()
    sys.path.insert(0, os.path.join(file_dir,
                                    os.pardir,
                                    os.pardir))

from clustered_SNe.analysis.parse import RunSummary

# ############# YT HELPERS ############


@yt.derived_field(("gas", "pressure"),
                  units="auto",
                  force_override=True,
                  dimensions=yt.units.dimensions.pressure,
                  sampling_type="cell")
def _pressure(field, data):
    return (gamma-1) * data[("gas", "thermal_energy")] * \
     data[("gas", "density")]


@yt.particle_filter(requires=["particle_velocity_magnitude"],
                    filtered_type='all')
def affected(pfilter, data):
    filter = data[(pfilter.filtered_type,
                  "particle_velocity_magnitude")].to("m/s") > 1
    return filter

# # Grackle-required fields
# # Make sure these units are correct!!!
rval, my_chemistry_particles = wrapped_initializer(
    float(yt.units.Msun.to("g").value),
    float(yt.units.pc.to("cm").value),
    float(yt.units.Myr.to("s").value),
    verbose=True,
    )
assert(rval == 1)

rval, my_chemistry_1D = wrapped_initializer(
    m_proton,  # grackle breaks if you try to use cgs units -_-
    1,
    1,
    verbose=True,
    )
assert(rval == 1)

# # bind `my_chemistry_particles` as a closure to `temperature`
# # Note: yt has problems with using lambdas; better to use a real function
def _temperature(field, data):
    return temperature(field, data, my_chemistry_particles)

yt.add_field(("all", "temperature"),
             function=_temperature,
             sampling_type="particle",
             units="K",
             force_override=True,
             )

# # bind `my_chemistry_particles` as a closure to `cooling_rate`
# # Note: yt has problems with using lambdas; better to use a real function
def _cooling_rate(field, data):
    return cooling_rate(field, data, my_chemistry_particles)
yt.add_field(("all", "cooling_rate"),
             function=_cooling_rate,
             sampling_type="particle",
             units="Msun * pc**2 / Myr**3 ",
             force_override=True,
             )

# #########################

snapshot_filename_format = "snapshot_???.hdf5"


def snapshot_filename_to_number(snapshot_filename):
    snapshot_number = int(os.path.basename(snapshot_filename).replace(".hdf5", "").replace("snapshot_", ""))
    return snapshot_number


def get_dirs(run_name):
    run_dir = os.path.join(os.path.pardir, "runs", run_name)

    inputs_dir  = os.path.join(run_dir, "inputs")
    outputs_dir = os.path.join(run_dir, "outputs")

    return inputs_dir, outputs_dir


def get_snapshot_filenames(outputs_dir,
                           snapshot_filename_format=snapshot_filename_format):

    return sorted(glob.glob(os.path.join(outputs_dir,
                                         snapshot_filename_format)))


def load_snapshots(outputs_dir):
    unit_base = {
        "UnitLength_in_cm": (pc),
        "UnitVelocity_in_cm_per_s": (pc / Myr),
        "UnitMass_in_g": (M_solar)
    }

    snapshot_filenames = get_snapshot_filenames(outputs_dir)

    n_files_ready = len(snapshot_filenames)
    if n_files_ready == 0:
        raise FileNotFoundError("No snapshots found in {}".format(outputs_dir))

    ts = yt.load(os.path.join(outputs_dir, snapshot_filename_format),
                 unit_base=unit_base)
    return ts


def load_ds_from_ts(ts, i):
    """gets the `ith` snapshot from `ts`.
    Allows automated binding of filters, etc
    """
    ds = ts[i]
    ds.add_particle_filter("affected")
    return ds


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


def map_to_all_snapshots(outputs_dir, mapped_function):
    snapshot_filenames = get_snapshot_filenames(outputs_dir)

    results = np.array(list(map(mapped_function, snapshot_filenames)), ndmin=1)

    return results


def yt_plot_saver(plot, plot_name, plots_dir="./"):
    """ Takes a *yt plot* object, and saves it as png, eps, pdf

    This can't be done simply by passing a matplotlib figure object,
    because yt's plotting doesn't play nice with matplotlib
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

    plot_filename = os.path.join(plots_dir, plot_name)
    plot.save(plot_filename + ".eps")
    plot.save(plot_filename + ".pdf")
    plot.save(plot_filename + ".png")


def mpl_plot_saver(fig, plot_name, plots_dir="./",
                   with_tight_layout=True):
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


# ####### High level plotting routines ########


def plot_projected_density(ts, snapshot_number, snapshot_number_to_index_map,
                           SN_times, plots_dir,
                           save_plot=True,
                           show_plot=True,
                           seaborn_style="white"):
    """ Creates [and optionally saves] a plot of projected density

    Inputs
    ------
    ts : yt.data_objects.time_series.DatasetSeries object
        - contains all the uncompressed snapshots in a directory
    snapshot_number : int
         - the specific snapshot you want to plot
    snpashot_number_to_index_map : dict (int -> int)
        - maps snapshot_number to an index of `ts`
    SN_times : np.ndarray (dtype: float)
        - the times that SNe occur for this simulation (in Myr)
    plots_dir : str
        - where plots should be saved, if they are to be saved
    save_plot : Optional(bool)
        - True if you want to save plot images
    show_plot : Optional(bool)
        - True if you want to show plots
         (you might want to disable this for non-interactive sessions)

    """

    i = snapshot_number_to_index_map[snapshot_number]
    ds = load_ds_from_ts(ts, i)

    p = yt.ProjectionPlot(ds, "x", ("gas", "density"))
    p.set_log('density', False)
    # p.set_cmap(field="density", cmap="viridis")
    p.annotate_timestamp(corner="upper_left", draw_inset_box=True)

    t = ds.current_time.convert_to_cgs().value / Myr
    N_SNe_so_far = np.sum(t > SN_times)
    p.annotate_text((.8, .94),
                    "$N_\mathrm{{SNe}}$: {}".format(N_SNe_so_far),
                    coord_system="axis",
                    inset_box_args={"facecolor": "darkslategray",
                                    "alpha": 0.9},
                    )
    if show_plot:
        with sns.axes_style(seaborn_style):
            p.show()

    if save_plot:
        subdir = "projected"
        if not os.path.exists(os.path.join(plots_dir, subdir)):
            os.mkdir(os.path.join(plots_dir, subdir))
        plot_name = os.path.join(subdir, "{}_snapshot_{:0>3}".format("density", snapshot_number))

        yt_plot_saver(p, plot_name, plots_dir)

    return p


field_type = {
    "density": "gas",
    "temperature": "gas",
    "pressure": "gas",
    "velocity_magnitude": "gas",
    "radius": "index",
    "metallicity": "gas"
}


def plot_sliced_field(ts, snapshot_number, snapshot_number_to_index_map, field,
                      SN_times, plots_dir,
                      save_plot=True,
                      show_plot=True,
                      seaborn_style="white",
                      add_magnetic_field_lines=False,
                      font_size=None,
                      SlicePlot_kwargs=dict()):
    """ Creates [and optionally saves] a plot of `field`, sliced at x=0

    Inputs
    ------
    ts : yt.data_objects.time_series.DatasetSeries object
        - contains all the uncompressed snapshots in a directory
    snapshot_number : int
        - the specific snapshot you want to plot
    snpashot_number_to_index_map : dict (int -> int)
        - maps snapshot_number to an index of `ts`
    field : str
        - a key of `field_type`
    SN_times : np.ndarray (dtype: float)
        - the times that SNe occur for this simulation (in Myr)
    plots_dir : str
        - where plots should be saved, if they are to be saved
    save_plot : Optional(bool)
        - True if you want to save plot images
    show_plot : Optional(bool)
        - True if you want to show plots
         (you might want to disable this for non-interactive sessions)
    seaborn_style : Optional(str)
        - a valid argument for `seaborn.axes_style`
    add_magnetic_field_lines : Optional(bool)
        - True if you want to try adding field lines to the plot
          (field lines must already be computed and saved to disk
           with filenames determined by `MHD.get_field_lines_filename_from_ds`)
    font_size : Optional(None or int)
        - font size that yt should use
    SlicePlot_kwargs : Optional(dict)
        - keyword arguments to be passed to `yt.SlicePlot`
    """
    i = snapshot_number_to_index_map[snapshot_number]
    ds = load_ds_from_ts(ts, i)
    s = yt.SlicePlot(ds, "x", (field_type[field], field), **SlicePlot_kwargs)

    if field == "density":
        s.set_unit(field, "amu/cm**3")
        # s.set_zlim(field, 1e-3, 1e1)
        s.set_colorbar_label(field, "Density $ [ m_\mathrm{H} \; \mathrm{cm}^{-3} ] $")

    # s.set_cmap(field=field, cmap="viridis")
    s.annotate_timestamp(corner="upper_left", draw_inset_box=True,
                         time_format="t = {time:.2f} {units}",
                         time_unit="Myr")
    t = ds.current_time.convert_to_cgs().value / Myr
    N_SNe_so_far = np.sum(t > SN_times)
    s.annotate_text((.8, .94),
                    "$N_\mathrm{{SNe}}$: {}".format(N_SNe_so_far),
                    coord_system="axis",
                    inset_box_args={"facecolor": "darkslategray",
                                    "alpha": 0.9},
                    )

    if font_size is not None:
        s.set_font({"size": font_size})

    if show_plot:
        with sns.axes_style(seaborn_style):
            s.show()

    if save_plot:
        subdir = "slice"
        if not os.path.exists(os.path.join(plots_dir, subdir)):
            os.mkdir(os.path.join(plots_dir, subdir))
        plot_name = os.path.join(subdir, "{}_snapshot_{:0>3}".format(field, snapshot_number))

        yt_plot_saver(s, plot_name, plots_dir)

    if add_magnetic_field_lines:
        filename = MHD.get_field_lines_filename_from_ds(ds, plots_dir)
        if os.path.exists(filename):

            df_lines = pd.read_csv(filename)
            p = s.plots[(field_type[field], field)]

            xlim = p.axes.get_xlim()
            unit_scaling = s.width[0].value / [xlim[1] - xlim[0]]

            MHD.plot_field_lines(df_lines, p.axes,
                                 unit_scaling=unit_scaling,
                                 replace_nans=False)

            if show_plot:
                with sns.axes_style(seaborn_style):
                    s.show()
            if save_plot:
                subdir = "slice"
                if not os.path.exists(os.path.join(plots_dir, subdir)):
                    os.mkdir(os.path.join(plots_dir, subdir))
                plot_name = os.path.join(subdir, "{}_snapshot_{:0>3}_with_field_lines".format(field, snapshot_number))
                yt_plot_saver(s, plot_name, plots_dir)

    return s


def get_pickle_filename_profile(run_name, snapshot_number, weight_field,
                                save_dir):
    pickle_filename = os.path.join(save_dir,
                                   "{}-snapshot_{:03}-{}.pickle".format(run_name,
                                                               snapshot_number,
                                                               weight_field))

    return pickle_filename


def is_weight_field_processed(run_name, snapshot_number, weight_field,
                              save_dir):
    """For use filtering `save_phase_diagram_data_3D`"""
    pickle_filename = get_pickle_filename_profile(run_name,
                                                  snapshot_number,
                                                  weight_field,
                                                  save_dir,
                                                  )
    return os.path.exists(pickle_filename)


def save_phase_diagram_data_3D(run_names,
                               snapshot_numbers,
                               weight_fields,
                               weight_units,
                               save_dir,
                               n_bins=450,
                               xlim_log=(-4, 3),
                               ylim_log=(1, 7),
                               verbose=True,
                               overwrite=False,
                               ):
    """Reads in a number of snapshots, and makes any phase diagrams that are
    missing from `save_dir`. Assumes the x-axis is density and the y-axis
    is temperature. Z-axis set by `weight_fields`

    Inputs
    ------
    run_names: iterable of strings
    snapshot_numbers: iterable of integers
        assumes these snapshots exist as uncompressed snapshots
    weight_fields: iterable of strings
        assumes a field_type of "affected"
    weight_units: dict
        must have keys for all elements of weight_fields
    save_dir: str formatted as a directory path
    n_bins: Optional(int)
    (x|y)lim_log: Optional(tuple of 2 numeric)
        the limits of your phase diagram, in log10 values.
        any data outside these limits are dropped from the phase diagram
    verbose: Optional(bool)
    overwrite: Optional(bool)
        overwrites existing save files if True

    Outputs
    -------
    None

    Side Effects
    ------------
    May save a pickled data file for future use.

    Notes
    -----
    Only checks for a file matching the phase diagram data name convention,
    set by `get_pickle_filename_profile`. Doesn't check if the data is correct,
    or using the same binning parameters.

    Be careful about choosing n_bins! If you want a lot of options with
    rebinning, then you'll need to chose a number with a diversity of prime
    factors. If you just choose a power of 2, you'll only be able to
    halve/double the resolution, which doesn't give you many options.
    """
    x_bins = np.logspace(*xlim_log, num=n_bins + 1)
    y_bins = np.logspace(*ylim_log, num=n_bins + 1)

    for run_name in run_names:
        for snapshot_number in snapshot_numbers:
            mask = np.array([not is_weight_field_processed(run_name,
                                                           snapshot_number,
                                                           weight_field,
                                                           save_dir)
                             for weight_field in weight_fields])
            if overwrite:
                weight_fields_remaining = weight_fields
            else:
                weight_fields_remaining = weight_fields[mask]
                if verbose:
                    print("skipping already-processed fields {} for snapshot {}, run {}".format(
                        weight_fields[~mask],
                        snapshot_number,
                        run_name,
                        ))
                if weight_fields_remaining.size == 0:
                        continue

            inputs_dir, outputs_dir = get_dirs(run_name)
            ts = load_snapshots(outputs_dir)
            snapshot_filenames = get_snapshot_filenames(outputs_dir)

            snapshot_number_to_index_map = {snapshot_filename_to_number(filename): i
                                            for i, filename in enumerate(snapshot_filenames)}

            uncompressed_snapshot_numbers = sorted(list(snapshot_number_to_index_map.keys()))

            ds = load_ds_from_ts(ts, snapshot_number_to_index_map[snapshot_number])
            dd = ds.all_data()

            x_in = dd["affected", "density"].to("amu / cm**3").value
            y_in = dd["affected", "temperature"].to("K").value

            for weight_field in weight_fields_remaining:
                pickle_filename = get_pickle_filename_profile(run_name,
                                                              snapshot_number,
                                                              weight_field,
                                                              save_dir)

                if verbose:
                    print("creating save file: {}".format(pickle_filename))

                Z, _, _ = np.histogram2d(
                    x_in, y_in,
                    weights=dd["affected", weight_field].to(weight_units[weight_field]).value,
                    bins=[x_bins, y_bins],
                )

                with open(pickle_filename, "wb") as pickle_file:
                    pickle.dump((x_bins, y_bins, Z), pickle_file)


def grackle_temperature_1D(df):
    """For use with `save_phase_diagram_data_3D`."""
    fc = FluidContainer(my_chemistry_1D, df.shape[0])

    fc["density"][:]    = df.Density.values / m_proton
    fc["metal"][:]      = df.Density.values * df.Z.values / m_proton
    fc["x-velocity"][:] = df.Velocity.values
    fc["y-velocity"][:] = 0
    fc["z-velocity"][:] = 0
    gamma = 5/3
    specific_internal_energy = (1 / (gamma - 1)) * df.Pressure.values \
                                                 / df.Density.values
    fc["energy"][:] = specific_internal_energy

    fc.calculate_temperature()

    return fc["temperature"]


def grackle_cooling_rate_1D(df):
    """For use with `save_phase_diagram_data_3D`."""
    fc = FluidContainer(my_chemistry_1D, df.shape[0])

    fc["density"][:]    = df.Density.values / m_proton
    fc["metal"][:]      = df.Density.values * df.Z.values / m_proton
    fc["x-velocity"][:] = df.Velocity.values
    fc["y-velocity"][:] = 0
    fc["z-velocity"][:] = 0
    gamma = 5/3
    specific_internal_energy = (1 / (gamma - 1)) * df.Pressure.values \
                                                 / df.Density.values
    fc["energy"][:] = specific_internal_energy

    fc.calculate_temperature()
    fc.calculate_cooling_time()

    dt = 1  # [second]
    fc.solve_chemistry(dt)

    de = (specific_internal_energy - fc["energy"])

    de_dt = de / dt

    dE_dt = de_dt * df.Mass.values * M_solar

    return dE_dt


def get_closest_1D_snapshot(time,
                            data_dir_1D="1D_data-high_time_res/",
                            data_id_1D="F5509BF1-3F9E-4008-B795-0482ECED199B",
                            verbose=False, very_verbose=False):
    """`time` should be in Myr, with t=0 as the first SN (like the 3D runs)."""
    run_summary = RunSummary(data_dir_1D, data_id_1D)
    corrected_times = run_summary.times - run_summary.times[0]
    corrected_times /= Myr

    if verbose:
        if very_verbose:
            print("corrected_times: ", corrected_times)

    i = np.argmin(np.abs(corrected_times - time))
    if verbose:
        print("using snapshot: {} at time {:.3f} Myr".format(i, corrected_times[i]))

    df_1D_snapshot = run_summary.df.loc[i]

    df_1D_snapshot["grackle_Temperature"] = grackle_temperature_1D(df_1D_snapshot)
    df_1D_snapshot["grackle_cooling_rate"] = grackle_cooling_rate_1D(df_1D_snapshot)

    return df_1D_snapshot


def save_phase_diagram_data_1D(snapshot_numbers,
                               weight_fields,
                               save_dir,
                               sample_run_name_3D,
                               n_bins=450,
                               xlim_log=(-4, 3),
                               ylim_log=(1, 7),
                               verbose=True,
                               get_closest_1D_snapshot_kwargs=dict(),
                               ):
    """1D analog of `save_phase_diagram_data_1D`.

    Minor differences:
     - Requires argument `sample_run_name_3D` (string). This is only needed
        for getting the snapshot times; actual snapshot data doesn't matter.
        the snapshots of `snapshot_numbers` should be uncompressed.
     - No need for `run_names`. Input data found using
        `get_closest_1D_snapshot_kwargs`. Output data saved using a dummy
        run_name of `1D`
     - weight_units not accepted. Assumes:
            particle_mass -> Msun
            cooling_rate -> erg / Myr
       Sorry that it's not more extensible =/
     - weight_fields elements can only be `particle_mass` or `cooling_rate`
       since they need to be hand-coded, unlike the 3D function

    """
    x_bins = np.logspace(*xlim_log, num=n_bins+1)
    y_bins = np.logspace(*ylim_log, num=n_bins+1)

    inputs_dir, outputs_dir = get_dirs(sample_run_name_3D)
    ts = load_snapshots(outputs_dir)
    snapshot_filenames = get_snapshot_filenames(outputs_dir)

    snapshot_number_to_index_map = {snapshot_filename_to_number(filename): i
                                    for i, filename in enumerate(snapshot_filenames)}

    uncompressed_snapshot_numbers = sorted(list(snapshot_number_to_index_map.keys()))

    for snapshot_number in snapshot_numbers:
        current_time = load_ds_from_ts(ts, snapshot_number_to_index_map[snapshot_number]
                                       ).current_time.value
        current_time = float(current_time)
        if verbose:
            print(snapshot_number, current_time)

        df_1D_snapshot = get_closest_1D_snapshot(
            current_time, verbose=True, **get_closest_1D_snapshot_kwargs)

        affected = (np.abs(df_1D_snapshot.Velocity) > 100)

        for weight_field in weight_fields:
            if weight_field == "particle_mass":
                Z, _, _ = np.histogram2d(
                    df_1D_snapshot["Density"][affected] / (m_proton),
                    df_1D_snapshot["grackle_Temperature"][affected],
                    weights=df_1D_snapshot["Mass"][affected],
                    bins=[x_bins, y_bins],
                )
            elif weight_field == "cooling_rate":
                # grackle incorrectly zero's out very low cooling rates :(
                # so I'm changing it back to small, but non-zero, so that it
                # won't appear "missing" in the plot
                weights = df_1D_snapshot["grackle_cooling_rate"][affected] * Myr
                weights[weights < .1] = .1
                Z, _, _ = np.histogram2d(
                    df_1D_snapshot["Density"][affected] / (m_proton),
                    df_1D_snapshot["grackle_Temperature"][affected],
                    weights=weights,
                    bins=[x_bins, y_bins],
                )
            else:
                raise NotImplementedError("Weight field {} not implemented in 1D".format(weight_field))

            pickle_filename = get_pickle_filename_profile("1D",
                                                          snapshot_number,
                                                          weight_field,
                                                          save_dir)
            with open(pickle_filename, mode="wb") as pf:
                pickle.dump((x_bins, y_bins, Z), pf)


def plot_phase_diagram(ts, snapshot_number, snapshot_number_to_index_map,
                       SN_times, plots_dir,
                       weight_field,
                       field_type="all",
                       save_plot=True,
                       show_plot=True,
                       seaborn_style="ticks",
                       bins=50):
    """Create [and optionally save] a density-temperature phase diagram.

    Inputs
    ------
    ts : yt.data_objects.time_series.DatasetSeries object
        - contains all the uncompressed snapshots in a directory
    snapshot_number : int
        - the specific snapshot you want to plot
    snpashot_number_to_index_map : dict (int -> int)
        - maps snapshot_number to an index of `ts`
    SN_times : np.ndarray (dtype: float)
        - the times that SNe occur for this simulation (in Myr)
    plots_dir : str
        - where plots should be saved, if they are to be saved
    weight_field : str
        - the field you'll be creating your 2D histogram of
          (e.g. `particle_mass`)
    field_type : Optional(str)
        - the yt field type (e.g. `all`, `deposit` or your own particle filter)
    save_plot : Optional(bool)
        - True if you want to save plot images
    show_plot : Optional(bool)
        - True if you want to show plots
         (you might want to disable this for non-interactive sessions)
    seaborn_style : Optional(str)
        - a valid argument for `seaborn.axes_style`
    bins : Optional(int)
        - number of bins in the x and y axes

    Returns
    -------
    ppp

    """
    i = snapshot_number_to_index_map[snapshot_number]
    ds = load_ds_from_ts(ts, i)
    dd = ds.all_data()

    ppp = yt.ParticlePhasePlot(dd,
                               (field_type, "density"),
                               (field_type, "temperature"),
                               z_fields=(field_type, weight_field),
                               x_bins=bins,
                               y_bins=bins,)

    if weight_field == "particle_mass":
        ppp.profile.set_field_unit("particle_mass", "Msun")
        ppp.plots[(field_type, "particle_mass")].zmin, ppp.plots[(field_type, "particle_mass")].zmax = (None, None)

    elif weight_field == "cooling_rate":
        ppp.profile.set_field_unit("cooling_rate", "erg / Myr")
        ppp.plots[(field_type, "cooling_rate")].zmin, ppp.plots[(field_type, "cooling_rate")].zmax = (None, None)

    ppp.profile.set_x_unit("amu/cm**3")
    ppp.set_xlabel("Density $ [ m_\mathrm{H} \; \mathrm{cm}^{-3} ] $")

    ppp.set_xlim(1e-4, 1e3)
    ppp.set_ylim(1e1, 1e7)

    ppp.set_log("density", True)
    ppp.set_log("temperature", True)
    ppp.set_log(weight_field, True)

    ppp.set_background_color((field_type, weight_field), color="white")

    if show_plot:
        with sns.axes_style(seaborn_style):
            ppp.show()

    if save_plot:
        subdir = "phase_plot"
        if not os.path.exists(os.path.join(plots_dir, subdir)):
            os.mkdir(os.path.join(plots_dir, subdir))
        plot_name = os.path.join(subdir, "snapshot_{:0>3}_{}-{}".format(snapshot_number, weight_field, field_type))
        yt_plot_saver(ppp, plot_name, plots_dir)

    return ppp

# ### Profile plotting (requires a bit more low-level programming)


def create_density_profile(ds, n_bins=20):
    dd = ds.all_data()
    r_max = ds.domain_width[0]/2

    dr = r_max / n_bins

    rs = np.linspace(0, r_max.value, num=n_bins+1)[1:]

    dmass = np.zeros(n_bins)
    ones = np.zeros(n_bins, dtype=int)

    for i in range(n_bins):
        r_i = dr*(i)
        r_o = dr*(i+1)

        mask =    (dd["all", "particle_position_spherical_radius"] >= r_i) \
                & (dd["all", "particle_position_spherical_radius"] <  r_o)

#         ones[i] = mask.sum()
        dmass[i] = dd["all", "Masses"][mask].sum().convert_to_cgs().value

    Vs = 4/3*np.pi*rs**3
    Vs = np.insert(Vs, 0, 0)
    dVs = Vs[1:] - Vs[:-1]

    densities = dmass / (dVs * pc**3)

    return rs, densities


field_y_labels = {
    "density": r"$\rho$ $[\mathrm{m_p}$ $\mathrm{cm}^{-3}]$",
    "temperature": r"$T$ $[\mathrm{K}]$",
    "pressure": r"$P$ $[\mathrm{ergs}$ $\mathrm{cm}^{-3}]$",
    "velocity_magnitude": r"$\|\mathbf{v}\|$ $[\mathrm{km}$ $\mathrm{s}^{-1}]$",
    "radial_velocity": r"$v_r$ $[\mathrm{km}$ $\mathrm{s}^{-1}]$",
    "Metallicity": r"$Z / Z_\odot$",
}

field_weight = {
    "temperature": "cell_mass",
    "pressure": "cell_volume",
    "velocity_magnitude": "cell_mass",
    "radial_velocity": "cell_mass",
    "Metallicity": "cell_mass",
}

field_units = {
    "density": m_proton,
    "temperature": 1,
    "pressure": 1,
    "velocity_magnitude": km / s,
    "radial_velocity": km / s,
    "Metallicity": 0.02,
}


def format_time_str(current_time):
    """`current_time` should be in Myr."""
    if current_time < 0:
        raise RuntimeError("Invalid time: {}".format(current_time))
    elif current_time < 1e-3:
        time = current_time
        time_units = "Myr"
        time_str = r"$t$ $= {:.1e}$ $\mathrm{{{}}}$".format(time, time_units)
    elif current_time < 1:
        time = current_time * 1e3
        time_units = "kyr"
        time_str = r"$t$ $= {:.0f}$ $\mathrm{{{}}}$".format(time, time_units)
    elif current_time < 10:
        time = current_time
        time_units = "Myr"
        time_str = r"$t$ $= {:.1f}$ $\mathrm{{{}}}$".format(time, time_units)
    else:
        time = current_time
        time_units = "Myr"
        time_str = r"$t$ $= {:.0f}$ $\mathrm{{{}}}$".format(time, time_units)

    return time_str


def plot_profile(ts, snapshot_number, snapshot_number_to_index_map, field,
                 rho_0, plots_dir,
                 save_plot=True,
                 show_plot=True):
    """Create [and optionally save] a radial profile plot of `field`

    Inputs
    ------
    ts : yt.data_objects.time_series.DatasetSeries object
        - contains all the uncompressed snapshots in a directory
    snapshot_number : int
        - the specific snapshot you want to plot
    snpashot_number_to_index_map : dict (int -> int)
        - maps snapshot_number to an index of `ts`
    field : str
        - a key of `field_type`
    rho_0 : float
        - background density in m_p / cm**3
        - adds a horizontal black dashed line at this value.
          If you don't want this value, set rho_0 outside of [1e-4, 1e1], e.g. 0
    plots_dir : str
        - where plots should be saved, if they are to be saved
    save_plot : Optional(bool)
        - True if you want to save plot images
    show_plot : Optional(bool)
        - True if you want to show plots
         (you might want to disable this for non-interactive sessions)

    """
    i = snapshot_number_to_index_map[snapshot_number]
    ds = load_ds_from_ts(ts, i)
    sp = ds.sphere(ds.domain_center, ds.domain_width[0]/2)

    if field is "density":
        rs, densities = create_density_profile(ds, n_bins=64)
        plt.plot(rs, densities / field_units[field])

        plt.ylim(ymin=1e-4)

        plt.axhline(rho_0, linestyle="dashed", color="k")

    else:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pp = yt.create_profile(sp,
                                   "radius", [field, "ones"],
                                   weight_field=field_weight[field],
                                   units={"radius": "pc"},
                                   logs={"radius": False},
                                   n_bins=64,
                                   )
        mask = pp["ones"] > 0.1  # filter out bins with no particles
        plt.plot(pp.x.value[mask],
                 pp[field][mask] / field_units[field])

    plt.yscale("log")
    plt.ylabel(field_y_labels[field])

    plt.xlabel(r"$R$ $[\mathrm{pc}]$")

    current_time = ds.current_time.convert_to_cgs().value / Myr
    title = format_time_str(current_time)

    plt.title(title)

    if save_plot:
        subdir = "profile"
        if not os.path.exists(os.path.join(plots_dir, subdir)):
            os.mkdir(os.path.join(plots_dir, subdir))
        plot_name = os.path.join(subdir, "{}_snapshot_{:0>3}".format(field, snapshot_number))
        mpl_plot_saver(plt.gcf(), plot_name, plots_dir)
