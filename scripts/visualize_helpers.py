import glob
import os
import h5py
import numpy as np
import pandas as pd
import yt
import seaborn as sns

from units import M_solar, m_proton, pc, yr, Myr, gamma, km, s
import MHD
from grackle_helpers import wrapped_initializer, temperature, cooling_rate

from snapshot_helpers import \
    total_mass_of_snapshot, \
    total_radial_momentum_of_snapshot, \
    total_kinetic_energy_of_snapshot, \
    total_internal_energy_of_snapshot, \
    total_energy_of_snapshot

from matplotlib import pyplot as plt

import warnings


# ############# YT HELPERS ############

@yt.derived_field(name="pressure", units="g  / s**2 / cm",
                  sampling_type="particle")
def _pressure(field, data):
    return (gamma-1) * data["thermal_energy"] * data["density"]


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
    p.set_cmap(field="density", cmap="viridis")
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
                      add_magnetic_field_lines=False):
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

    """
    i = snapshot_number_to_index_map[snapshot_number]
    ds = load_ds_from_ts(ts, i)
    s = yt.SlicePlot(ds, "x", (field_type[field], field))

    if field == "density":
        s.set_unit(field, "amu/cm**3")
        # s.set_zlim(field, 1e-3, 1e1)
        s.set_colorbar_label(field, "Density $ [ m_\mathrm{H} \; \mathrm{cm}^{-3} ] $")

    s.set_cmap(field=field, cmap="viridis")
    s.annotate_timestamp(corner="upper_left", draw_inset_box=True)
    t = ds.current_time.convert_to_cgs().value / Myr
    N_SNe_so_far = np.sum(t > SN_times)
    s.annotate_text((.8, .94),
                    "$N_\mathrm{{SNe}}$: {}".format(N_SNe_so_far),
                    coord_system="axis",
                    inset_box_args={"facecolor": "darkslategray",
                                    "alpha": 0.9},
                    )
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


def plot_phase_diagram(ts, snapshot_number, snapshot_number_to_index_map,
                       SN_times, plots_dir,
                       weight_field,
                       field_type="all",
                       save_plot=True,
                       show_plot=True,
                       seaborn_style="ticks"):
    """ Creates [and optionally saves] a density-temperature phase diagram

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
                               x_bins=100,
                               y_bins=100,)

    if weight_field == "particle_mass":
        ppp.profile.set_field_unit("particle_mass", "Msun")
        ppp.plots[(field_type, "particle_mass")].zmin, ppp.plots[(field_type, "particle_mass")].zmax = (None, None)

    ppp.set_unit("density", "amu/cm**3")
    ppp.set_xlabel("Density $ [ m_\mathrm{H} \; \mathrm{cm}^{-3} ] $")

    ppp.set_xlim(1e-4, 1e3)
    ppp.set_ylim(1e1, 1e7)

    ppp.set_log("density", True)
    ppp.set_log("temperature", True)
    ppp.set_log("particle_mass", True)

    ppp.set_background_color(weight_field, color="white")

    if show_plot:
        with sns.axes_style(seaborn_style):
            ppp.show()

    if save_plot:
        subdir = "phase_plot"
        if not os.path.exists(os.path.join(plots_dir, subdir)):
            os.mkdir(os.path.join(plots_dir, subdir))
        plot_name = os.path.join(subdir, "snapshot_{:0>3}_{}".format(snapshot_number, weight_field))
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
    """ `current_time` should be in Myr """
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
    """ Creates [and optionally saves] a radial profile plot of `field`

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
