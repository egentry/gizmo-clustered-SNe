from __future__ import print_function, division

import os
import numpy as np
import pandas as pd

from scipy import interpolate, integrate

from matplotlib import pyplot as plt

# ######### FUNCTIONS


def get_field_lines_filename_from_ds(ds, plots_dir):
    filename = os.path.join(plots_dir,
                            "field_lines",
                            "field_lines_{}.csv".format(os.path.splitext(ds.basename)[0]))

    return filename


def visualize_magnetic_field_layout(interpolator, box_size=200):
    plt.figure()
    y = np.linspace(-box_size, box_size, num=400)
    z = np.linspace(-box_size, box_size, num=400)
    x = np.array([0])

    xx, yy, zz = np.meshgrid(x, y, z)

    interpolated = interpolator(xx, yy, zz).squeeze()

    By_Bz = interpolated[:, :, 1] / interpolated[:, :, 2]

    min_value = -10
    max_value = +10

    By_Bz[By_Bz < min_value] = np.nan
    By_Bz[By_Bz > max_value] = np.nan

    print(interpolated.shape)

    cmap = plt.get_cmap("bwr")
    cmap.set_bad(color="black", alpha=1.)

    plt.pcolormesh(yy.squeeze(), zz.squeeze(),
                   By_Bz,
                   cmap=cmap,
                   vmin=-2, vmax=+2,
                   )

    plt.colorbar(label="$B_y$ / $B_z$")
    plt.xlabel("$y$ [pc]")
    plt.ylabel("$z$ [pc]")
    plt.show()


def mag_derivs(y, z, interpolator):
    # to be integrated using scipy.integrate.odeint
    pos = np.array([0, y, z])
    B = interpolator(*pos)
    deriv = B[1] / B[2]
    deriv = np.clip(deriv, -15, 15)
    return deriv


def plot_field_lines(df_lines, ax,
                     plot_kwargs={"color": "white"},
                     unit_scaling=1,
                     replace_nans=False):
    for line_id in df_lines.line_id.drop_duplicates():
        df_tmp = df_lines.set_index("line_id").loc[line_id].copy()

        y_final = df_tmp.dropna().sort_values("z").iloc[-1].y
        if replace_nans:
            df_tmp.loc[df_tmp.y.isnull(), "y"] = y_final

        ax.plot(df_tmp.y/unit_scaling, df_tmp.z/unit_scaling,
                **plot_kwargs)


def calculate_magnetic_field_lines(ds, plots_dir, make_plots=True):
    dd = ds.all_data()
    coords = (dd["PartType0", "Coordinates"] - ds.domain_center).value
    magnetic_field = dd["PartType0", "MagneticField"].value

    box_size = ds.domain_width.value[0]

    slice_filter = np.abs(coords[:, 0]) < 40
    interpolator = interpolate.LinearNDInterpolator(
        coords[slice_filter],
        magnetic_field[slice_filter],
    )

    if make_plots:
        visualize_magnetic_field_layout(interpolator, box_size=box_size/4)

    # # INTEGRATE
    y_0s = np.linspace(-box_size/2, box_size/2, num=16)[1:-1]

    max_coord = box_size/2

    ys_integrated = [None]*y_0s.size

    for i, y_0 in enumerate(y_0s):

        zs = np.arange(-max_coord, max_coord, step=.25)[1:-1]
        y_integrated = integrate.odeint(mag_derivs, y_0, zs,
                                        args=(interpolator,))
        ys_integrated[i] = y_integrated

    # # SAVE RESULTS
    df_lines = pd.concat(
        [pd.DataFrame(
            data={"line_id": i, "y": y_integrated.squeeze(), "z": zs})
         for i, y_integrated in enumerate(ys_integrated)])

    filename = get_field_lines_filename_from_ds(ds, plots_dir)

    os.makedirs(os.path.dirname(filename), exist_ok=True)
    df_lines.to_csv(filename, index=False)

    if make_plots:
        plot_field_lines(df_lines, plt.axes(aspect="equal"),
                         plot_kwargs={"color": "black"})

    return df_lines
