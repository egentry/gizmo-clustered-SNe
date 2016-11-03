#!/usr/bin/python

import os, sys, glob, h5py
import numpy as np



from injection_helpers import create_snapshot_with_new_SN, \
    create_restart_params_wrapped, \
    get_ith_snapshot_file_in_dir, \
    get_SNe, \
    create_restart_type_file, \
    Params, \
    RestartType


from units import M_solar, m_proton, pc, yr, Myr, km, s, gamma

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("end_file", help="location of 'end' sentinel file, relative to scripts dir")
parser.add_argument("inputs_dir", help="dir of ICs and params, relative to scripts dir")
parser.add_argument("--from-checkpoint", 
    help="Force a restart from checkpoint files, even if we could have tried using restartfiles",
    action="store_true")
args = parser.parse_args()


if __name__ == "__main__":

    # roughly: check if last_time =? TimeMax
    # - if last_time = TimeMax, inject SN, create restart params
    # - if last time > TimeMax, throw error
    # - if last time > TimeMax, do nothing
    #   - create bash file with RESTART_FLAG
    restart_params_filename = glob.glob(os.path.join(args.inputs_dir, "*params.restart"))
    assert(len(restart_params_filename)==1)
    restart_params_filename = restart_params_filename[0]
    current_params = Params.from_filename(restart_params_filename)
    t_max = current_params.TimeMax
    t_min = current_params.TimeBegin


    output_dir = args.inputs_dir.replace("ICs", "output")
    last_snapshot_filename = get_ith_snapshot_file_in_dir(output_dir, -1)
    with h5py.File(last_snapshot_filename, mode="r") as f:
        t_current = f["Header"].attrs["Time"]
        
    SNe = get_SNe(args.inputs_dir)
    SN_times = np.array([SN.time for SN in SNe])

    t_tol = .5e-3

    if t_current >= 40:
        # don't restart
        print("Simulation already reached nominal end. Not restarting.")


    elif np.isclose(t_current, 0, atol=t_tol, rtol=0):
        raise NotImplementedError("Can't yet handle initially starting simulations")


    elif t_current + t_tol < t_max:
        # probably stopped during normal operation
        print("Restarting without changing params or adding snapshots")

        if t_current < t_min:
            raise RuntimeError("Can't restart as is, if t_current<t_min (t_min={:.3}, t_current={:.3}".format(t_min, t_current))

        if np.isclose(t_current, t_min, rtol=0, atol=t_tol):
            restart_type = RestartType.SNAPSHOTS
        else:
            restart_type = RestartType.RESTARTFILE

        if args.from_checkpoint:
            restart_type = RestartType.SNAPSHOTS

            current_params.InitCondFile = os.path.join(
                os.path.dirname(current_params.InitCondFile),
                os.path.splitext(os.path.basename(last_snapshot_filename))[0]
                )
            current_params.TimeBegin = t_current
            current_params.TimeOfFirstSnapshot = current_params.TimeBegin + current_params.TimeBetSnapshot
            with open(restart_params_filename, mode="w") as params_file:
                current_params.to_file(params_file)


        create_restart_type_file(args.inputs_dir, restart_type)
        os.remove(args.end_file)


    elif t_max + t_tol < t_current:
        print("Creating new params file (but no new snapshots), then restarting")

        # probably has created a snapshot with a new SNe, but hasn't created params file
        # but let's check that

        if np.sum((t_current - 1e-3 >= SN_times) & (t_current > SN_times)) != 1:
            raise RuntimeError("Can't find a single SN about to occur")

        create_restart_params_wrapped(args.inputs_dir)

        create_restart_type_file(args.inputs_dir, RestartType.SNAPSHOTS)
        os.remove(args.end_file)


    elif np.isclose(t_current, t_max, rtol=0, atol=t_tol):
        print("Creating new snapshot, then new params file, then restarting")

        # probably needs a new snapshot (with SN added)
        # but let's check that first

        if np.sum((t_current + 1e-3 >= SN_times) & (t_current < SN_times)) != 1:
            raise RuntimeError("Can't find a single SN about to occur")

        create_snapshot_with_new_SN(args.inputs_dir)
        create_restart_params_wrapped(args.inputs_dir)

        create_restart_type_file(args.inputs_dir, RestartType.SNAPSHOTS)
        os.remove(args.end_file)


    else:
        raise RuntimeError("Timing confusing; take a closer look at this") 

