#!/usr/bin/env python
"""Module for CpHMD replica exchange on GPUs

Used to run CpHMD replica exchange on multiple GPUs
with the CpHMD GPU Amber code
"""

import multiprocessing
import subprocess
import random
import os
import math
import MDAnalysis as mda
import f90nml
import glob
from shutil import copyfile


class PrintFreqError(Exception):

    """Exception for checking if print frequency less than swap frequency"""

    if __name__ == "__main__":
        print('The print frequency must be <= the swap frequency.')


class NoRestartError(Exception):

    """Exception for checking if the user provided a rst file on an initial run"""

    if __name__ == "__main__":
        print('On an initial run the user must provide a rst7 file.')


def run_replica(md_set, index, swap_attempts, cuda_devices, phs, output_directory,
                template_phmdin_file, template_mdin, initial_restart_file, amberdir, library_path,
                amber_parm_file, phmd_parm_file, restart):

    """Internal function for replica exchange -- DO NOT CALL"""

    job_index = md_set * len(cuda_devices) + index
    swap_count = math.floor(job_index / len(phs))
    if swap_count < swap_attempts:
        rep = job_index - len(phs) * swap_count
        if swap_count == 0:
            local_ph = phs[rep]
            mdin_nml = f90nml.read(template_mdin)
            save_traj_freq = mdin_nml['cntrl']['ntwx']
            md_steps_per_swap = mdin_nml['cntrl']['nstlim']
            if (save_traj_freq >= md_steps_per_swap) and (save_traj_freq % md_steps_per_swap != 0):
                print("The template file {} has improper value for 'ntwx'.".format(template_mdin))
                raise Exception("InputError: bad template file.")
            elif save_traj_freq < md_steps_per_swap:
                mdin_nml['cntrl']['ntwx'] = save_traj_freq
            if save_traj_freq > md_steps_per_swap:
                mdin_nml['cntrl']['ntwx'] = 0
            mdin_nml['cntrl']['solvph'] = float(local_ph)
            mdin_nml.write('{}/rep{}.mdin'.format(output_directory, rep))
            phmdin_nml = f90nml.read(template_phmdin_file)
            phmdin_nml['phmdin']['qphmdstart'] = bool(not restart)
            phmdin_nml.write('{}/input{}.phmdin'.format(output_directory, rep))
            if not restart:
                amber_rst_file = initial_restart_file
            else:
                amber_rst_file = '{}/rep{}.rst7'.format(output_directory, rep)
                phmdrst_string = '{}/rep{}.phmdstrt'.format(output_directory, rep)
            lambda_file = '{}/rep{}.lambda'.format(output_directory, rep)
            nc_file = '{}/rep{}.nc'.format(output_directory, rep)
            save_mdout_freq = mdin_nml['cntrl']['ntpr']
            if save_mdout_freq == 0:
                mdout_file = '/dev/null'
            else:
                mdout_file = '{}/rep{}.mdout'.format(output_directory, rep)
        else:
            mdin_nml = f90nml.read('{}/rep{}.mdin'.format(output_directory, rep))
            tmp_mdin_nml = f90nml.read(template_mdin)
            save_traj_freq = tmp_mdin_nml['cntrl']['ntwx']
            md_steps_per_swap = mdin_nml['cntrl']['nstlim']
            if swap_count == 1:
                phmdin_nml = f90nml.read(template_phmdin_file)
                phmdin_nml['phmdin']['qphmdstart'] = False
                os.remove('{}/input{}.phmdin'.format(output_directory, rep))
                phmdin_nml.write('{}/input{}.phmdin'.format(output_directory, rep))
            if (swap_count + 1) % (save_traj_freq // md_steps_per_swap) == 0:
                mdin_nml['cntrl']['ntwx'] = md_steps_per_swap
                os.remove('{}/rep{}.mdin'.format(output_directory, rep))
                mdin_nml.write('{}/rep{}.mdin'.format(output_directory, rep))
            else:
                mdin_nml['cntrl']['ntwx'] = 0
                os.remove('{}/rep{}.mdin'.format(output_directory, rep))
                mdin_nml.write('{}/rep{}.mdin'.format(output_directory, rep))
            amber_rst_file = '{}/rep{}.rst7'.format(output_directory, rep)
            phmdrst_string = '{}/rep{}.phmdstrt'.format(output_directory, rep)
            lambda_file = '{}/rep{}_current.lambda'.format(output_directory, rep)
            nc_file = '{}/rep{}_current.nc'.format(output_directory, rep)
            save_mdout_freq = mdin_nml['cntrl']['ntpr']
            if save_mdout_freq == 0:
                mdout_file = '/dev/null'
            else:
                mdout_file = '{}/rep{}_current.mdout'.format(output_directory, rep)
        phmdrst_out_file = '{}/rep{}.phmdrestrt'.format(output_directory, rep)
        rst_out_file = '{}/rep{}_current.rst7'.format(output_directory, rep)
        phmdin_file = '{}/input{}.phmdin'.format(output_directory, rep)
        with open('{}/out.log'.format(output_directory), 'a') as log_file, \
             open('{}/out.err'.format(output_directory), 'a') as error_file:
            if swap_count == 0 and not restart:
                subprocess.run(['{}/bin/pmemd.cuda'.format(amberdir), '-O', '-i',
                                 '{}/rep{}.mdin'.format(output_directory, rep), '-p', '{}'.format(amber_parm_file),
                                 '-c', '{}'.format(amber_rst_file), '-ref', '{}'.format(amber_rst_file), '-phmdparm',
                                 '{}'.format(phmd_parm_file), '-phmdin', '{}'.format(phmdin_file),
                                 '-phmdout', '{}'.format(lambda_file), '-o',
                                 '{}'.format(mdout_file), '-r', '{}'.format(rst_out_file), '-x',
                                 '{}'.format(nc_file), '-phmdrestrt',
                                 '{}'.format(phmdrst_out_file)],
                                stdout=log_file, stderr=error_file,
                                env={'CUDA_VISIBLE_DEVICES': cuda_devices[index],
                                     'LD_LIBRARY_PATH': library_path})
            else:
                subprocess.run(['{}/bin/pmemd.cuda'.format(amberdir), '-O', '-i',
                                 '{}/rep{}.mdin'.format(output_directory, rep), '-p', '{}'.format(amber_parm_file),
                                 '-c', '{}'.format(amber_rst_file), '-ref', '{}'.format(amber_rst_file), '-phmdparm',
                                 '{}'.format(phmd_parm_file), '-phmdin', '{}'.format(phmdin_file),
                                 '-phmdstrt', '{}'.format(phmdrst_string), '-phmdout',
                                 '{}'.format(lambda_file), '-o', '{}'.format(mdout_file), '-r',
                                 '{}'.format(rst_out_file), '-x', '{}'.format(nc_file),
                                 '-phmdrestrt', '{}'.format(phmdrst_out_file)],
                                stdout=log_file, stderr=error_file,
                                env={'CUDA_VISIBLE_DEVICES': cuda_devices[index],
                                     'LD_LIBRARY_PATH': library_path})


def append_files(md_set, cuda_devices, phs, replica_steps, output_directory, concatenate_parm_file, template_mdin, swap_attempts):

    """Internal function for replica exchange -- DO NOT CALL"""

    for i in range(len(cuda_devices)):
        job_index = md_set * len(cuda_devices) + i
        if job_index == len(phs) * swap_attempts:
            break
        swap_count = math.floor(job_index / len(phs))
        ph_index = job_index - swap_count * len(phs)
        if swap_count > 0:
            os.remove('{}/rep{}.rst7'.format(output_directory, ph_index))
        if os.path.exists('{}/rep{}_current.rst7'.format(output_directory, ph_index)):
            os.rename('{}/rep{}_current.rst7'.format(output_directory, ph_index),
                    '{}/rep{}.rst7'.format(output_directory, ph_index))
        with open('{}/rep{}.phmdrestrt'.format(output_directory, ph_index), 'r') as infile:
            with open('{}/rep{}.phmdstrt'.format(output_directory, ph_index), 'w') as outfile:
                for line in infile:
                    outfile.write(line.replace('PHMDRST', 'PHMDSTRT'))
        if swap_count != 0:
            mdin_nml = f90nml.read(template_mdin)
            save_traj_freq = mdin_nml['cntrl']['ntwx']
            md_steps_per_swap = mdin_nml['cntrl']['nstlim']
            with open('{}/rep{}_current.lambda'.format(output_directory, ph_index), 'r') as in_file:
                with open('{}/rep{}.lambda'.format(output_directory, ph_index), 'a') as outfile:
                    for line in in_file:
                        if line[0] != '#':
                            steps = int(line.split()[0])
                            steps += replica_steps * swap_count
                            print_line = '{:8}'.format(steps) + line[8:]
                            outfile.write(print_line)
            if (swap_count + 1) % (save_traj_freq // md_steps_per_swap) == 0:
                if (swap_count + 1) // (save_traj_freq // md_steps_per_swap) == 1:
                    os.rename('{}/rep{}_current.nc'.format(output_directory, ph_index),
                              '{}/rep{}.nc'.format(output_directory, ph_index))
                else:
                    universe = mda.Universe(concatenate_parm_file, '{}/rep{}.nc'.format(output_directory, ph_index),
                                            format='NCDF')
                    all_selection = universe.select_atoms('all')
                    with mda.Writer('{}/temp{}_{}.ncdf'.format(output_directory, ph_index, swap_count),
                                    all_selection.n_atoms) as writer:
                        for dummy in universe.trajectory:
                            writer.write(all_selection)
                        universe = mda.Universe(concatenate_parm_file,
                                                '{}/rep{}_current.nc'.format(output_directory, ph_index), format='NCDF')
                        all_selection = universe.select_atoms('all')
                        for dummy in universe.trajectory:
                            writer.write(all_selection)
                    os.remove('{}/rep{}.nc'.format(output_directory, ph_index))
                    os.rename('{}/temp{}_{}.ncdf'.format(output_directory, ph_index, swap_count),
                            '{}/rep{}.nc'.format(output_directory, ph_index))

            save_mdout_freq = mdin_nml['cntrl']['ntpr']
            if save_mdout_freq != 0:        
                with open('{}/temp{}_{}.mdout'.format(output_directory, ph_index, swap_count), 'w') as outfile:
                    with open('{}/rep{}.mdout'.format(output_directory, ph_index), 'r') as mdout:
                        for line in mdout:
                            if line.find('A V E R A G E S') != -1:
                                break
                            if len(line) > 4 and line[0:4] == ' WTH':
                                continue
                            outfile.write(line)
                    printing = False
                    with open('{}/rep{}_current.mdout'.format(output_directory, ph_index), 'r') as mdout:
                        for line in mdout:
                            if line.find('RESULTS') != -1:
                                printing = True
                                line = mdout.readline()
                                line = mdout.readline()
                                line = mdout.readline()
                            if line.find('A V E R A G E S') != -1:
                                break
                            if printing:
                                outfile.write(line)
                os.remove('{}/rep{}.mdout'.format(output_directory, ph_index))
                os.rename('{}/temp{}_{}.mdout'.format(output_directory, ph_index, swap_count),
                        '{}/rep{}.mdout'.format(output_directory, ph_index))


def do_swaps(md_set, phs, cuda_devices, log_file, total_swaps, replicas, output_directory):

    """Internal function for replica exchange -- DO NOT CALL"""

    for i in range(len(cuda_devices)):
        index = md_set * len(cuda_devices) + i
        swap_index = math.floor(index / len(phs))
        ph_index = index - swap_index * len(phs)
        local_ph = phs[ph_index]
        if swap_index % 2 != ph_index % 2 and ph_index != 0 and swap_index < total_swaps:
            other_ph_index = ph_index - 1
            other_ph = phs[other_ph_index]
            if swap_index == 0:
                lambda_file = '{}/rep{}.lambda'.format(output_directory, ph_index)
            else:
                lambda_file = '{}/rep{}_current.lambda'.format(output_directory, ph_index)
            line = ''
            with open(lambda_file, 'r') as in_file:
                for line in in_file:
                    if line.find('itauto') != -1:
                        tautos = [int(tauto) for tauto in line.split()[2:]]
            proton_count = (sum([1 - float(mylambda) for j, mylambda in
                            enumerate(line.split()[1:]) if (tautos[j] == 0 or
                                                            tautos[j] == 1 or
                                                            tautos[j] == 3)]))
            if swap_index == 0:
                other_lambda_file = '{}/rep{}.lambda'.format(output_directory, other_ph_index)
            else:
                other_lambda_file = '{}/rep{}_current.lambda'.format(output_directory, other_ph_index)
            line = ''
            with open(other_lambda_file, 'r') as in_file:
                for line in in_file:
                    pass
            other_proton_count = (sum([1 - float(mylambda) for j, mylambda in
                                  enumerate(line.split()[1:]) if (tautos[j] == 0 or
                                                                  tautos[j] == 1 or
                                                                  tautos[j] == 3)]))
            delta = math.log(10.0) * ((proton_count - other_proton_count) *
                                      (float(other_ph) - float(local_ph)))
            if delta < 0 or math.exp(-delta) > random.random():
                os.rename('{}/rep{}.phmdstrt'.format(output_directory, ph_index),
                          '{}/temp.phmdstrt'.format(output_directory))
                os.rename('{}/rep{}.phmdstrt'.format(output_directory, other_ph_index),
                          '{}/rep{}.phmdstrt'.format(output_directory, ph_index))
                os.rename('{}/temp.phmdstrt'.format(output_directory),
                          '{}/rep{}.phmdstrt'.format(output_directory, other_ph_index))
                os.rename('{}/rep{}.rst7'.format(output_directory, ph_index), '{}/temp.rst7'.format(output_directory))
                os.rename('{}/rep{}.rst7'.format(output_directory, other_ph_index),
                          '{}/rep{}.rst7'.format(output_directory, ph_index))
                os.rename('{}/temp.rst7'.format(output_directory),
                          '{}/rep{}.rst7'.format(output_directory, other_ph_index))
                temp = replicas[other_ph_index]
                replicas[other_ph_index] = replicas[ph_index]
                replicas[ph_index] = temp
        if ph_index == len(phs) - 1 and swap_index < total_swaps - 1:
            with open(log_file, 'a') as out_file:
                out_file.write('swap {}'.format(swap_index + 1))
                for replica in replicas:
                    out_file.write(' {}'.format(replica))
                out_file.write('\n')


def prepare_restart(phs, output_directory, restart_directory):

    """Function for restarting replica exchange -- DO NOT CALL"""

    mdin_files = glob.glob('{}/rep*.mdin'.format(restart_directory))
    old_number_of_replicas = max([int(mdin_file[(mdin_file.find('rep') + 3):mdin_file.find('.mdin')])
                                  for mdin_file in mdin_files]) + 1

    old_phs = []
    for i in range(old_number_of_replicas):
        mdin_file = f90nml.read('{}/rep{}.mdin'.format(restart_directory, i))
        old_phs.append(mdin_file['cntrl']['solvph'])

    for i, ph in enumerate(phs):
        ph_differences = [abs(float(ph) - old_ph) for old_ph in old_phs]
        rep = ph_differences.index(min(ph_differences))
        copyfile('{}/rep{}.rst7'.format(restart_directory, rep),
                 '{}/rep{}.rst7'.format(output_directory, i))
        copyfile('{}/rep{}.phmdstrt'.format(restart_directory, rep),
                 '{}/rep{}.phmdstrt'.format(output_directory, i))


def do_replica_exchange(phs, swap_attempts, template_mdin, amber_parm_file,
                        phmd_parm_file, template_phmdin_file,
                        output_directory, restart_directory=None, concatenate_parm_file=None,
                        restart=False, initial_restart_file=None):

    """Main function for coordinating replica exchange

    Keyword arguments:
    phs -- list of pH values for replicas, must be in order
    swap_attempts -- how many swaps to attempt
    template_mdin -- template Amber mdin file specifying how
                     to run a single replica DO NOT USE NUMEXCHG
                     ALWAYS USE IREST = 1, NTX = 5
    amber_parm_file -- Amber parm7 file
    phmd_parm_file -- CpHMD Amber parameter file
    template_phmdin_file -- CpHMD phmdin control file
    initial_restart_file -- initial amber rst7 file
                            MUST CONTAIN VELOCITIES
    output_directory -- new directory to store output of calculation
                        ensures files not overwritten
    concatenate_parm_file -- if neccessary, parm file for concatenation
                             of trajectories with MDAnalysis. Eg.,
                             if running with CHARMM force field,
                             this should be psf file passed to
                             Chamber
    restart_directory -- directory containing previous run
                         only used if restart = True
    restart -- set to True if you are restarting a replica-exchange run
    """

    cuda_device_string = os.environ['CUDA_VISIBLE_DEVICES']

    amberdir = os.environ['AMBERHOME']

    try:
        library_path = os.environ['LD_LIBRARY_PATH']
    except KeyError:
        pass

    if concatenate_parm_file is None:
        concatenate_parm_file = amber_parm_file

    replica_log = '{}/reps.log'.format(output_directory)
    os.mkdir(output_directory)
    if restart:
        prepare_restart(phs, output_directory, restart_directory)

    mdin_nml = f90nml.read(template_mdin)
    replica_steps = mdin_nml['cntrl']['nstlim']

    # print_steps = mdin_nml['cntrl']['ntwx']
    # # if print_steps > replica_steps:
    # #     raise PrintFreqError

    if initial_restart_file is None and not restart:
        raise NoRestartError

    replicas = [i for i, ph in enumerate(phs)]

    cuda_devices = cuda_device_string.split(',')
    total_md_runs = len(phs) * swap_attempts
    md_sets = math.ceil(total_md_runs / len(cuda_devices))
    with open(replica_log, 'w') as replica_file:
        replica_file.write('swap 0')
        for replica in replicas:
            replica_file.write(' {}'.format(replica))
        replica_file.write('\n')
    for md_set in range(md_sets):
        with multiprocessing.Pool(processes=len(cuda_devices)) as pool:
            arguments = [[md_set, i, swap_attempts, cuda_devices, phs, output_directory, template_phmdin_file,
                          template_mdin, initial_restart_file, amberdir, library_path,
                          amber_parm_file, phmd_parm_file, restart]
                         for i, device in enumerate(cuda_devices)]
            pool.starmap(run_replica, arguments)
        if swap_attempts > 1:
            append_files(md_set, cuda_devices, phs, replica_steps, output_directory, concatenate_parm_file, template_mdin, swap_attempts)
            do_swaps(md_set, phs, cuda_devices, replica_log, swap_attempts, replicas, output_directory)
        else:
            append_files(md_set, cuda_devices, phs, replica_steps, output_directory, concatenate_parm_file, template_mdin, swap_attempts)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description="""A Python program to run asynchronous CpHMD replica exchange simulations.""")
    parser.add_argument(
        '-p', '--pHs', default='7.0 7.5 8.0 8.5', dest='pHs', metavar='PHS',
        help="""A list of pH conditions to run CpHMD replica exchange simulations.
            [default: '7.0 7.5 8.0 8.5']""")
    parser.add_argument(
        '-s', '--swap_attempts', default=2000, dest='attemps', metavar='ATTEMPS',
        help="""Number of swaps to attempt during replica exchange.
            [default: 2000]""")
    parser.add_argument(
        '-tm', '--template_mdin', default='template_mdin', dest='mdin', metavar='MDIN',
        help="""Template Amber mdin file specifying how to run a single replica DO NOT USE NUMEXCHG
            ALWAYS USE IREST = 1, NTX = 5.
            [default: 'template_mdin']""")
    parser.add_argument(
        '-tp', '--template_phmdin_file', default='template_phmdin_file', dest='phmdin', metavar='PHMDIN',
        help="""Template Amber mdin file specifying how to run a single replica DO NOT USE NUMEXCHG
            ALWAYS USE IREST = 1, NTX = 5.
            [default: 'template_phmdin']""")
    parser.add_argument(
        '-ap', '--amber_parm_file', default='protein.parm7', dest='parm7', metavar='PAMM7',
        help="""Amber parm7 file for the system.
            [default: 'protein.parm7']""")
    parser.add_argument(
        '-pp', '--phmd_parm_file', default='input.parm', dest='phmdparm', metavar='PHMDPARM',
        help="""CpHMD Amber parameter file.
            [default: 'input.parm']""")
    parser.add_argument(
        '-rf', '--initial_restart_file', default='equil.rst7', dest='rst7', metavar='RST7',
        help="""Initial Amber rst7 file. MUST CONTAIN VELOCITIES.
            [default: 'equil.rst7']""")
    parser.add_argument(
        '-od', '--output_directory', default='.', dest='dir', metavar='DIR',
        help="""New directory to store output of calculation ensures files not overwritten.
            [default: '.']""")
    parser.add_argument(
        '-cp', '--concatenate_parm_file', default=None, dest='concparm', metavar='CONCPARM',
        help="""If neccessary, parm file for concatenation of trajectories with MDAnalysis. Eg.,
            if running with CHARMM force field, this should be psf file passed to Chamber.
            [default: None]""")
    parser.add_argument(
        '-rd', '--restart_directory', default='prev_restart', dest='rdir', metavar='RDIR',
        help="""Directory containing previous run only used if restart = True.
            [default: 'prev_restart']""")
    parser.add_argument(
        '-rst', '--restart', default=False, dest='restart', metavar='RESTAET',
        help="""Set to True if you are restarting a replica-exchange run
            [default: 'prev_restart']""")
    args = parser.parse_args()
    do_replica_exchange(args.pHs.split(' '), args.attemps, args.mdin, args.parm7, args.phmdparm, args.phmdin,
                        args.dir, restart_directory=args.rdir, concatenate_parm_file=args.concparm,
                        restart=args.restart, initial_restart_file=args.rst7)
