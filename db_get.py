#! /usr/bin/env python
from __future__ import print_function
import os
from argparse import ArgumentParser
from matgendb import QueryEngine
import pymatgen as pmg
from pymatgen.io.vasp.outputs import VaspParserError
from pymatgen.electronic_structure.plotter import BSPlotter


def getSummary(t_id, query_engine=QueryEngine()):
    from prettytable import PrettyTable
    c = query_engine.query_one(criteria={"task_id": t_id},
                               properties=['output.final_energy',
                                           'output.crystal',
                                           'calculations.output'])
    struct = pmg.Structure.from_dict(c['output.crystal'])
    # is something always stored for f? may need to adjust if not
    mag_moms = PrettyTable(field_names=['Element',
                                        's', 'p', 'd', 'f',
                                        'total'])
    mags = c['calculations.output'][0]['outcar']['magnetization']
    for i, m in enumerate(mags):
        mag_moms.add_row([struct[i].species_string,
                          m['s'], m['p'], m['p'], m['d'],
                          m['tot']])
    run_stats = PrettyTable(header=False)
    run_stats_dict = c['calculations.output'][0]['outcar']['run_stats']
    for k, v in sorted(run_stats_dict.iteritems()):
        run_stats.add_row([k, v])
    force_table = PrettyTable(field_names=['Element',
                                           'position', 'x', 'y', 'z'])
    forces = c['calculations.output'][0]['ionic_steps'][-1]['forces']
    for i, f in enumerate(forces):
        force_table.add_row([struct[i].species_string]
                            + [str(struct[i].coords)]
                            + f)
    summary = (str(struct)
               + '\nTotal Energy: {} eV\n'.format(c['output.final_energy'])
               + 'Band gap: {} eV\n'.format(
                   c['calculations.output'][0]['bandgap'])
               + 'Magnetic moments:\n {} \n'.format(mag_moms)
               + 'Total Magnetization: {}\n'.format(
                   c['calculations.output'][0]['outcar']['total_magnetization'])
               + 'Forces:\n{}\n'.format(force_table)
               + 'Run Statistics:\n'.format(run_stats)
               )
    return summary


def getInputs(t_id, query_engine=QueryEngine()):
    """
    Args:
        t_id:
            task_id of the run to get input objects for
    Returns:
        A tuple of (Incar, Kpoints, Poscar, PspInfo)
            all but PspInfo are pymatgen objects
        PspInfo is a list of dicts with keys 'hash' and 'titel'(pymatgen typo?)
            hash is a hash of the POTCAR file, titel is a string of functional,
            POTCAR symbol, and date it was made
    """
    c = query_engine.query_one(criteria={"task_id": t_id},
                               properties=['calculations.input',
                                           'input.potcar_spec'])
    inp = c['calculations.input'][0]
    incar = inp['incar']
    Incar = pmg.io.vasp.Incar.from_dict(incar)
    kpts = inp['kpoints']
    # pymatgen does not recognize Monkhorst-Pack as a generation_stlye
    # there may be other cases where vasp is more flexible
    if kpts['generation_style'] == 'Monkhorst-Pack':
        kpts['generation_style'] = 'Monkhorst'
    Kpoints = pmg.io.vasp.Kpoints.from_dict(kpts)
    crystal = inp['crystal']
    struct = pmg.Structure.from_dict(crystal)
    Poscar = pmg.io.vasp.Poscar(struct)
    PspInfo = c['input.potcar_spec']
    return (Incar, Kpoints, Poscar, PspInfo)


def plotBands(t_id, filename=" ", query_engine=QueryEngine(),
              band_range=None, kpoints_file=None):
    # TODO: make this work for ISPIN=1 calculations
    res = query_engine.query_one(criteria={"task_id": t_id},
                                 properties=['calculations.output.eigenvalues',
                                             'calculations.output.crystal.lattice',
                                             'calculations.output.efermi',
                                             'calculations.input.kpoints'])
    eigs = res['calculations.output.eigenvalues'][0]
    lat = pmg.Lattice.from_dict(res['calculations.output.crystal.lattice'][0])
    efermi = res['calculations.output.efermi'][0]

    # pymatgen object that knows about high-sym pts
    Kpoints = pmg.io.vasp.Kpoints.from_dict(
        res['calculations.input.kpoints'][0])
    if not Kpoints.style.value == 3:
        raise VaspParserError(
            'Can only plot bands for a run which was preformed with Line_mode kpoints')
    if kpoints_file:
        Kpoints = pmg.io.vasp.Kpoints.from_file(kpoints_file)
    # raw data of actual kpoints
    kpts = res['calculations.input.kpoints'][0]['actual_points']
    # list of actual kpoints
    kpoints = [kpt['abc'] for kpt in kpts]

    # reorganize eigenvalues because the database way of storing them is no good for plotting
    kpt_is = [str(key) for key in sorted([int(k) for k in eigs.keys()])]
    eigs_up = [[eigs[kpt_i]['1'][band_i][0] for kpt_i in kpt_is]
               for band_i in range(len(eigs['0']['1']))]
    eigs_dn = [[eigs[kpt_i]['-1'][band_i][0] for kpt_i in kpt_is]
               for band_i in range(len(eigs['0']['-1']))]
    eigenvals = {pmg.Spin.up: eigs_up, pmg.Spin.down: eigs_dn}

    have_sym_labels = False
    if len(set(Kpoints.labels)) == 1 and next(iter(set(Kpoints.labels))) == 'unknown':
        # want set(Kpoints.kpts) but can't get set of lists,
        # there is faster solution with itertools
        # but this list is small so I go for simplicity
        unique_kpts = []
        for k in Kpoints.kpts:
            if k not in unique_kpts:
                unique_kpts.append(k)
        kpt_labels = {str(v): v for i, v in enumerate(unique_kpts)}
    else:
        kpt_labels = {k: v for k, v in zip(Kpoints.labels, Kpoints.kpts)}
        have_sym_labels = True

    bs = pmg.io.vasp.BandStructureSymmLine(kpoints, eigenvals, lat,
                                           efermi, labels_dict=kpt_labels)
    bsp = BSPlotter(bs)
    bs_dict = bsp.bs_plot_data()
    fig = plt.figure()
    bands_sp = fig.add_subplot(111)
    n_bands = len(bs_dict['energy'][0]['1'])
    for d in range(len(bs_dict['distances'])):
        for i in range(n_bands):
            bands_sp.plot(bs_dict['distances'][d], bs_dict['energy'][d]['1'][i], 'b-')
            bands_sp.plot(bs_dict['distances'][d], bs_dict['energy'][d]['-1'][i], 'r--')
    x_max = bs_dict['distances'][-1][-1]
    plt.xlim(0, x_max)
    if band_range:
        plt.ylim(band_range[0], band_range[1])
    bands_sp.tick_params(labelsize=25)
    bands_sp.set_xticks(bs_dict['ticks']['distance'])
    if have_sym_labels:
        bands_sp.set_xticklabels(bs_dict['ticks']['label'], size=30)
    else:
        bands_sp.set_xticklabels(bs_dict['ticks']['label'], rotation=90, size=12)
    for pt in bs_dict['ticks']['distance']:
        bands_sp.axvline(pt, color='k')
    # bands_sp.axhline(0., color='k', linewidth=2)
    bands_sp.set_ylabel('$E-E_F$ (eV)', size=30)
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    arg_parser = ArgumentParser()
    arg_parser.add_argument("t_id", metavar="t_id", type=int,
                            help="task_id of the database entry to gather data from")
    arg_parser.add_argument("-s", "--summary", action="store_true", required=False,
                            help="print a summary of the calculation results")
    arg_parser.add_argument("-i", "--inputs", nargs="?", type=str, const=" ", required=False,
                            help=("recreate the VASP input files used for this calculation "
                                  "in the directory provided, if no directory is provided "
                                  "all files will be printed"))
    arg_parser.add_argument("-b", "--bands", nargs="?", type=str, const=" ", required=False,
                            help=("plot the band structure "
                                  "for a run with k-points along lines, "
                                  "optionally give a filename to save the plot"))
    arg_parser.add_argument("-br", "--band_range", nargs=2, type=float, required=False,
                            help=("Range of energies over which to plot band structure "
                                  "in eV give two numbers separated by a space, "
                                  "default is entire eigenvalue range"))
    arg_parser.add_argument("-bl", "--band_labels", type=str, required=False,
                            help=("Provide a KPOINTS file with high symmetry points labeled "
                                  "so these labels can be used in the band structure plot"))
    arg_parser.add_argument("-si", "--structure_initial", nargs="?",
                            type=str, const=" ", required=False,
                            help=("With no arguments will print the initial structure, "
                                  "optionally give a file name to create a file containing "
                                  "the starting structure of the calculation, "
                                  "filetype is based on the file name given here, "
                                  "i.e. 'BTO.cif' will produce a cif file, "
                                  "'POSCAR' will produce a vasp input file"))
    arg_parser.add_argument("-sf", "--structure_final", nargs="?",
                            type=str, const=" ", required=False,
                            help=("With no arguments will print the final structure, "
                                  "optionally give a file name to create a file containing "
                                  "the final structure of the calculation, "
                                  "filetype is based on the file name given here, "
                                  "i.e. 'BTO.cif' will produce a cif file, "
                                  "'POSCAR' will produce a vasp input file"))
    arg_parser.add_argument("-c", "--config_file", required=False,
                            default="/home/public/perovskite2016/DB/db.json",
                            help="location of db.json config file")
    args = arg_parser.parse_args()

    qe = QueryEngine.from_config(args.config_file)
    if args.summary:
        print(getSummary(args.t_id, query_engine=qe))
    if args.inputs:
        input_data = getInputs(args.t_id, query_engine=qe)
        if args.inputs.strip(" "):
            if not os.path.exists(args.inputs):
                os.makedirs(args.inputs)
            input_data[0].write_file(args.inputs+'/INCAR')
            input_data[1].write_file(args.inputs+'/KPOINTS')
            input_data[2].write_file(args.inputs+'/POSCAR')
            with open(args.inputs+'/pspInfo', 'w') as f:
                for psp in input_data[3]:
                    f.write(psp['titel']+'\t'+psp['hash']+'\n')
        else:
            print('INCAR:')
            print(str(input_data[0]), "\n")
            print('KPOINTS:')
            print(str(input_data[1]), "\n")
            print('POSCAR:')
            print(str(input_data[2]), "\n")
            print('pspInfo:')
            for psp in input_data[3]:
                print(psp['titel']+'\t'+psp['hash'])
    if args.bands:
        from matplotlib import pyplot as plt
        fig = plotBands(args.t_id, filename=args.bands, query_engine=qe,
                        band_range=args.band_range, kpoints_file=args.band_labels)
        if args.bands.strip(" "):
            plt.savefig(args.bands)
        else:
            plt.show()
    if args.structure_initial:
        struct = qe.get_structure_from_id(args.t_id, final_structure=False)
        if args.structure_initial.strip(" "):
            struct.to(filename=args.structure_initial)
        else:
            print(struct)
    if args.structure_final:
        struct = qe.get_structure_from_id(args.t_id, final_structure=True)
        if args.structure_final.strip(" "):
            struct.to(filename=args.structure_final)
        else:
            print(struct)
