import argparse
from utils import validate_input, get_structure
from molecular_dynamics import xtb_md
from optimization import xtb_optim


def main(args: argparse.Namespace) -> None:
    if not validate_input(args=args):
        return

    structure, file_extension, output_name = get_structure(args.structure, args.pbc, args.cell)
    n_atoms = len(structure)

    # Extract cell dimensions and create a string from them
    cell_dimensions = ''.join(map(str, map(int, args.cell)))  # Converts [12, 12, 12] to '121212'

    if args.task == 'md':

        md_output_name = (
            f"{output_name}_{args.method}_T_{args.T:.0f}_st_{args.step:.1f}_ac_{args.accu:.1e}_fr_{args.fric:.1f}_mdamp_{args.mdamp:.1f}")
        run_time = xtb_md(args.method, structure, file_extension, md_output_name, T=args.T, N_steps=args.N_steps,
                          accu=args.accu, fric=args.fric, step=args.step, mdamp=args.mdamp)

        print(
            f"""\t------------------------------------------------
            |              MD run completed!               |
            ------------------------------------------------""")
        print(f'\t* MD Time: {run_time:.2f} s, MD Time: {run_time / 60:.2f} min')
        print(f'\t* MD Steps: {args.N_steps}')
        print(f'\t* N Atoms: {n_atoms}')

    elif args.task == 'opt':

        opt_output_name = (f"{output_name}_{args.method}_{args.optimizer}_fm_{args.fmax:.2f}"
                           f"_mst_{args.maxstep:.2f}_ac_{args.scf_accuracy:.1e}_mdamp_{args.mixer_damping:.1f}")

        energy, opt_time, n_steps = xtb_optim(method=args.method,
                                              fmax=args.fmax,
                                              steps=args.steps,
                                              maxstep=args.maxstep,
                                              structure=structure,
                                              file_extension=file_extension,
                                              output_name=opt_output_name,
                                              mixer_damping=args.mixer_damping,
                                              scf_accuracy=args.scf_accuracy,
                                              max_iterations=args.max_iterations,
                                              verbosity=args.verbosity,
                                              charge=args.charge,
                                              multiplicity=args.multiplicity,
                                              optimizer=args.optimizer)

        print(
            f"""\t------------------------------------------------
        |           Optimization completed!            |
        ------------------------------------------------""")
        print(f'\t* Energy: {energy:.3f} kcal/mol')
        print(f'\t* Optim Time: {opt_time:.2f} s, Optim Time: {opt_time / 60:.2f} min')
        print(f'\t* Optim Steps: {n_steps}')
        print(f'\t* N Atoms: {n_atoms}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", choices=['GFN1-xTB', 'GFN2-xTB'], default='GFN1-xTB',
                        help="Select xTB method for structure optimization.")
    parser.add_argument("--task", choices=['md', 'opt'], default='opt',
                        help="Select task (Molecular Dynamics or Optimization).")
    parser.add_argument("--structure", default=None, help="Provide a structure for optimization.")
    parser.add_argument("--pbc", action='store_true', help="Select if structure is periodic.")
    parser.add_argument("--cell", nargs='+', type=float, default=None, help="Provide periodic cell dimensions.")
    parser.add_argument("--T", type=float, default=300.0, help="Provide temperature in K for MD.")
    parser.add_argument("--fric", type=float, default=10, help="Provide friction value for MD.")
    parser.add_argument("--N_steps", type=int, default=50000, help="Provide number of MD steps.")
    parser.add_argument("--accu", type=float, default=1, help="Set the accuracy for SCF convergence for MD.")
    parser.add_argument("--step", type=float, default=1, help="Set the timestep value for MD.")
    parser.add_argument("--mdamp", type=float, default=1, help="Set the mixer damping value for MD.")
    parser.add_argument("--fmax", type=float, default=0.15,
                        help="Set the maximum force for optimization termination for Opt.")
    parser.add_argument("--steps", type=int, default=5000, help="Provide number of optimization steps for Opt.")
    parser.add_argument("--maxstep", type=float, default=0.05,
                        help="Set the maximum step length in optimization for Opt.")
    parser.add_argument("--mixer_damping", type=float, default=0.4,
                        help="Set the mixer damping value for SCF convergence.")
    parser.add_argument("--scf_accuracy", type=float, default=1, help="Set the SCF accuracy for convergence.")
    parser.add_argument("--max_iterations", type=int, default=5000, help="Set the maximum number of SCF iterations.")
    parser.add_argument("--verbosity", type=int, default=1, help="Set the verbosity level for SCF convergence.")
    parser.add_argument("--charge", type=int, default=0, help="Set the total charge of the system.")
    parser.add_argument("--logfile", type=str, default='opt.log', help="Log file name for the  FIRE optimizer output")
    parser.add_argument("--multiplicity", type=int, default=1, help="Set the spin multiplicity of the system.")
    parser.add_argument("--optimizer", choices=['FIRE', 'BFGS'], default='FIRE',
                        help="Select optimizer for structure optimization.")
    parser.add_argument("--restart", type=int, default=None, help="Whether to restart the MD simulation.")
    args = parser.parse_args()

    main(args)
