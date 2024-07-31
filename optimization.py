import os
import time
from ase import Atoms
from tblite.interface import Calculator
from tblite.ase import TBLite
from ase.optimize import BFGS, FIRE
from ase.io.trajectory import Trajectory
from utils import traj_to_extxyz
from ase.units import Bohr, Hartree

def xtb_optim(method: str = None,
              fmax: float = 0.15,
              steps: int = 5000,
              maxstep: float = 0.05,
              structure: Atoms = None,
              file_extension: str = None,
              output_name: str = None,
              output_dir: str = '',
              mixer_damping: float = 1,
              scf_accuracy: float = 1,
              max_iterations: int = 5000,
              verbosity: int = 1,
              charge: int = None,
              multiplicity: int = None,
              optimizer: str = 'BFGS') -> (float, float, int):

    # Generate a timestamp to append to the directory name
    timestamp = time.strftime('%H:%M:%S_%d-%m-%y')

    # Fetch the job ID from the environment variable set by SLURM or a default if not found
    job_id = os.getenv('SLURM_JOB_ID', 'no_job_id')
    job_name = os.getenv('SLURM_JOB_NAME', 'default_job_name')

    # Define the directory path using the method, optimization keyword, timestamp, and job ID
    results_dir = f'results/{output_dir}GO__{job_id}__{output_name}_{job_name}/'
    os.makedirs(results_dir, exist_ok=True)

    # Save the initial structure in the newly created directory
    structure.write(f'{results_dir}{output_name}_init{file_extension}')

    # Set up the ASE-compatible TBLite calculator for the optimization
    ase_calc = TBLite(method=method, max_iterations=max_iterations, verbosity=verbosity,
                      scf={'mixer_damping': mixer_damping, 'accuracy': scf_accuracy},
                      charge=charge, multiplicity=multiplicity)
    structure.calc = ase_calc

    # Initialize the optimizer based on the provided argument
    if optimizer == 'FIRE':
        opt = FIRE(structure, trajectory=f'{results_dir}{output_name}.traj', maxstep=maxstep)
    elif optimizer == 'BFGS':
        opt = BFGS(structure, trajectory=f'{results_dir}{output_name}.traj', maxstep=maxstep)
    else:
        raise ValueError(f"Unsupported optimizer: {optimizer}")

    # Perform the optimization
    t_start = time.time()
    opt.run(fmax=fmax, steps=steps)
    t_end = time.time()
    opt_time = t_end - t_start

    # Save the final optimized structure and convert trajectory to extended XYZ format
    structure.write(f'{results_dir}{output_name}_optimized{file_extension}')
    traj_to_extxyz(f'{results_dir}{output_name}')
    n_steps = len(Trajectory(f'{results_dir}{output_name}.traj'))
    energy = structure.get_total_energy() * 23.060  # Conversion factor for energy unit if needed

    # Perform a single-point calculation with the interface Calculator on the optimized structure
    # The structure now contains the optimized positions
    numbers = structure.get_atomic_numbers()
    positions = structure.get_positions() / Bohr
    lattice = structure.cell / Bohr if structure.cell is not None else None
    periodic = structure.get_pbc()

    calc = Calculator(method, numbers, positions, lattice=lattice, periodic=periodic)
    calc.set("max-iter", max_iterations)
    calc.set("verbosity", verbosity)
    calc.set("mixer-damping", mixer_damping)
    calc.set("accuracy", scf_accuracy)

    res = calc.singlepoint()
    orbital_energies = res.get("orbital-energies") * Hartree
    orbital_occupations = res.get("orbital-occupations")

    # Extract and print orbital energies and occupations
    with open(f'{results_dir}{output_name}_orbital_info.txt', 'w') as f:
        f.write("Orbital Energies and Occupations:\n")
        f.write(" #   Occupation    Energy/Eh    Energy/eV\n")
        f.write("-------------------------------------------\n")
        for i, (energy, occupation) in enumerate(zip(orbital_energies, orbital_occupations), start=1):
            f.write(f"{i:<4} {occupation:<12.6f} {energy/Hartree:<12.6f} {energy:<12.6f}\n")

    return energy, opt_time, n_steps