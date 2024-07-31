import os
import ase
import time
from tblite.ase import TBLite
from ase.io.trajectory import Trajectory
from utils import traj_to_extxyz
from ase.md.langevin import Langevin
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.andersen import Andersen
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.npt import NPT
from ase.md.nptberendsen import NPTBerendsen
from ase import units

def xtb_md(method: str = None,
           structure: ase.Atoms = None,
           file_extension: str = None,
           output_name: str = None,
           output_dir: str = '',
           step: float = 1,
           accu: float = 1,
           fric: float = 10,
           mdamp: float = 1,
           T: float = 300.0,
           restart: int = 0,
           N_steps: int = 10000) -> (float, float, int):

    # Generate a timestamp to append to the directory name
    timestamp = time.strftime('%H:%M:%S_%d-%m-%y')

    # Fetch the job ID from the environment variable set by SLURM or a default if not found
    job_id = os.getenv('SLURM_JOB_ID', 'no_job_id')
    job_name = os.getenv('SLURM_JOB_NAME', 'default_job_name')

    # Construct the results directory name with the method, temperature, timestamp, and job ID
    results_dir = f'results/{output_dir}MD__{job_id}__{output_name}_{job_name}/'
    os.makedirs(results_dir, exist_ok=True)

    # Save the initial structure in the newly created directory
    structure.write(f'{results_dir}{output_name}_{job_id}_init{file_extension}')

    if method == 'GFN1-xTB':
        calc = TBLite(method="GFN1-xTB", max_iterations=5000, verbosity=2, accuracy=accu, mixer_damping=mdamp)
    elif method == 'GFN2-xTB':
        calc = TBLite(method="GFN2-xTB", max_iterations=5000, verbosity=2, accuracy=accu, mixer_damping=mdamp)

    structure.calc = calc

    if restart == 1:
        # No Maxwell-Boltzmann Distribution
        pass
    else:
        # Apply Maxwell-Boltzmann Distribution
        MaxwellBoltzmannDistribution(structure, temperature_K=T)

    # Set up and run the molecular dynamics simulation
    dyn = Langevin(structure, step * units.fs, temperature_K=T, friction=fric, logfile=f'{results_dir}_{job_id}_{output_name}.log')
    # dyn = VelocityVerlet(structure, timestep=step * units.fs, logfile=f'{results_dir}_{job_id}_{output_name}.log')
    # dyn = Andersen(structure, 5 * units.fs, T, 0.02, logfile=f'{results_dir}_{job_id}_{output_name}.log')
    # dyn = NVTBerendsen(structure, 0.1 * units.fs, T, taut=0.5 * 1000 * units.fs, logfile=f'{results_dir}_{job_id}_{output_name}.log')
    # dyn = NPTBerendsen(structure, timestep=step * units.fs, temperature_K=T, taut=100 * units.fs, pressure_au=1.01325 * units.bar, taup=1000 * units.fs, compressibility_au=4.57e-5 / units.bar, logfile=f'{results_dir}_{job_id}_{output_name}.log')

    traj = Trajectory(f'{results_dir}{job_id}_{output_name}.traj', 'w', structure)
    dyn.attach(traj.write, interval=1)

    t_start = time.time()
    dyn.run(N_steps)
    t_end = time.time()
    run_time = t_end - t_start

    # Convert trajectory to extended XYZ format and save
    traj_to_extxyz(f'{results_dir}{job_id}_{output_name}')

    return run_time
