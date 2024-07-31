# xTB Molecular Dynamics and Optimization Scripts

# xTB Molecular Dynamics and Optimization Scripts

This repository contains Python scripts for performing geometry optimizations and molecular dynamics (MD) simulations using the xTB semi-empirical quantum chemistry methods (`GFN1-xTB` and `GFN2-xTB`). These scripts utilize the Atomic Simulation Environment (ASE) to set up, run, and analyze simulations. The methods can handle both isolated molecules and periodic boundary conditions (PBC), making them suitable for a wide range of computational chemistry tasks.

## Dependencies

To run these scripts, you need to have the following dependencies installed:

- Python 3.6+
- ASE (Atomic Simulation Environment)
- tblite
- numpy
- argparse

## Installation

You can install these dependencies using the following commands:

```bash
pip install numpy
pip install ase
pip install tblite
```
## Dependencies

To run these scripts, you need to have the following dependencies installed:

- Python 3.6+
- ASE (Atomic Simulation Environment)
- tblite
- numpy
- argparse

## Installation
You can install these dependencies using the following commands:

```
pip install numpy
pip install ase
pip install tblite
...
```
## Usage
The main script is main.py, which you can use to run either geometry optimization or molecular dynamics simulations.


To run a geometry optimization, use the --task opt option. Below is an example command:
```
python main.py --task opt --structure your_structure_file.xyz --method GFN1-xTB --fmax 0.15 --steps 5000 --maxstep 0.05 --scf_accuracy 1e-6 --mixer_damping 0.4 --charge 0 --multiplicity 1 --optimizer FIRE
```



To run a molecular dynamics simulation, use the --task md option. Below is an example command:

```
python main.py --task md --structure your_structure_file.xyz --method GFN1-xTB --T 300 --N_steps 50000 --step 1 --accu 1e-6 --fric 10 --mdamp 0.4
```


## Command-Line Arguments

Here are the descriptions of the command-line arguments you can use:
```
--method: Choose the xTB method (GFN1-xTB or GFN2-xTB).
--task: Select the task (md for Molecular Dynamics or opt for Optimization).
--structure: Provide a structure file for optimization or MD simulation.
--pbc: Select if the structure is periodic.
--cell: Provide periodic cell dimensions.
--T: Temperature in K for MD.
--fric: Friction value for MD.
--N_steps: Number of MD steps.
--accu: Accuracy for SCF convergence for MD.
--step: Timestep value for MD.
--mdamp: Mixer damping value for MD.
--fmax: Maximum force for optimization termination.
--steps: Number of optimization steps.
--maxstep: Maximum step length in optimization.
--mixer_damping: Mixer damping value for SCF convergence.
--scf_accuracy: SCF accuracy for convergence.
--max_iterations: Maximum number of SCF iterations.
--verbosity: Verbosity level for SCF convergence.
--charge: Total charge of the system.
--multiplicity: Spin multiplicity of the system.
--optimizer: Optimizer for structure optimization (FIRE or BFGS).
--restart: Whether to restart the MD simulation (1 to restart, 0 otherwise).
```