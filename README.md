# QMC Final

This project implements the concepts put forth in the paper XXXX. 

## Setup

### Download the code

From the commandline, do the following:
`
% cd <a dir where you want this code>
`
`
% git clone git@github.com:cstarke3/QMC_Final.git
`

### Create and Activate your local environment

To create a local environment to run this simulation, and assuming `python >= 3.11` is installed, do the following on the commandline:

`python -m venv qmc`

To activate this environment on a Mac/Linux machine: `source qmc/bin/activate`

To activate this environment on a Windows machine: `qmc\Scripts\activate`

Next, install all of the required packages using the following command:  
`pip install -r requirements.txt`

or upgrade, if you have already installed once: 
`pip install --upgrade --force-reinstall -r requirements.txt`

## Running the simulation

The primary controller (in MVC parlance) is `src/qmc_cli.py`, and once your environment is set up, it is invoked with

`
python src/qmc_cli.py
`

