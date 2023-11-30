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
` cd ./QMC_Final`

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
`Example: python src/qmc_cli.py -n 3 -s 10000 -t`

options:
  -h,               --help            
                        show this help message and exit
  -n PARTICLES,     --particles PARTICLES
                        the number of particles to simulate (default: 2)
  -m MIN_REPLICAS,  --min_replicas MIN_REPLICAS
                        the minimum number of replicas to use during the simulation (default: 500))
  -x MAX_REPLICAS,  --max_replicas MAX_REPLICAS
                        the maximum number of replicas to use during the simulation (default: 2000))
  -s STEPS,         --steps STEPS
                        the number of timesteps to use during the simulation (default: 1000)))
  -r RANDOM,        --random RANDOM
                        set the random seed value (default: 42))))
  -t,               --trandom         
                        set the random seed based on the current timestamp (default: varies))))
  -b BINS,          --bins 
                        the number of spatial “boxes” (nb) for sorting the replicas during their sampling (default: 100))
  -p,               --plot_it         
                        plot the data (default: False))
  -d,               --debug           
                        plot the data (default: False))
