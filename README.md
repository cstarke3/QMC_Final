# QMC Final

This project implements the concepts put forth in the paper. 

## Setup

### Download the code

From the commandline, in a directory of your choosing, do the following:

```bash
git clone https://github.com/cstarke3/QMC_Final.git
cd ./QMC_Final
```


### Create and Activate your local environment

To create a local environment to run this simulation, and assuming `python >= 3.11` is installed, do the following:
* on a **Mac/Linux machine**:

```bash
python -m venv qmc
source qmc/bin/activate
```

* on a **Windows machine**:

```bash
python -m venv qmc
qmc\Scripts\activate
```

Next, install all of the required packages using the following command:  

```bash
pip install -r requirements.txt
```

or upgrade, if you have already installed once: 

```bash
pip install --upgrade --force-reinstall -r requirements.txt`
```

## Running the simulation

The primary entry-point is `qmc_cli.py`, and once your environment is set up, it is invoked like this

```bash
python src/qmc_cli.py
```

which will spit out a brief help screen. To get a more detailed help output, use the `-h` option:

```bash
python src/qmc_cli.py -h
```

To run the simulation with 1000 replicas at 100 steps for 2 particles, you would simply do:

```bash
python src/qmc_cli.py -n 2
```

### Other noteable args:
```

  -n PARTICLES  the number of particles to simulate (default: 2)
  -m MIN_REPLICAS the minimum number of replicas to use during the simulation (default: 500)
  -x MAX_REPLICAS the maximum number of replicas to use during the simulation (default: 3000)
  -s STEPS the number of timesteps to use during the simulation (default: 100)
  -r RANDOM set the random seed value (default: 42)
  -t set the random seed based on the current timestamp (default: varies)
  -p plot the data (default: False)
  -d print out a bunch of stuff each time through the loop (default: False)
  -a ALPHA modify the rate at which N/N_0 impacts potential calculation (default: 0.13)
  -l loop through the algorihm for n=2-10 (default: False)
  ```
  