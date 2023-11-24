import numpy as np
import scipy.constants as const
from pydantic import BaseModel, conlist
from typing import List, Union, Callable, Any, Optional

hbar = 1


class QMC(BaseModel):
    V: Optional[Callable[[float], float]] = None     # potential function V(x) associated with the specific quantum system we are modeling
    n_part: int = 2                # number of particles to simulate
    dim: int = 1                   # number of dimensions (number of spatial dimensions = 1 x number of particles)
    min_replicas: int = 500        # minimum number of replicas
    max_replicas: int = 2000       # maximum number of replicas
    max_steps: int = 1000          # maximum number of time steps to run the simulation (τ0 = 1000)
    delta_tau: float = 0.1           # time step size (Δτ = 0.1)
    xmin: float = -20              # minimum value of the spatial coordinate (xmin = −20)
    xmax: float = 20               # maximum value of the spatial coordinate (xmax = 20)
    bins: int = 200                # number of spatial bins for sorting the replicas (only used during 'Counting' to plot the ground state wave function)
    seed: int = 42                 # seed value for the random number generators (for repeatability)
    mass: float = 1                # mass of the particle (m = 1)
    E_ref: float = 0               # reference energy (E_ref = 0)
       
    replicas: Optional[List[List[Union[bool, float]]]] = None
    N: int = 500                   # the count of ALIVE replicas; initially equal to min_replicas
    N_prev: int = 500              # the count of ALIVE replicas from the previous step

    def initialize(self):
        """ Initialize the simulation based on the input values. """

        print(f"Initializing simulation with {self.n_part} particles.")
        xs_array = [0.0 for _ in range(self.n_part-1)]
        print(f"xs_array: {xs_array}")
        self.replicas = [[True] + xs_array for _ in range(self.min_replicas)]
        for ii,replica in enumerate(self.replicas):
            print(f"replica[{ii}]: {replica}")
        np.random.seed(seed=self.seed)
        # The initial value of the reference energy E_ref is the potential energy at the initial position of the replicas.
        self.E_ref = self.Average_Potential()
        print(f"Initial E_ref: {self.E_ref}")
    
    def replica_tot_pot(self, xs):
        """ Calculates the total potential energy of the system for a replica."""
        if self.n_part < 2: raise ValueError("There should be at least 2 particles")

        V_tot = 0.0

        # Calculate potential for directly stored distances
        for x in xs: V_tot += self.V(x)

        # Calculate potential for derived distances
        for i in range(self.n_part-1):
            for j in range(i+1, self.n_part-1):
                # Derived distances calculated using differences between distances in xs
                derived_distance = abs(xs[i] - xs[j])
                V_tot += self.V(derived_distance)

        return V_tot

    def Sort(self):
        """ 
        Sorts live replicas to the front of the array
        True values (alive) will come before False values (dead)
        """
        self.replicas = sorted(self.replicas, key=lambda x: x[0], reverse=True)
        
    def CountAliveReplicas(self):
        """ Counts the number of alive replicas. """
        self.N_prev = self.N    # Save the previous count
        self.N = sum([1 for replica in self.replicas if replica[0]])
        
    def Average_Potential(self):
        """ 
        Calulates the average potential across all alive replicas.
        ASSUMES that the replicas array has been sorted.
        For each replica we calculate the potential energy between all particles in a replica:
            replicas[i].V_tot = V(x41) + V(x42) + V(x43) + V(x31) + V(x32) + V(x21) 
        and then take the average of each V_tot across the replicas
            V_avg = (1/N_replicas ) sum
        """
        tot_pot_sum = 0.0
        
        if self.N == 0 or self.N is None:
            raise ValueError("No alive replicas found")

        for replica in self.replicas[:self.N]:
            xs_array = replica[1:]  # Exclude the first alive/dead element
            tot_pot_sum += self.replica_tot_pot(xs_array)

        avg_pot = tot_pot_sum / self.N
        
        self.E_ref = avg_pot

    def Walk(self):
        """ 
        This walks the relative positions between the particles for each ALIVE replica. 
        ASSUMES that the replicas array has been sorted.
        """
        prefactor = np.sqrt(self.delta_tau)
        for replica in self.replicas[:self.N]:
            # Iterate over indices and values in xs_array
            for i, x in enumerate(replica[1:]):
                # Update the value in the replica array
                replica[i + 1] += prefactor * np.random.randn()  # eqn 3.4

    def print_replicas(self):
        for ii,replica in enumerate(self.replicas):
            print(f"replica[{ii}]: {replica}")
                
    def Branch(self):
        """ 
        This is the Birth/Death decision branch for each replica. 
        m_n = min[W(x), u_n]
        """
        # prefactor = np.sqrt(self.delta_tau)
        # W = 1 - ((self._V(self.replicas[i][0,k]) - self.Energy[-1]/self.hbar))*self._delta_t
        # W = int(W + np.random.uniform())
        # m = min(W, 3)
        # NOTE: this is a guess, but we're trying to decide whether or not to kill a replica
        # based on its potential energy relative to the reference energy
        prefactor = self.delta_tau/self.hbar
        
        for replica in self.replicas[:self.N]:
            xs_array = replica[1:]  # Exclude the first alive/dead element
            tot_pot_sum += self.replica_tot_pot(xs_array)
            W = np.exp(-prefactor*(self.replica_tot_pot(xs_array) - self.E_ref))
            m_n = min(int(W+np.random.uniform()), 3)

        pass
        
    def Binning(self):
        """ 
        Divides the range [x_min,x_max] into bin_count equally sized bins. 
        For each replica, we compute the average position of this replica. 
        
        This is to be done after the system has stabilized. 
        """
        pass
    
    def step(self):
        """ Steps the simulation forward 1 delta-t step and returns <V> and N. """
        self.Sort()
        self.CountAliveReplicas()
        self.Walk()
        self.Average_Potential()
        self.print_replicas()
        self.Branch()
        return self.E_ref, self.N, self.N_prev