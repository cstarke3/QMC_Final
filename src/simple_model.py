import numpy as np
import scipy.constants as const
from pydantic import BaseModel, conlist
from typing import List, Union, Callable, Any, Optional

hbar = 1

class QMC(BaseModel):
    V: Callable[[float], float] = None     # potential function V(x) associated with the specific quantum system we are modeling
    n_part: int = 2                # number of particles to simulate
    dim: int = 1                   # number of dimensions (number of spatial dimensions = 1 x number of particles)
    min_replicas: int = 500        # minimum number of replicas
    max_replicas: int = 2000       # maximum number of replicas
    delta_tau: float = 0.1         # time step size (Δτ = 0.1)
    xmin: float = -20              # minimum value of the spatial coordinate (xmin = −20)
    xmax: float = 20               # maximum value of the spatial coordinate (xmax = 20)
    bins: int = 200                # number of spatial bins for sorting the replicas (only used during 'Counting' to plot the ground state wave function)
    seed: int = 42                 # seed value for the random number generators (for repeatability)
    mass: float = 1                # mass of the particle (m = 1)
    E_ref: float = None            # reference energy (E_ref = 0)
       
    replicas: Optional[List[List[Union[bool, float]]]] = None
    
    N: int = 500                   # the count of ALIVE replicas; initially equal to min_replicas
    N_prev: int = 500              # the count of ALIVE replicas from the previous step
    alpha: float = 0.1             # Used in our modified E_ref calculation

    def __init__(self, **data):
        """ initialize the simulation based on input values (things are implicitly set via the super class __init__)"""
        super().__init__(**data)  # Call the super class __init__

        if self.V is None: raise ValueError("Potential function V(x) must be provided, dammit")
        np.random.seed(seed=self.seed)

        rep_pos = [0.0] # tracks the position of the last particle for each replica
        xs_array = [0.0 for _ in range(self.n_part-1)]

        # sets the first min_replicas to be alive, the rest are dead
        self.replicas = [[True] + rep_pos + xs_array if i < self.min_replicas else [False] + rep_pos + xs_array for i in range(self.max_replicas)]

        
        # making sure that our initial count is equal to the minimum number of replicas
        self.N = self.min_replicas
        self.N_prev = self.min_replicas
            
        # The initial value of the reference energy E_ref is the potential energy at the initial position of the replicas.
        self.Calculate_E_ref()
    
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
        
    def Calculate_E_ref(self):
        """ 
        Calulates the average potential across all alive replicas.
        ASSUMES that the replicas array has been sorted.
        
        First calculate the total potential for each replica, 
           and then find the average potential across all replicas.
        """
        V_tot = 0.0
        
        if self.N == 0 or self.N is None:
            raise ValueError("No alive replicas found")

        # if self.E_ref is None: # according to eqn 2.36, after the first calc, E_ref depends only on the number of alive replicas between steps
        #     for replica in self.replicas[:self.N]:
        #         xs_array = replica[2:]  # Exclude the alive/dead element and the last particle position
        #         V_tot += self.replica_tot_pot(xs_array)

        #     self.E_ref = V_tot / self.N
        # else:
        #     self.E_ref += (hbar / self.delta_tau) * (1 - self.N / self.N_prev)


        # alpha = (hbar / self.delta_tau)
        for replica in self.replicas[:self.N]:
            xs_array = replica[2:]  # Exclude the alive/dead element and the last particle position
            V_tot += self.replica_tot_pot(xs_array)
        V_avg = V_tot / self.N
        self.E_ref = V_avg + self.alpha * (1 - self.N / self.N_prev)

    def recover_absolute_positions(xs_array):
        n_particles = len(xs_array) + 1
        absolute_positions = np.zeros(n_particles)

        # Calculate absolute positions
        for i in range(1, n_particles):
            absolute_positions[i] = absolute_positions[i-1] + xs_array[i-1]

        return absolute_positions
    
    def Walk(self):
        """ 
        This walks the relative positions between the particles for each ALIVE replica. 
        ASSUMES that the replicas array has been sorted.
        """
        prefactor = np.sqrt(self.delta_tau)
        for replica in self.replicas[:self.N]: # walk all of the alive replicas
            # add a random amount to the position of the last particle AND all of the relative distances between particles
            for i, x in enumerate(replica[2:]): 
                replica[i + 2] += prefactor * np.random.normal()  # eqn 3.4

    def print_replicas(self):
        for ii,replica in enumerate(self.replicas):
            print(f"replica[{ii}]: {replica}")
                
    def Branch(self):
        """ 
        This is the Birth/Death decision branch for each replica. 
        m_n = min[W(x), u_n]
        """
        dtau_over_hbar = self.delta_tau/hbar
        
        for replica in self.replicas[:self.N]:
            xs_array = replica[2:]  # Exclude the alive/dead element and the last particle position
            W = np.exp(-dtau_over_hbar * (self.replica_tot_pot(xs_array) - self.E_ref)) # eqn 2.16
            m_n = min(int(W+np.random.uniform()), 3)
            if m_n == 0:  # Kill the replica
                replica[0] = False
            else:  # For m_n == 2 or 3, replicate the current replica
                count_copies = m_n - 1  # Number of copies to make
                for ii, target_replica in enumerate(self.replicas):
                    if count_copies == 0:
                        break
                    if not target_replica[0]:  # If the target replica is dead
                        self.replicas[ii] = list(replica)  # Make a copy
                        count_copies -= 1


    def centroid(self, xs_array):
        """
        Calculate the centroid of a system based on relative offsets from the last particle.

        :param xs_array: Relative offsets of all but the last particle from the last particle.
                        These offsets can be floating point numbers.
        :return: The centroid of the system as a floating point number.
        """
        n = len(xs_array) + 1  # Total number of particles
        # Position of the last particle is 0. Positions of other particles are their negative offsets.
        positions = [-x for x in xs_array] + [0.0]
        return sum(positions) / n

        
    def Binning(self):
        """ 
        Divides the range [x_min,x_max] into bin_count equally sized bins. 
        For each replica, we compute the average position of this replica. 
        
        This is to be done after the system has stabilized. 
        """
        # roll through all of the replicas and calculate the centroid position of each 
        centroid_array = []
        for replica in self.replicas[:self.N]:
            xs_array = replica[2:]
            centroid_array.append(self.centroid(xs_array))
            
        hist_array = np.histogram(centroid_array, bins=self.bins, range=[self.xmin, self.xmax])
        return hist_array, centroid_array
    
    def step(self):
        """ Steps the simulation forward 1 delta-t step and returns <V> and N. """
        self.Walk()
        self.Calculate_E_ref()
        self.Branch()
        self.Sort()
        self.CountAliveReplicas()
        # nothing is returned, but class variables have been updated
