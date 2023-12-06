import numpy as np
import scipy.constants as const
from pydantic import BaseModel, conlist
from typing import List, Union, Callable, Any, Optional

hbar = 1
docs = """
11/30/2023 Changes (mostly for readability):
   * Added Replica class that keeps track of the state of each replica
   * Removed fixed array of replicas, now has initial value of min_replicas (target replicas)
   * Added N_target to be explicit about
   * Cull dead replicas from the array at the end of each step
"""
class Replica(BaseModel):
    alive: bool = False
    xs_array: conlist(float, min_length=0) = []
    xs_bumps: conlist(float, min_length=0) = []
    def __str__(self):
        return f"alive: {self.alive}, xs_array: {self.xs_array}"
    def __repr__(self):
        return f"alive: {self.alive}, xs_array: {self.xs_array}"
    
class QMC(BaseModel):
    V: Callable[[float], float] = None     # potential function V(x) associated with the specific quantum system we are modeling
    particle_count: int = 2                # number of particles to simulate
    dim: int = 1                   # number of dimensions (number of spatial dimensions = 1 x number of particles)
    min_replicas: int = 1000       # minimum number of replicas
    max_replicas: int = 20000      # maximum number of replicas
    delta_tau: float = 0.1         # time step size (Δτ = 0.1)
    xmin: float = -20              # minimum value of the spatial coordinate (xmin = −20)
    xmax: float = 20               # maximum value of the spatial coordinate (xmax = 20)
    bins: int = 200                # number of spatial bins for sorting the replicas (only used during 'Counting' to plot the ground state wave function)
    seed: int = 42                 # seed value for the random number generators (for repeatability)
    mass: float = 1                # mass of the particle (m = 1)
    E_ref: float = None            # reference energy (E_ref = 0)
       
    DEBUG: bool = False            # debug flag
    # replicas: Optional[List[List[Union[bool, float]]]] = None
    replicas: List[Replica] = None
    
    N: int = 500                   # the count of ALIVE replicas; initially equal to min_replicas
    N_prev: int = 500              # the count of ALIVE replicas from the previous step
    N_target: int = 500            # the target number of replicas (initially equal to min_replicas)
    alpha: float = 0.2             # Used in our modified E_ref calculation

    def __init__(self, **data):
        """ initialize the simulation based on input values (things are implicitly set via the super class __init__)"""
        super().__init__(**data)  # Call the super class __init__

        if self.V is None: raise ValueError("Potential function V(x) must be provided, dammit")
        np.random.seed(seed=self.seed)

        xs_array = [0.0 for _ in range(self.particle_count-1)]

        # replicas array contains live replicas (the dead are culled at the end of each step)
        self.replicas = [Replica(alive=True, 
                                 last_particle_pos=0.0, 
                                 xs_array=xs_array,
                                 xs_bumps=xs_array) for i in range(self.min_replicas)]
        if self.DEBUG: self.print_replicas()
        
        # making sure that our initial count is equal to the minimum number of replicas
        self.N = self.min_replicas # initially equal to min_replicas
        self.N_prev = self.min_replicas 
        self.N_target = self.min_replicas # let's do this once and hold onto the initial value forever (don't update down below)
            
        # The initial value of the reference energy E_ref is the potential energy at the initial position of the replicas.
        self.Calculate_E_ref()
    
    def replica_tot_pot(self, xs_array):
        """ 
        Calculates the total potential energy of the system for a replica.
        The xs_array contains the relative distances between the particles.
        """
        if self.particle_count < 2: raise ValueError("There should be at least 2 particles")

        V_tot = 0.0

        # LOOP1: Calculate potential for directly stored relative distances
        for x in xs_array: V_tot += self.V(x)

        # LOOP2: Calculate potential for 'derived relative distances'
        for i in range(self.particle_count-1):
            for j in range(i+1, self.particle_count-1):
                # Derived distances calculated using differences between distances in xs_array
                derived_distance = xs_array[i] - xs_array[j]
                V_tot += self.V(derived_distance)

        return V_tot
    
    def CountReplicas(self):
        """ NOTE: RUN THIS AFTER CullDeadReplicas .... Counts the number of alive replicas. """
        self.N_prev = self.N    # Save the previous count
        self.N = len(self.replicas) # simple count of the number of replicas (which are all alive)
        if self.DEBUG: print(f"CountReplicas: N: {self.N}  N_prev: {self.N_prev}  N_target: {self.N_target}")
        
    def CullDeadReplicas(self):
        """ Reap all of the dead replicas by creating a new replica array using only those that are alive. """
        self.replicas = [replica for replica in self.replicas if replica.alive]

    def Calculate_V_avg(self):
        """ 
        Calulates the average potential for each replicas.
        
        First calculate the total potential for each replica, 
           and then find the average potential across all replicas.
        """
        V_tot = 0.0
        for replica in self.replicas:
            V_tot += self.replica_tot_pot(replica.xs_array)
        V_avg = V_tot / self.N
        return V_avg    
    
    def Calculate_KineticEnergy(self):
        """ Charlie made me do it! """        
        # for replica in self.replicas:
        #     xs_bumps = replica.xs_bumps
        #     KE = 0.0
        #     for i in range(len(xs_bumps)):
        #         KE += .5 * self.mass * (xs_bumps[i]/self.delta_tau) **2
        #     KE_tot += KE
            
        xs_tot = 0.0
        delta_tot = [0.0 for _ in range(self.particle_count-1)]

        for replica in self.replicas:
            xs_bumps = replica.xs_bumps
            for ii in range(len(delta_tot)):
                delta_tot[ii] += xs_bumps[ii]
        
        delta_tot = np.multiply(delta_tot,  (1 / self.N))
        KE_tot = 0.0
        for ii in range(len(delta_tot)):
            KE_tot += .5 * self.mass * (delta_tot[ii]/self.delta_tau) **2

        return KE_tot 

    def Calculate_E_ref(self):
        """ 
        Calulates the reference energy of the system -- updates each step.
        """
        
        if self.N == 0 or self.N is None:
            raise ValueError("No alive replicas found")

        # to calculate E_ref, we either:
        #    1. use eqn 2.33 E_ref = <V> for the the original E_ref, and then use 2.36 to update E_ref, OR 
        #    2. use eqn 2.33 E_ref = <V> for the the original E_ref, and then using alpha instead of (hbar / self.delta_tau)
        
        
# to calculate the next value of E_ref we can either 
# using an updated value of V average every time, and "tune" E_ref based on the ratio between N and N_prev, or
# calculate V_average at the very beginning and "tune" E_ref based on the ratio between N and N_target

        
        if self.E_ref is None: 
            V_avg = self.Calculate_V_avg()
            self.E_ref = V_avg    # eqn 2.33
        else: 
            use_eqn_2_35 = True # select between 2.35 and 2.36 to compute E_ref
            if use_eqn_2_35:
                V_avg = self.Calculate_V_avg()
                alpha = (const.hbar / self.delta_tau)
                self.E_ref = V_avg - self.alpha * (1 - self.N / self.N_prev) # eqn 2.35
            else:
                self.E_ref += self.alpha * (1 - self.N / self.N_target)  # eqn 2.36

            self.E_ref += self.Calculate_KineticEnergy()
    
    def Walk(self):
        """ 
        This walks the relative positions between the particles for all replicas. 
        """
        prefactor = np.sqrt(self.delta_tau)
        for replica in self.replicas: # iterate over all replicas
            # add a random amount to the relative distances between particles
            for i, x in enumerate(replica.xs_array): 
                replica.xs_bumps[i] = prefactor * np.random.normal()
                replica.xs_array[i] += replica.xs_bumps[i]

    def print_replicas(self):
        for ii,replica in enumerate(self.replicas):
            print(f"replica[{ii}]: {replica}")
                
                        
    def Branch(self):
        """ 
        This is the Birth/Death decision branch for each replica. 
        m_n = min[W(x), u_n]
        """
        dtau_over_hbar = self.delta_tau/hbar
        
        for i in range(self.N): # iterate over all replicas
            replica = self.replicas[i]
#           W = np.exp(-dtau_over_hbar * (self.replica_tot_pot(xs_array) - self.E_ref)) # eqn 2.16
            W = 1 - ((self.replica_tot_pot(replica.xs_array) - self.E_ref) * dtau_over_hbar) # eqn 2.29
            
            m_n = min(int(W + np.random.uniform()), 3)
            if m_n == 0:  # Kill the replica
                if self.DEBUG: print(f" Killing replica {i}  W: {W}  m_n: {m_n}  E_ref: {self.E_ref:.4f}")
                replica.alive = False

            else:  # For m_n == 2 or 3, replicate the current replica
                count_copies = m_n - 1  # Number of copies to make
                while count_copies > 0:
                    if self.DEBUG: print(f" Replica {i} - Duplicate {count_copies} times  W: {W}  m_n: {m_n}  E_ref: {self.E_ref:.4f}")
                    # Create a new replica and add it to the list
                    new_replica = Replica(alive=replica.alive, 
                                          xs_array=replica.xs_array.copy(),
                                          xs_bumps=replica.xs_bumps.copy())
                    self.replicas.append(new_replica)
                    count_copies -= 1
                        
                if count_copies > 0:
                    raise ValueError("Not enough dead replicas to make copies")
        
    def Binning(self):
        """ 
        Divides the range [x_min,x_max] into bin_count equally sized bins. 
        For each replica, we compute the average position of this replica. 
        
        This is to be done after the system has stabilized. 
        """
        # roll through all of the replicas and calculate the centroid position of each 
        centroid_array = []
        for replica in self.replicas: # iterate over all replicas
            centroid_array.append(self.centroid(replica.xs_array))
            
        hist_array = np.histogram(centroid_array, bins=self.bins, range=[self.xmin, self.xmax])
        return hist_array, centroid_array
    
    def step(self):
        """ Steps the simulation forward 1 delta-t step and returns <V> and N. """
        self.Walk()
        self.Calculate_E_ref()
        self.Branch()
        self.CullDeadReplicas()
        self.CountReplicas()
        if self.DEBUG: print("-"*80)
        # nothing is returned, but class variables have been updated
