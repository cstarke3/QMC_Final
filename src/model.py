import numpy as np
import scipy.constants as const
from pydantic import BaseModel

hbar = 1

class QMC(BaseModel):
    V: float               # potential function V(x) associated with the specific quantum system we are modeling
    particles: int         # number of particles to simulate
    dim: int = 1           # number of dimensions (D=1 for this exercise)
    min_replicas: int = 500 # minimum number of replicas
    max_replicas: int = 2000 # maximum number of replicas
    max_steps: int = 1000 # maximum number of time steps to run the simulation (τ0 = 1000)
    delta_t: float = 0.1 # time step size (Δτ = 0.1)
    xmin: float = -20 # minimum value of the spatial coordinate (xmin = −20)
    xmax: float = 20 # maximum value of the spatial coordinate (xmax = 20)
    bins: int = 200 # number of spatial bins for sorting the replicas (only used during 'Counting' to plot the ground state wave function)
    seed: int = 42 # seed value for the random number generators (for repeatability)
    mass: float = 1 # mass of the particle (m = 1)

    self.replicas = dict()
    self.Energy = []
    self.total_count = [min_replicas*particles]
    self.count = np.zeros(particles)

    def initialize(self, particles):
        """ Initialize the simulation. """
        partial_Energy = []

        for i in range(self.particles): # Create replica arrays assigned to individual particles
            self.replicas[i] = np.zeros((2, self._max_replicas))
            self.replicas[i][1,0:self._min_replicas] = 1

        for i in range(self.particles):
            self.count[i] = sum(self.replicas[i][1,:])

        for i in range(self.particles): # Finds the Energy from Eq, 2.32 and 2.33
            partial_Energy.append(self.Average_Potential(i))

        self.Energy.append(np.sum(partial_Energy))


    def Average_Potential(self, i): 
        """ Calulates the average potential for one set of replicas from """
        sum = 0
        for k in range(int(self.count[i])): #I might be indexing wrong here, might be count-1; Could also switch this to N_max and not have to do a sort function, right now we need one though
            sum += self._V(self.replicas[i][0, k])
        return (1/self.count[i])*sum

    def Sort(self, i): 
        """ Sorts live replicas to the front of the array"""
        for k in range(int(self.count[i])): #sets dead replicas postions to 0
            if self.replicas[i][1,k] == 0:
                self.replicas[i][0,k] = 0

        transposed_array = list(zip(*self.replicas[i]))
        def custom_sort(item):
          return item[1]

        sorted_transposed_array = sorted(transposed_array, key=custom_sort)
        dummy_array = list(zip(*sorted_transposed_array))
        final_array = [dummy_array[0][::-1], dummy_array[1][::-1]]

        for k in range(self._max_replicas):
          self.replicas[i][0,k] = final_array[0][k]
          self.replicas[i][1,k] = final_array[1][k]#sorts in descending order to group live replicas to the front of the array


    def Count_func(self): 
        """ Counts the number of alive particles in the replica set associated with a particular particle """
        for k in range(self.particles):
          self.count[k] = sum(self.replicas[k][1,:])



    def Walk(self, i): 
        """ Diffusion step: Walks every replica associated with particle i according to Eq. 2.30"""
        for k in range(int(self.count[i])):
            self.replicas[i][0,k] = self.replicas[i][0,k] + np.sqrt(self._delta_t)*np.random.randn()

    def Branch(self, i): 
        """ The Birth/Death process - conducts the branching of the replicas"""
        Zed = 0
        Two = 0
        Three = 0
        index = 1
        for k in range(int(self.count[i])):
            W = 1 - ((self._V(self.replicas[i][0,k]) - self.Energy[-1]/self.hbar))*self._delta_t
            W = int(W + np.random.uniform())
            m = min(W, 3)
            if m == 0:
                """ kill the replica """
                self.replicas[i][1, k] = 0
                self.replicas[i][0, k] = 0
                Zed += 1
            elif m == 2:
                """ create 1 replica with the same position as the parent replica """
                self.replicas[i][1, int(self.count[i]+index)] = 1
                self.replicas[i][0, int(self.count[i]+index)] = self.replicas[i][0,k]
                index += 1
                Two += 1
            elif m == 3:
                """ create 2 replicas with the same position as the parent replica"""
                self.replicas[i][1, int(self.count[i]+index)] = 1
                self.replicas[i][0, int(self.count[i]+index)] = self.replicas[i][0,k]
                self.replicas[i][1, int(self.count[i]+index+1)] = 1
                self.replicas[i][0, int(self.count[i]+index+1)] = self.replicas[i][0,k]
                index += 2
                Three += 1
        #print(Zed, Two, Three)
        #print(Two - Zed)

    def Energy_Step(self, i): #finds the next energy value based on the previous energy value
        self.total_count.append(np.sum(self.count))
        self.Energy.append(self.Energy[i-1] + (self.hbar/self._delta_t)*(1-(self.total_count[i]/self.total_count[i-1])))
        print(self.Energy[-1], self.total_count[-1])

    def Test(self, i): #only run this if i is greater than a set number, 10 maybe?
        flag = 0
        sum = 0
        for k in range(10):
            sum += np.absolute(self.total_count[i] - self.total_count[k]) #Sum the count differences for the current count compared to the last 9
        avg = sum/10 #Average the count difference
        if avg <= 5: #5 here is the tolerance, not sure what a good value is or if there is a just a better way in general to run this test
            flag = 1
        return flag

    def Ground_State_Energy(self):
        self.initialize(self.particles)
        flag = 0
        if flag == 0:
            for t in range(1, self._max_steps):
              for k in range(self.particles):
                if t > 10:
                  flag = self.Test(k)
                self.Walk(k)
                self.Branch(k)
                self.Count_func()
                self.Sort(k)
              self.Energy_Step(t)

        return self.Energy[-1], self.Energy