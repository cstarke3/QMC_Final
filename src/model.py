import numpy as np
import scipy.constants as const

class QMC: #I'm sure I am missing stuff on referencing inside the class and whatnot, but the structure should be right

    def __init__(self, V, particles, dim=1, min_replicas=500, max_replicas=2000, max_steps=1000, delta_t=0.1, xmin=-20, xmax=20, bins=200):
        self._V = V # potential function V(x) associated with the specific quantum system we are modeling
        self._dim = dim # number of dimensions (D=1 for this exercise)
        self._min_replicas = min_replicas # minimum number of replicas
        self._max_replicas = max_replicas # maximum number of replicas
        #self._seed = seed # seed value for the random number generator
        self._max_steps = max_steps # maximum number of time steps to run the simulation (τ0 = 1000)
        self._delta_t = delta_t # time step size (Δτ = 0.1)
        self._xmin = xmin # minimum value of the spatial coordinate (xmin = −20)
        self._xmax = xmax # maximum value of the spatial coordinate (xmax = 20)
        self._bins = bins # number of spatial bins for sorting the replicas during their sampling
        self.mass = 1
        self.hbar = 1 #changed to one

        self.particles = particles # number of particles
        self.replicas = dict()
        self.Energy = []
        self.total_count = [min_replicas]
        self.count = np.zeros(self.particles)



    def initialize(self, particles):
      partial_Energy = []

      for i in range(self.particles): # Create replica arrays assigned to individual particles
          self.replicas[i] = np.zeros((2, self._max_replicas))
          self.replicas[i][1,0:self._min_replicas] = 1

      for i in range(self.particles):
          self.count[i] = sum(self.replicas[i][1,:])

      for i in range(self.particles): # Finds the Energy from Eq, 2.32 and 2.33
          partial_Energy.append(self.Average_Potential(i))

      self.Energy.append(np.sum(partial_Energy))


    def Average_Potential(self, i): #Calulates the average potential for one set of replicas from
        sum = 0
        for k in range(int(self.count[i])): #I might be indexing wrong here, might be count-1; Could also switch this to N_max and not have to do a sort function, right now we need one though
            sum += self._V(self.replicas[i][0, k])
        return (1/self.count[i])*sum

    def Sort(self, i): # Sorts live replicas to the front of the array
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


    def Count_func(self): #Counts the number of alive particles in the replica set associated with a particular particle
        for k in range(self.particles):
          self.count[k] = sum(self.replicas[k][1,:])

        self.total_count.append(np.sum(self.count))


    def Walk(self, i): #Walks every replica associated with particle i according to Eq. 2.30
        for k in range(int(self.count[i])):
            self.replicas[i][0,k] = self.replicas[i][0,k] + np.sqrt(self.hbar*self._delta_t/self.mass)*np.random.randn()

    def Branch(self, i): #conducts the branching of the replicas
        for k in range(int(self.count[i])):
            index = 1
            W = np.absolute(1 - ((self._V(self.replicas[i][0,k]) - self.Average_Potential(i))/self.hbar)*self._delta_t)
            W = int(W + np.random.uniform())
            m = min(W, 3)
            if m == 0:
               self.replicas[i][1, k] = 0
            elif m == 2:
                self.replicas[i][1, int(self.count[i]+index)] = 1
                self.replicas[i][0, int(self.count[i]+index)] = self.replicas[i][0,k]
                index += 1
            elif m == 3:
                self.replicas[i][1, int(self.count[i]+index)] = 1
                self.replicas[i][0, int(self.count[i]+index)] = self.replicas[i][0,k]
                self.replicas[i][1, int(self.count[i]+index+1)] = 1
                self.replicas[i][0, int(self.count[i]+index+1)] = self.replicas[i][0,k]
                index += 2

    def Energy_Step(self, i): #finds the next energy value based on the previous energy value
        self.Energy.append(self.Energy[i-1] + (self.hbar/self._delta_t)*(1-(self.total_count[i]/self.total_count[i-1])))

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
                self.Sort(k)
                self.Walk(k)
                self.Branch(k)
                self.Count_func()
              self.Energy_Step(t)

        return self.Energy[-1], self.Energy
