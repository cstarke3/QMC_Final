import numpy as np
import matplotlib as plt
import scipy

class QMC: #I'm sure I am missing stuff on referencing inside the class and whatnot, but the structure should be right

    def __init__(self, V, particles, steps):
        self.V = V # potential function
        self.particles = particles # number of particles
        self.steps = steps # number of monte carlo steps to take
        self.replicas = replicas
        

    def initialize(self, particles):
        N_max = 2000
        N_0 = 500
        self.replicas = dict()
        Energy = []
        for i in range(self.particles): # Create replica arrays assigned to individual particles
            self.replicas[i] = np.zeros((2, N_max))
            self.replicas[i][1,0:N_0] = 1
        
        for i in range(self.particles): # Finds the Energy from Eq, 2.32 and 2.33
            Energy[i] = Average_Potential(i, count)
        
        Total_Energy = np.sum(Energy)

    
    def Average_Potential(self, i, count): #Calulates the average potential for one set of replicas from 
        sum = 0
        for k in range(count): #I might be indexing wrong here, might be count-1; Could also switch this to N_max and not have to do a sort function, right now we need one though
            sum += self.V(self.replicas[i][0, k])
        return (1/count)*sum
    
    def Sort(): # Haven't though about this one yet, but going to be easier to sort the replica matrices after every branch


    def Count(self,i): #Counts the number of alive particles in the replica set associated with a particular particle
        count = 0
        count = sum(self.replicas[i][1,:])

    def Walk(self, i): #Walks every replica associated with particle i according to Eq. 2.30
        for k in range(count)
            self.replica[i][0,k] = self.replica[i][0,k] + np.sqrt(hbar*delta_t/mass)*np.randn()

    def Branch(self, ): #conducts the branching of the replicas
    
    def Energy_Step(self, ): #finds the next energy value based on the previous energy value
