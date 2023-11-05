import numpy as np
import matplotlib as plt
import scipy

class QMC: #I'm sure I am missing stuff on referencing inside the class and whatnot, but the structure should be right

    def __init__(self, V, particles, steps):
        self.V = V # potential function
        self.particles = particles # number of particles
        self.steps = steps # number of monte carlo steps to take
        self.replicas = replicas
        self.Energy = Energy
        self.total_counts = total_counts
        self.hbar = const.hbar
        self.delta_t = 0.1
        

    def initialize(self, particles):
        N_max = 2000
        N_0 = 500
        self.replicas = dict()
        for i in range(self.particles): # Create replica arrays assigned to individual particles
            self.replicas[i] = np.zeros((2, N_max))
            self.replicas[i][1,0:N_0] = 1
        
        for i in range(self.particles): # Finds the Energy from Eq, 2.32 and 2.33
            partial_Energy[i] = Average_Potential(i, count)
        
        self.Energy[0] = np.sum(partial_Energy)

    
    def Average_Potential(self, i, count): #Calulates the average potential for one set of replicas from 
        sum = 0
        for k in range(count): #I might be indexing wrong here, might be count-1; Could also switch this to N_max and not have to do a sort function, right now we need one though
            sum += self.V(self.replicas[i][0, k])
        return (1/count)*sum
    
    def Sort(self, i): # Sorts live replicas to the front of the array
        for k in range(count): #sets dead replicas postions to 0
            if self.replica[i][1,k] == 0
                self.replica[i][0,k] = 0
        self.replica = np.sort(self.replica[i], reverse = True) #sorts in descending order to group live replicas to the front of the array

    def Count(self,i): #Counts the number of alive particles in the replica set associated with a particular particle
        count = 0
        count = sum(self.replicas[i][1,:])

    def Walk(self, i): #Walks every replica associated with particle i according to Eq. 2.30
        for k in range(count)
            self.replica[i][0,k] = self.replica[i][0,k] + np.sqrt(self.hbar*delta_t/mass)*np.random.randn()

    def Branch(self, i): #conducts the branching of the replicas
        for k in range(count):
            index = 1
            W = 1 - ((self.V(self.replicas[i][0,k]) - E_r)/self.hbar)*delta_t
            W = int(W + np.random.uniform())
            m = min(W, 3)
            if m == 0:
               self.replicas[i][1, k] = 0
            elif m == 2:
                self.replicas[i][1, count+index] = 1
                self.replicas[i][0, count+index] = self.replicas[i][0,k]
                index += 1
            elif m == 3:
                self.replicas[i][1, count+index] = 1
                self.replicas[i][0, count+index] = self.replicas[i][0,k]
                self.replicas[i][1, count+index+1] = 1
                self.replicas[i][0, count+index+1] = self.replicas[i][0,k]
                index += 2          
    
    def Energy_Step(self, i): #finds the next energy value based on the previous energy value
        self.Energy[i] = self.Energy[i-1] +(self.hbar/self.delta_t)*(1-(self.total_count[i]/self.total_count[i-1])

    def Test(self, i): #only run this if i is greater than a set number, 10 maybe?
        flag = 0
        sum = 0
        for k in range(10):
            sum += np.absolute(self.total_count[i] - self.total_count[k]) #Sum the count differences for the current count compared to the last 9
        avg = sum/10 #Average the count difference
        if avg <= 5: #5 here is the tolerance, not sure what a good value is or if there is a just a better way in general to run this test 
            flag = 1
        return flag
        
    def Ground_State_Energy(self, V, particles, steps):
        initialize(particles)
        flag = 0
        if flag == 0:
            for t in range(steps):
                Walk(t)
                Branch(t)
                Count(t)
                Sort(t)
                Energy_Step(t)
                flag = Test(t)
        return self.Energy[-1]


        
