from model import QMC

def potential(x):
    return 2*x

val = QMC(potential, 2, 10)
Energy = val.Ground_State_Energy()
print(Energy)
