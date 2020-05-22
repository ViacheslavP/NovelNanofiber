import numpy as np
#from matplotlib import pyplot as plt
from novelnnf import atomic_state, new_echain, new_mchain, get_temporal_metrics, to_mathematica

echain = new_echain(100, 1/2)
mirror = new_mchain(100, 1/2)
absorber = new_mchain(100, 1/2, random=True)

resonator = mirror - (echain - mirror)
dumped_res = mirror + (echain + mirror)
env_res = absorber + (echain + absorber)
large_resonator = mirror * (echain * mirror)

time = np.linspace(-1,8,1000)
f_max = 100
nof = 1000


print('\n Dicke chain')
Dicke, _ = get_temporal_metrics(echain, time, nof, f_max)
print('\n Rabi resonator')
ResSq, ResP = get_temporal_metrics(resonator, time, nof, f_max)
print('\n Dumped resonator')
DurSq, DurP = get_temporal_metrics(dumped_res, time, nof, f_max)
print('\n Chain i environment')
EnvSq, EnvP = get_temporal_metrics(env_res, time, nof, f_max)
print('\n Echo resonator')
LarSq, LarP = get_temporal_metrics(large_resonator, time, nof, f_max)

to_mathematica('figures', time, ResSq, ResP, DurSq, DurP, EnvSq, EnvP, LarSq, LarP, Dicke, csvpath='data/18.05.20/')