import numpy as np
#from matplotlib import pyplot as plt
from novelnnf import atomic_state, new_echain, new_mchain, get_metrics, to_mathematica

echain = new_echain(10, 1/2)
mirror = new_mchain(10, 1/2)
absorber = new_mchain(10, 1/2, random=True)

resonator = mirror - (echain - mirror)
dumped_res = mirror + (echain + mirror)
env_res = absorber + (echain + absorber)
large_resonator = mirror * (echain * mirror)

time = np.linspace(-1,8,100)
f_max = 100
nof = 1000

ResSq, ResP = get_metrics(resonator, time, nof, f_max)
DurSq, DurP = get_metrics(dumped_res, time, nof, f_max)
EnvSq, EnvP = get_metrics(env_res, time, nof, f_max)
LarSq, LarP = get_metrics(large_resonator, time, nof, f_max)

to_mathematica('figures', time, ResSq, ResP, DurSq, DurP, EnvSq, EnvP, LarSq, LarP, csvpath='data/18.05.20/')