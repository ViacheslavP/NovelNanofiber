import numpy as np
#from matplotlib import pyplot as plt
from novelnnf import atomic_state, new_echain, new_mchain, get_temporal_metrics, to_mathematica

echain = new_echain(100, 1/2)
mirror = new_mchain(100, 1/2)
absorber = new_mchain(100, 1/2, random=True)
bmirror = new_mchain(500, 1/2)

resonator = mirror - (echain - mirror)
dumped_res = mirror + (echain + mirror)
env_res = absorber + (echain + absorber)
large_resonator = bmirror * (echain * bmirror)
large_resonator_of = bmirror / (echain / bmirror)

openres = echain - (mirror + mirror)
openres_dumped = echain + (mirror + mirror)
openres_env = echain + (absorber + absorber)
openres_large = echain * (mirror + mirror)
openres_large_of= echain / (mirror + mirror)
time = np.linspace(-1,8,2000)
f_max = 100
nof = 1000


print('\n Dicke chain')
#Dicke, _, _ = get_temporal_metrics(echain, time, nof, f_max)
print('\n Rabi resonator')
#ResSq, ResP, ResIn = get_temporal_metrics(resonator, time, nof, f_max)
#to_mathematica('figures_rabiexp_res', time, ResSq, ResP, ResIn , csvpath='data/06.06.20/')
print('\n Dumped resonator')
#DurSq, DurP = get_temporal_metrics(dumped_res, time, nof, f_max)
print('\n Chain in environment')
#EnvSq, EnvP, EnvIn = get_temporal_metrics(env_res, time, nof, f_max)
#to_mathematica('res_rnd', time, EnvSq, EnvP, EnvIn, csvpath='data/29.05.20/')
print('\n Echo resonator')
#LarSq, LarP = get_temporal_metrics(large_resonator, time, nof, f_max)

#to_mathematica('figures', time, ResSq, ResP, DurSq, DurP, EnvSq, EnvP, LarSq, LarP, Dicke, csvpath='data/18.05.20/')

print('\n Dicke chain')
#Dicke, _ = get_temporal_metrics(openres, time, nof, f_max)
print('\n Rabi resonator')
#ResSq, ResP, ResIn = get_temporal_metrics(openres, time, nof, f_max)
#to_mathematica('figures_rabiexp_mir', time, ResSq, ResP, ResIn , csvpath='data/06.06.20/')
print('\n Dumped resonator')
#DurSq, DurP = get_temporal_metrics(openres_dumped, time, nof, f_max)
print('\n Chain i environment')
#EnvSq, EnvP, EnvIn = get_temporal_metrics(openres_env, time, nof, f_max)
#to_mathematica('openres_rnd', time, EnvSq, EnvP, EnvIn, csvpath='data/29.05.20/')
print('\n Echo resonator')
#LarSq, LarP = get_temporal_metrics(openres_large, time, nof, f_max)

#to_mathematica('figures_openres', time, ResSq, ResP, DurSq, DurP, EnvSq, EnvP, LarSq, LarP, Dicke, csvpath='data/18.05.20/')


LarSq, LarP, LarIn = get_temporal_metrics(large_resonator, time, nof, f_max)
LarSq1, LarP1, LarIn1 = get_temporal_metrics(large_resonator_of, time, nof, f_max)

to_mathematica('figures_distexp_5', time, LarSq, LarP, LarIn, LarSq1, LarP1, LarIn1, csvpath='data/06.06.20/')