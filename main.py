import numpy as np
#from matplotlib import pyplot as plt
from novelnnf import atomic_state, new_echain, new_mchain, get_temporal_metrics, to_mathematica, get_pulses

noe = 100
nom = 100

echain = new_echain(noe, 1/2)
mirror = new_mchain(nom, 1/2)
absorber = new_mchain(nom, 1/2, random=True)
bmirror = new_mchain(10*nom, 1/2)
babsorber = new_mchain(10*nom, 1/2, random=True)

resonator = mirror - (echain - mirror)
dumped_res = mirror + (echain + mirror)
env_res = absorber + (echain + absorber)
large_resonator = bmirror * (echain * bmirror)
large_resonator_of = bmirror / (echain / bmirror)

openres = echain - (mirror + mirror)
openres2 = echain - bmirror
openres3 = echain / babsorber
openres_dumped = echain + (mirror + mirror)
openres_env = echain + (absorber + absorber)
openres_large = echain * (mirror + mirror)
openres_large_of= echain / (mirror + mirror)
time = np.linspace(-1,8,2000)
f_max = 400
nof = 4000


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
print('\n Chain in environment')
#EnvSq, EnvP, EnvIn = get_temporal_metrics(openres_env, time, nof, f_max)
#to_mathematica('openres_rnd', time, EnvSq, EnvP, EnvIn, csvpath='data/29.05.20/')
print('\n Echo resonator')
#LarSq, LarP = get_temporal_metrics(openres_large, time, nof, f_max)

#to_mathematica('figures_openres', time, ResSq, ResP, DurSq, DurP, EnvSq, EnvP, LarSq, LarP, Dicke, csvpath='data/18.05.20/')


#LarSq, LarP, LarIn = get_temporal_metrics(large_resonator, time, nof, f_max)
#LarSq1, LarP1, LarIn1 = get_temporal_metrics(large_resonator_of, time, nof, f_max)

#to_mathematica('figures_distexp_5', time, LarSq, LarP, LarIn, LarSq1, LarP1, LarIn1, csvpath='data/06.06.20/')

#Pulses

print('\n Dicke chain')
LeftDicke, RightDicke = get_pulses(echain, time, nof, f_max)

print('\n Collective')
LeftCollective, RightCollective = get_pulses(dumped_res, time, nof, f_max)

print('\n Two Wells')
LeftTW, RightTW = get_pulses(resonator, time, nof, f_max)

print('\n Mirror')
LeftMir, RightMir = get_pulses(openres, time, nof, f_max)

to_mathematica('figures_pulses_left', time, LeftDicke, LeftCollective, LeftTW, LeftMir, csvpath='data/06.06.20/')
to_mathematica('figures_pulses_right', time, RightDicke, RightCollective, RightTW, RightMir, csvpath='data/06.06.20/')