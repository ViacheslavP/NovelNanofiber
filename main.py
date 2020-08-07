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
Dicke, DickeP, DickeIn = get_temporal_metrics(echain, time, nof, f_max)
LeftDicke, RightDicke = get_pulses(echain, time, nof, f_max)
to_mathematica('Dicke', time, Dicke, DickeP, DickeIn, LeftDicke, RightDicke, csvpath='data/07.08.20/')


print('\n With collective state')
CollSq, CollP, CollIn = get_temporal_metrics(dumped_res, time, nof, f_max)
LeftColl, RightColl = get_pulses(dumped_res, time, nof, f_max)
to_mathematica('Coll', time, CollSq, CollP, CollIn, LeftColl, RightColl, csvpath='data/07.08.20/')

print('\n Single defect')
SdSq, SdP, SdIn = get_temporal_metrics(openres, time, nof, f_max)
LeftSd, RightSd = get_pulses(openres, time, nof, f_max)
to_mathematica('Single_defect', time, SdSq, SdP, SdIn, LeftSd, RightSd, csvpath='data/07.08.20/')

print('\n Two defects')
TdSq, TdP, TdIn = get_temporal_metrics(resonator, time, nof, f_max)
LeftTd, RightTd = get_pulses(resonator, time, nof, f_max)
to_mathematica('Two_defects', time, TdSq, TdP, TdIn, LeftTd, RightTd, csvpath='data/07.08.20/')

print('\n Absorbing from one side')
AfosSq, AfosP, AfosIn = get_temporal_metrics(openres_env, time, nof, f_max)
LeftAfos, RightAfos = get_pulses(openres_env, time, nof, f_max)
to_mathematica('Afos', time, AfosSq, AfosP, AfosIn, LeftAfos, RightAfos, csvpath='data/07.08.20/')

print('\n Absorbing env')
AenvSq, AenvP, AenvIn = get_temporal_metrics(env_res, time, nof, f_max)
LeftAenv, RightAenv = get_pulses(env_res, time, nof, f_max)
to_mathematica('Aenv', time, AenvSq, AenvP, AenvIn, LeftAenv, RightAenv, csvpath='data/07.08.20/')

