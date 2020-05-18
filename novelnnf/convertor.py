from .atomic_states import np

CSVPATH = 'data/fig18.05/'


def to_mathematica(fname, *argv, csvpath=CSVPATH):
    import os
    if not os.path.exists(csvpath):
        os.makedirs(csvpath)
    toCsv = np.column_stack(argv)
    np.savetxt(csvpath + fname + '.csv', toCsv, delimiter=',', fmt='%1.8f')
