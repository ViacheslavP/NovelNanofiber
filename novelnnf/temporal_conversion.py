
from .wave_pack import convolution2

from .atomic_states import np


def get_psi_temporal(psi0: object, time, nof, f_max, solver='4'):
    if solver=='4':
        from .dyson_solvers import solver4_step as solver
    elif solver=='3':
        pass
    elif solver=='2':
        pass
    elif solver=='1':
        pass
    else:
        raise ValueError("Solver hasn't been found")

    freq = np.linspace(-f_max, f_max, nof)
    noa = len(psi0.zpos)
    dim = noa + 2*noa*(noa-1)
    Decay = np.zeros((dim, nof), dtype=np.complex)
    for i, om in enumerate(freq):
        Decay[:, i] = solver(psi0, om)

    DecayTemp = np.zeros((dim, len(time)), dtype=np.complex)

    for i in range(dim):
        DecayTemp[i, :] = convolution2(time, freq, Decay[i, :], np.ones_like(freq))
    return time, DecayTemp


def get_metrics(setup: object, time2, nof, f_max):
    t, setup_td = get_psi_temporal(setup, time2, nof, f_max)
    noa = len(setup.zpos)
    instate = np.pad(setup.campl, (0, 2 * noa * (noa - 1)), 'constant', constant_values=(0, 0))
    setup_square = np.real(np.dot(np.conj(np.transpose(setup_td)), setup_td).diagonal())
    setup_init = abs(np.dot(np.conj(np.transpose(setup_td)), instate)) ** 2

    return setup_square, setup_init
