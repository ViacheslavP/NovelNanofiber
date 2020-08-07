from scipy.sparse import lil_matrix, csr_matrix, eye
from scipy.sparse.linalg import spsolve
from .utills import counter
from .atomic_states import np
from .atomic_states import atomic_state

"""
Four solvers:
1. w/o dispersion and limited basis
2. w/o dispersion  and extended basis
3. w/ dispersion and limited basis
4. w/ dispersion  and extended basis
"""

# Initializion

hbar = 1  # dirac's constant = 1 => enrergy = omega
c = 1  # speed of light c = 1 => time in cm => cm in wavelenght
gd = 1.  # vacuum decay rate The only
lambd = 1  # atomic wavelength (lambdabar)
k = 1 / lambd
a = 2 * np.pi * 200 / 780 * lambd  # nanofiber radius in units of Rb wavelength
n0 = 1.45  # (silica)
R = 1.5 * a  # position of a chain

# dipole moments (Wigner-Eckart th.), atom properies

d01 = np.sqrt(hbar * gd / 4 * (lambd ** 3))  # |F = 0, n = 0> <- |F0 = 1> - populated level
d10 = d01
d00 = np.sqrt(hbar * gd / 4 * lambd ** 3)
d01m = np.sqrt(hbar * gd / 4 * lambd ** 3)
d1m0 = d01m

# Validation (all = 1 iff ful theory, except SINGLE_RAMAN)

RADIATION_MODES_MODEL = 0  # = 1 iff assuming our model of radiation modes =0 else
L = 2 * np.pi
DDI = 1
VACUUM_DECAY = 0  # = 0 iff assuming only decay into fundamental mode, =1 iff decay into fundamental and into radiation
PARAXIAL = 1  # = 0 iff paraxial, =1 iff full mode
SINGLE_RAMAN = True
FIX_RANDOM = True
SEED = 5
RABI = 1e-16  # 4#4+1e-16 #Rabi frequency
DC = 0  # Rabi detuning
SHIFT = 0

FULL = 0
FIRST = 1  # = 0 iff assuming full subtraction
SECOND = 1

EPS = (2 * np.pi * 6.065) / (384 * 1e6)  # \gamma / \omega
from .bmode import exact_mode

m = exact_mode(1 / lambd, n0, a)
VG = m.vg  # group velocity
VPH = m.vp
K = m.k
E = m.E
dE = m.dE
Ez = m.Ez


def sigma_matrix(chain: object, omega):
    zs = 2 * np.pi * chain.zpos
    noa = len(zs)
    rpos = np.asarray([[abs(zs[i] - zs[j]) for i in range(noa)] for j in range(noa)], dtype=np.complex)
    rposf = (1/1.09) * rpos / (2 * np.pi) # (lambda_0 / lambda_wg)
    D1 = RADIATION_MODES_MODEL * ((DDI * 1 - 1j * rposf - rposf ** 2) / ((rposf + np.identity(noa)) ** 3) \
                                  * np.exp(1j * 2 * np.pi * rposf)) * (np.ones(noa) - np.identity(noa))
    D2 = RADIATION_MODES_MODEL * -1 * hbar * (
            (DDI * 3 - 3 * 1j * rposf - rposf ** 2) / (((rposf + np.identity(noa)) ** 3) * np.exp(1j * 2 * np.pi * rposf)))

    neigh = np.array(list(map(lambda x: x < L, rpos)), dtype=np.complex)
    Dz = +1 * 2j * np.pi * np.exp(1j * (1 + VPH / VG * omega * EPS) * rpos) * (1 / VG)
    DzSub = -1 * 2j * np.pi * np.exp(2j * np.pi * (rposf)) * (RADIATION_MODES_MODEL * neigh)
    Di = np.zeros([noa, noa, 3, 3], dtype=np.complex)

    em = -E(R)
    ep = PARAXIAL * -dE(R) * 1
    ez = PARAXIAL * Ez(R)

    emfi = np.array([em, ez, ep], dtype=np.complex)
    emfjc = np.conjugate(np.array([-em, ez, -ep], dtype=np.complex))

    epfi = np.conjugate(np.array([-ep, ez, -em], dtype=np.complex))
    epfjc = np.array([ep, ez, em], dtype=np.complex)

    embi = np.array([-em, ez, -ep], dtype=np.complex)
    embjc = np.conjugate(np.array([em, ez, ep], dtype=np.complex))

    epbi = np.conjugate(np.array([ep, ez, em], dtype=np.complex))
    epbjc = np.array([-ep, ez, -em], dtype=np.complex)

    # Perp mode - projected onto polarization vectors

    emfjcSub = np.conjugate(np.array([-em, 0, 0], dtype=np.complex))
    epfjcSub = np.array([0, 0, em], dtype=np.complex)

    embjcSub = np.conjugate(np.array([em, 0, 0], dtype=np.complex))
    epbjcSub = np.array([0, 0, -em], dtype=np.complex)

    """
    Forward and Backward (and its sub-duality!) matrices are symmetric. 
    Moreover: they're cylindricaly symmetric!
    """
    forwardSub = np.outer(emfi, emfjcSub) + np.outer(epfi, epfjcSub)
    backwardSub = np.outer(embi, embjcSub) + np.outer(epbi, epbjcSub)

    forward = np.outer(emfi, emfjc) + np.outer(epfi, epfjc)
    backward = np.outer(embi, embjc) + np.outer(epbi, epbjc)

    gd = 1 + 8 * d00 * d00 * np.pi * ((1 / VG - 1 ) * (abs(em) ** 2 + abs(ep) ** 2 + abs(ez) ** 2) +  abs(ez) ** 2)


    for i in range(noa):
        for j in range(noa):
            """
            ________________________________________________________
            Guided modes interaction (and self-energy for V atom)
            ________________________________________________________

            Using next notations:
                epfjc - vector guided mode of Electric field 
                                           in Plus polarisation (c REM in 176)
                                  propagating Forward
                                           in J-th atom position
                                          and Conjugated components

            forward - part of Green's function for propagating forward (!)
            backward - guess what


                                               e+            e0           e-

            """

            # Full modes - exact solutions


            zi = zs[i]
            zj = zs[j]

            if abs(zi - zj) < 1e-12:
                Di[i, j, :, :] = 0.5 * d00 * d00 * (
                        (forward + backward) * Dz[i, i] - (forwardSub + backwardSub) * DzSub[i, j])

            elif zi < zj:
                Di[i, j, :, :] = (forward * Dz[i, j] - forwardSub * DzSub[i, j]) * d00 * d00

            elif zi > zj:
                Di[i, j, :, :] = (backward * Dz[i, j] - backwardSub * DzSub[i, j]) * d00 * d00

            """
            _________________________________________________________
            +Vacuum interaction (No self-energy)
            _________________________________________________________
            """

            Di[i, j, 0, 0] += d01m * d1m0 * (D1[i, j] - 0 * D2[i, j])

            Di[i, j, 0, 1] += d01m * d00 * (0 * D2[i, j]);

            Di[i, j, 0, 2] += d01m * d10 * (0 * D2[i, j]);

            Di[i, j, 1, 0] += d00 * d1m0 * (0 * D2[i, j]);

            Di[i, j, 1, 1] += d00 * d00 * (D1[i, j] + 1 * D2[i, j]);

            Di[i, j, 1, 2] += d00 * d10 * (0 * D2[i, j]);

            Di[i, j, 2, 0] += d01 * d1m0 * (-0 * D2[i, j]);

            Di[i, j, 2, 1] += d01 * d00 * (-0 * D2[i, j]);

            Di[i, j, 2, 2] += d01 * d10 * (D1[i, j] - 0 * D2[i, j])


            if i==j:
                Di[i, i, :, :] = -0.5j*gd

    return Di

def super_matrix_ext(chain: object, omega):
    dm = sigma_matrix(chain, omega)
    noa = len(chain.zpos)
    dim = noa + 2 * noa * (noa - 1)

    sm = lil_matrix((dim, dim), dtype=np.complex)
    for initial in range(dim):
        sm[initial, initial] = 1.
        if initial < noa:  # If initial state has no raman atoms (i.e. |+++e+++>)
            for s in range(noa - 1):  # s is final excited atom (rel. position)

                sr = noa + 2 * (noa - 1) * initial + 2 * s

                ne = s
                if s >= initial:
                    ne += 1  # ne is excited atom position (true)

                sm[initial, ne] = dm[initial, ne, 2, 2]  # Transition to rayleigh
                sm[initial, sr] = dm[initial, ne, 2, 1]  # Transition to Raman with m=0
                sm[initial, sr + 1] = dm[initial, ne, 2, 0]
                # assert initial != sr
                # assert initial != sr+1
        elif initial >= noa:
            k = (initial - noa) // (2 * (noa - 1))  # Raman atom
            ni = (initial - noa) % (2 * (noa - 1)) // 2  # Excited atom rel pos
            gk = (initial - noa) % (2 * (noa - 1)) % 2  # Raman state (0 or 1)
            assert initial == noa + 2 * (noa - 1) * k + 2 * ni + gk

            if ni >= k:
                ni += 1

            for s in range(noa - 1):
                ne = s
                if s >= k:
                    ne += 1

                # w/ same raman atom
                sr = noa + 2 * (noa - 1) * k + 2 * s + gk
                if sr != initial:
                    sm[initial, sr] = dm[ni, ne, 2, 2]

            # transfer exc to the Raman atom
            if ni < k:
                s = k - 1
            else:
                s = k

            sr = noa + 2 * (noa - 1) * ni + 2 * s
            assert initial != sr
            assert initial != sr + 1

            sm[initial, k] = dm[ni, k, gk, 2]  # The only way to transfer back to elastics
            sm[initial, sr] = dm[ni, k, gk, 0]
            sm[initial, sr + 1] = dm[ni, k, gk, 1]
            continue
    return csr_matrix(sm)


def super_matrix_small(chain: object, omega):
    dm = sigma_matrix(chain, omega)
    noa = len(chain.zpos)
    dim = noa

    sm = np.empty((noa, noa), dtype=np.complex)
    for initial in range(dim):
        sm[initial, initial] = 1.
        if initial < noa:  # If initial state has no raman atoms (i.e. |+++e+++>)
            for s in range(noa - 1):  # s is final excited atom (rel. position)
                ne = s
                if s >= initial:
                    ne += 1  # ne is excited atom position (true)

                sm[initial, ne] = dm[initial, ne, 2, 2]  # Transition to rayleigh
    return sm

def solver4(chain: object, freqs, rabi=0, dc=0.00001):
    noa = len(chain.zpos)
    dim = noa + 2*noa*(noa-1)
    nof = len(freqs)
    instate = np.pad(chain.campl, (0, 2*noa*(noa-1)), 'constant', constant_values=(0, 0))
    scV = np.empty((dim, nof), dtype=np.complex)
    oned = eye(dim, dim, format='csr', dtype=np.complex)
    lc = (chain.zpos[0] - chain.zpos[-1]) > 1e4
    if not lc:
        Sigma = super_matrix_ext(chain, 0)
    for i,om in enumerate(freqs):
        if lc:
            Sigma = super_matrix_ext(chain, om)
        omg = -om + rabi ** 2 / (4 * (om - dc)) - 0.5j
        resolvent = (omg - 1.) * oned + Sigma
        scV[:, i], exitCode = spsolve(resolvent, instate)
        try:
            assert exitCode == 0
        except AssertionError:
            if exitCode > 0:
                print(f'Convergence not achieved. Step {i}, Number of iterations {exitCode} \n Continue ...')
            elif exitCode < 0:
                print('Something bad happened. Exitting...')
                assert exitCode == 0
        counter(i, nof)

    return scV

def solver4_step(chain: object, om, rabi=0, dc=0.00001):
    noa = len(chain.zpos)
    dim = noa + 2*noa*(noa-1)
    instate = np.pad(chain.campl, (0, 2*noa*(noa-1)), 'constant', constant_values=(0, 0))
    assert len(instate) == dim
    scV = np.empty((dim), dtype=np.complex)
    oned = eye(dim, dim, format='csr', dtype=np.complex)

    Sigma = super_matrix_ext(chain, om)
    omg = -om + rabi ** 2 / (4 * (om - dc)) - 0.5j
    resolvent = (omg - 1.) * oned + Sigma
    scV = spsolve(resolvent, instate)

    return scV

def solver3_step(chain: object, om, rabi=0, dc=0.00001):
    noa = len(chain.zpos)
    instate = chain.campl
    oned = np.eye(noa, dtype=np.complex)

    Sigma = super_matrix_small(chain, om)
    omg = -om + rabi ** 2 / (4 * (om - dc)) - 0.5j
    resolvent = (omg - 1.) * oned + Sigma
    scV = np.linalg.solve(resolvent, instate)

    return scV
