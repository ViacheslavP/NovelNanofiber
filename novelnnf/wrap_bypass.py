def get_solution_pairs(dim, nof, Sigma, ddRight, freq, gamma, rabi, dc):
    import numpy as np
    from scipy.sparse import lil_matrix, csr_matrix
    from scipy.sparse.linalg import lgmres as spsolve
    import sys
    scV = np.empty((dim, nof), dtype=np.complex)
    freq_scaled = (-freq + rabi[0] ** 2 / (4 * (freq - dc)) - 0.5j * gamma)
    oned = lil_matrix((dim, dim), dtype=np.complex)

    for i in range(dim):
        oned[i,i] = 1.

    Sigma = csr_matrix(Sigma, dtype=np.complex)
    oned = csr_matrix(oned, dtype=np.complex)


    for i,om in enumerate(freq_scaled):
        resolvent = (om-1.)*oned + Sigma
        scV[:,i], exitCode = spsolve(resolvent, ddRight)
        try:
            assert exitCode == 0
        except AssertionError:
            if exitCode > 0:
                print(f'Convergence not achieved. Step {i}, Number of iterations {exitCode} \n Continue ...')
            elif exitCode < 0:
                print('Something bad happened. Exitting...')
                assert exitCode == 0
        ist = 100 * (i+1) / nof
        sys.stdout.write("\r%d%%" % ist)
        sys.stdout.flush()
        sys.stdout.write("\033[K")

    return scV

def get_solution(dim, nof, noa, Sigma, ddRight, freq, gamma, rabi, dc):

    import numpy as np

    scV = np.empty((dim, nof), dtype=np.complex)

    from scipy.sparse import csr_matrix
    from scipy.sparse.linalg import spsolve

    print("Memory Error was found. Using iterations:")
    for i, om in enumerate(freq):
        resolvent = csr_matrix(np.eye(dim) * (om + rabi[0] ** 2 / (4 * (om - dc)) - 0.5j * gamma) + Sigma)
        scV[:, i] = spsolve(resolvent, ddRight)

    return scV



