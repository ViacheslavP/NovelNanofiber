import numpy as np

# methods and procedures for creating various chains
class atomic_state(object):
    def __init__(self, noa: int, zpos, campl):
        if (len(zpos) != noa) or (len(campl) != noa):
            raise TypeError("Wrong Atomic State")
        self.noa = noa
        self.zpos = np.asarray(zpos, dtype=np.float64)
        if not np.vdot(campl, campl) == 0:
            self.campl = campl / np.vdot(campl, campl)
        else:
            self.campl = campl

    def __add__(self, other):
        return self.merge_atomic_states(other, distance=0.5)

    def __sub__(self, other):
        return self.merge_atomic_states(other, distance=0.25)

    def __mul__(self, other):
        return self.merge_atomic_states(other, distance=np.floor(2 * 1e8 / 780)/2 + 1/4)

    def __truediv__(self, other):
        return self.merge_atomic_states(other, distance=np.floor(2 * 1e8 / 780)/2)

    def merge_atomic_states(self, bstate: object, distance, to_end=True) -> object:
        n1, n2 = self.noa, bstate.noa
        add_dist = 0
        if to_end:
            add_dist = self.zpos[-1]
        zpos = np.concatenate((self.zpos, bstate.zpos + distance + add_dist), axis=None)
        campl = np.concatenate((self.campl, np.exp(2j * (distance + add_dist) * np.pi) * bstate.campl), axis=None)
        return atomic_state(n1 + n2, zpos, campl)

# Creating a simple chain of noa atoms with period d
def new_chain(noa: int, d, random=False):
    zpos = d * np.arange(noa)
    if random:
        zpos = noa * d * np.sort(np.random.rand(noa))
    campl = np.ones_like(zpos, dtype=np.complex)
    return atomic_state(noa, zpos, campl)


# Create an excited chain with shared excitation
def new_echain(noa: int, d) -> object:
    chain = new_chain(noa, d)
    for i in range(noa):
        chain.campl[i] *= np.exp(-2j * np.pi * chain.zpos[i])
    return chain


# Creating an chain that has no excitation
def new_mchain(noa: int, d, random=False) -> object:
    chain = new_chain(noa, d, random)
    chain.campl = np.zeros_like(chain.campl, dtype=np.complex)
    return chain
