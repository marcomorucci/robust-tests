import numpy as np
from tabulate import tabulate


class SparsePmf:
    def __init__(self, array, lo, hi):
        self.idx = np.squeeze(np.dstack(np.where(array > 1e-20)))
        self.vals = [array[tuple(i)] for i in self.idx]
        self.shape = tuple(np.max(self.idx, 0) + 1)
        self.lo = np.array(lo)
        self.hi = np.array(hi)
        self.domain = []

    def __getitem__(self, key):
        if(len(key) != np.shape(self.idx)[1]):
            raise IndexError("Requested %d dimensions but %d needed." %
                             (len(key), np.shape(self.idx)[1]))
        r = tuple(np.array(key) - self.lo)
        a = np.where(np.squeeze(np.apply_along_axis(lambda x: np.all(x == r), 1, self.idx)))[0]
 
        if(len(a) == 0):
            return 0.0
        else:
            return self.vals[int(a)]

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        tab = [list(self.idx[i] + self.lo) + [np.round(self.vals[i], 3)]
               for i in xrange(len(self.idx))]
        if len(tab[0]) == 5:
            return tabulate(tab, headers=["TE-", "SD-", "TE+", "SD+", "P"], tablefmt='orgtbl')
        elif len(tab[0]) == 3:
            return tabulate(tab, headers=["Z-", "Z+", "P"], tablefmt='orgtbl')
        else:
            return tabulate(tab, tablefmt='orgtbl')

    def __iter__(self):
        self.i = 0
        return self

    def next(self):
        if self.i < len(self.idx):
            res = self.idx[self.i] + self.lo
            self.i += 1
            return res
        else:
            raise StopIteration

    def to_np(self):
        arr = np.zeros(self.hi - self.lo + 1)
        locs = [self.idx[:, c] for c in xrange(np.shape(self.idx)[1])]
        arr[locs] = self.vals
        return arr


class SparseSingle:
    def __init__(self):
        self.shape = (1, 1, 1, 1)

    def to_np(self):
        return np.array([[[[1.0]]]])
