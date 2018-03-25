#!/usr/bin/python


# A bunch of modular matrix arithmetic.
# This is pretty disorganized for now.


import collections
from sympy import Matrix


def finiteInverse(a, p):
    """Compute 1/a mod p using the Euclidean algorithm."""
    x = [a % p, 1, 0]
    y = [p, 0, 1]
    while x[0] != 1:
        z = divmod(y[0], x[0])
        x, y = [z[1], y[1] - z[0] * x[1], y[2] - z[0] * x[2]], x
    return x[1] % p


def modularized(m, p):
    return m.applyfunc(lambda x: x % p)


def modularInverse(m, p):
    return modularized(m.adjugate() * finiteInverse(m.det(), p), p)


class GLFinite(collections.Iterator):
    """
    Iterator for everything in GL(n, F_p), where p is prime.

    This doesn't actually check that p is prime.
    Ensuring that is foisted on the user.

    """

    def __init__(self, n, p):
        self.n = n
        self.n2 = self.n ** 2
        self.p = p
        self._current = Matrix(n, n, lambda i, j: 0)

    def next(self):
        # There's almost certainly a smarter way to do this.
        while True:
            self._current[0] += 1
            for i in range(self.n2):
                if self._current[i] == self.p:
                    self._current[i] = 0
                    try:
                        self._current[i + 1] += 1
                    except IndexError:
                        raise StopIteration(
                            "GL(%s,F_%s) exhausted" %
                            (self.n, self.p))
                else:
                    break
            if self._current.det() % self.p != 0:
                return Matrix(self._current)
