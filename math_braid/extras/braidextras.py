#!/usr/bin/python

"""Extra functions for braids."""

from __future__ import print_function, division
from functools import reduce
import random
import sys
import math

from ..braid import Braid

#############################
# General-purpose functions #
#############################


class _lineout(object):
    """Produces a function like print, except it replaces the last line."""

    def __init__(self, dest=sys.stdout, end='\n'):
        self.dest = dest
        self.end = end
        self.lastlen = 0

    def __call__(self, s):
        print('\b \b' * self.lastlen, s, sep='', end='', file=self.dest)
        self.dest.flush()
        self.lastlen = len(s) - s.rfind(self.end) - 1

    def nextline(self):
        print(self.end, end='', file=self.dest)
        self.dest.flush()
        self.lastlen = 0


lineout = _lineout()
wordout = _lineout(end=' ')

########################
# Arithmetic functions #
########################


def tally(things):
    ans = {}
    for thing in things:
        ans[thing] = ans.get(thing, 0) + 1
    return ans


def product(lst):
    """Idiom for multiplying all the elements in a list."""
    if not lst:
        return 1
    return reduce(lambda x, y: x * y, lst)


def mean(lst):
    """
    Calculate the average of the entries in a list.

    >>> mean([1, 2, 3, 4])
    2.5

    """
    return sum(lst) / len(lst)


def stdev(lst, bessel=True):
    """
    Calculate standard deviation of the entries in a list.
    By default, assume it's a sample and apply Bessel's correction.

    >>> stdev([1, 2], bessel=False)
    0.5

    """
    m = mean(lst)
    denom = len(lst)
    if bessel:
        denom -= 1
    return math.sqrt(sum([(x - m) ** 2 for x in lst]) / denom)

#########################
# Extra braid functions #
#########################


def complexity_canonical(lst):
    """
    Complexity of a list of braids

    measured as the number of canonical factors plus abs(the exponent of D).

    """
    return sum(abs(b.p) + b.k for b in lst)


def complexity_transpositions(lst):
    """
    Complexity of a list of braids

    measured as the number of transpositions in its normal form.

    """
    return sum((b.n - 1) * abs(b.p) + sum([a.numTranspositions() for a in b.a]) for b in lst)


def complexity_mixed(lst):
    """
    Complexity of a list of braids

    measured as the number of transpositions in mixed canonical form
    (see D. Epstein, _Word Processing In Groups_ (1992), p. 198)

    """
    return sum(b.numMixedTranspositions() for b in lst)


def randomBraid(n=None):
    """Returns a random braid with 5-20 strands and 1-100 twists."""
    if n is None:
        n = random.randint(5, 20)
    k = random.randint(1, 100)
    return Braid([random.choice([1, -1]) * random.randrange(1, n)
                  for i in range(0, k)], n)

####################################################
# Geometry - braid monodromy factorization methods #
####################################################


def numComponents(factorization):
    """
    Return number of connected components in a factorization.

    Two strings are connected if some braid replaces one with the other.
    For a braid monodromy factorization corresponding to a surface S in R^4,
    the return value is the number of connected components of S.

    >>> numComponents([Braid([1],3)])
    2
    >>> numComponents([Braid([1],4)])
    3
    >>> numComponents([Braid([1],4), Braid([2],4), Braid([3],4)])
    1

    """
    if len(factorization) == 0:
        return False
    # We assign n component numbers at the beginning.
    # A component number n+i means "use the component number of string i"
    n = factorization[0].n
    components = [i for i in range(n)]

    def _getComponent(i):
        x = n + i
        while x >= n:
            x = components[x - n]
        return x
    # Loop through braids in the factorization, unioning components
    for braid in factorization:
        p = list(braid.getPermutation())
        for i in range(n):
            c1 = _getComponent(i)
            c2 = _getComponent(p[i])
            if c1 != c2:
                components[c2] = c1 + n
    return len(set(_getComponent(i) for i in range(n)))


def numBoundaryComponents(factorization):
    """
    The number of boundary components in a factorization.

    If the product is closed up by joining ends of strings together, the result
    is a collection of linked circles. This counts the number of circles.

    For a braid monodromy factorization corresponding to a surface S in R^4,
    this is the number of boundary components of S.

    >>> numBoundaryComponents([Braid([1],3)])
    2
    >>> numBoundaryComponents([Braid([1],4)])
    3
    >>> numBoundaryComponents([Braid([1],2), Braid([1],2)])
    2

    """
    if len(factorization) == 0:
        return 0
    p = list(product(factorization).getPermutation())
    n = len(p)
    components = [i for i in range(n)]
    for i in range(n):
        if components[i] >= n:
            continue
        components[i] = n + i
        x = p[i]
        while x != i:
            components[x] = n + i
            x = p[x]
    return len(set(components))


def getTwist(main_twist, conjugate_by, n):
    """Return the braid conjugate_by main_twist conjugate_by^{-1}."""
    return Braid(conjugate_by + [main_twist] +
                 [-x for x in reversed(conjugate_by)], n)


def getLoops(main_twist, conjugate_by):
    """
    Returns two (homotopy classes of) loops for a twist.

    Pass in a main_twist (integer) and conjugate_by (list of integers);
    the braid you should construct is:
        conjugate_by main_twist conjugate_by^{-1}
    This braid corresponds to an element [f] of the mapping class group
    MCG(D) of a punctured disk D. MCG(D) acts on the fundamental group
    of D (with some basepoint fixed below the punctures).

    The return value encodes the homotopy classes of [f](r_1) and [f](r_2)
    where r_1 and r_2 are loops enclosing the punctures corresponding to
    the two strings involved in main_twist.

    The return value is formatted as [loop1, loop2].
    Each loop is written as a list of gaps (in +/- 1 thru n+1) through
    which it passes. Gap i lies between punctures i and i+1.

    A positive sign indicates an approach from below.
    As an artifact of the format and the choice of basepoint,
    each loop will alternate signs and begin with a positive entry.

    There is some ambiguity about whether the action described above
    is a left action or right action. I have chosen to take it as a left
    action in this implementation, but as long as you can apply an
    appropriate automorphism of the braid group, it shouldn't matter.

    """
    paths = [[main_twist + 1, -main_twist], [main_twist + 2, -main_twist - 1]]
    for twist in reversed(conjugate_by):
        for i in range(len(paths)):
            a = []
            # Modify passes
            for x in paths[i]:
                # The cases happened to collapse just right to use abs(x)
                if abs(x) == abs(twist) + 1:
                    # Flip this inequality to reflect the picture
                    # twist < 0 : positive is a clockwise twist
                    # twist > 0 : positive is an anticlockwise twist
                    if twist < 0:
                        a.extend([x + 1, -x, x - 1])
                    else:
                        a.extend([x - 1, -x, x + 1])
                else:
                    a.append(x)
            paths[i] = a
    return paths


def loopToWord(loop):
    """
    Write a loop in gap notation using the standard generators of pi_1(D).

    The return format is a list of integers with the following meaning:
        [1, -2] = x_1 x_2^{-1}
    ...where x_i (1 <= i <= n) is the anticlockwise loop around puncture i.

    By convention, path composition is left-to-right. That is, "x_1 x_2"
    means "follow x_1, then follow x_2".

    """
    # We need the first pass to be positive (approaching from below)
    # since we assume the basepoint is below the punctures.
    if loop[0] < 0:
        raise NotImplementedError
    # By construction, passes alternates in sign and has even length.
    # So we pull off elements in pairs:
    i = iter(loop)
    ans = []
    while True:
        try:
            up = i.next()
        except StopIteration:
            break
        down = -i.next()
        # Enforce up > down
        if down == up:
            continue
        elif down > up:
            down = -down
            up = -up
        else:
            down -= 1
            up -= 1
        # Write just one call to range()
        ans.extend(range(up, down, -1))
    return ans


def getComplementGroup(twists, n):
    """
    twists -> a surface -> complement group -> presentation

    Each twist should have the format [main_twist, conjugate_by]
    (see braidextras.py)

    The presentation is given as [generators, relations]
    Each generator is an integer; i denotes a generator x_i.
    Each relation is a list of integers;
        [1, -2] means x_1 x_2^{-1} = id

    """
    if len(twists) == 0:
        return [[], []]
    generators = range(1, n + 1)
    relations = []
    for twist in twists:
        loops = getLoops(*twist)
        relations.append(loopToWord(
            loops[0]) + [-x for x in reversed(loopToWord(loops[1]))])
    return [generators, relations]

################################
# Generate some factorizations #
################################


def bigdelta(n):
    """ Return the standard factorization of big Delta. """
    return sum(([Braid([x], n) for x in range(1, k)]
                for k in range(n, 1, -1)), [])


def bigdelta_alt(n):
    """ Big Delta can also be factorized backwards. """
    return sum(([Braid([x], n) for x in range(k, 0, -1)]
                for k in range(1, n)), [])


def bigdelta2(n):
    """ Return Delta^2 in a slightly different factorization. """
    return n * [Braid([x], n) for x in range(1, n)]


if __name__ == '__main__':
    import doctest
    doctest.testmod()
