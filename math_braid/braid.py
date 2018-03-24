#!/usr/bin/python

import re
from functools import reduce
from .canonical_band_permutation import CanonicalBandPermutation

class Braid:
    """
    Implements a braid in left-greedy normal form, e.g. D^p A_1 A_2 ... A_k

    Properties:
        n: braid width (number of strands)
        p: power of fundamental element D
        k: number of canonical factors (implemented via a getter)
        a: list of canonical factors

    Some tests that are a little ridiculous to write out:
        >>> Braid([1, 2, 1, 3, 2, 1, 4, 3, 2, 1, 5, 4, 3, 2, 1], 6) ** 3 * Braid([3, 1, -4, 1, -5, 2, 5, 3, 5], 6)**3 == Braid([1, 2, 3, 2, 4, 3, 2, 1, 5, 4, 3, 2, 1, 1, 2, 1, 4, 3, 2, 1, 5, 1, 2, 3, 2, 4, 5, 4, 3, 2, 1, 3, 2, 4, 5, 4, 2, 1, 3, 4, 3, 5, 4, 3, 2, 1, 2, 3, 4, 3, 2, 5, 4, 3, 2, 1, 1, 2, 3, 5], 6)
        True
        >>> Braid([1, 2, 1, 3, 2, 1, 4, 3, 2, 1, 5, 4, 3, 2, 1], 6) ** 2 *  Braid([3, 1, -4, 1, -5], 6)**3 == Braid([1, 2, 1, 3, 2, 1, 4, 3, 5, 4, 1, 2, 4, 3, 2, 1, 5, 4, 3, 2, 1, 1, 4, 5, 1, 5, 1, 5, 4, 3, 1, 1, 1], 6)
        True
        >>> Braid([1, 2, 1, 3, 2, 1, 4, 3, 2, 1, 5, 4, 3, 2, 1], 6) ** 2 *  Braid([3, 1, -4, 1, -5], 6)**3 == Braid([1, 2, 1, 3, 2, 1, 4, 3, 5, 4, 1, 2, 4, 3, 2, 1, 5, 4, 3, 2, 1, 1, 4, 5, 1, 5, 1, 5, 4, 3, 1, 1, 2], 6)
        False
    """

    CanonicalFactor = CanonicalBandPermutation

    def __init__(self, obj=None, n=None, p=None, *args, **kwargs):
        """
        Initialize a braid.

        Input can be any of:
            * Another braid to copy (obj)
            * A power of D (p) and a list of canonical factors (obj)
            * A list of Artin generators, given as integers (obj)
            * A list of band generators, given as 2-element lists (obj)
            * The result of str(some braid)

        Properties:
            n = number of strands
            p = exponent for D in left canonical form
            a = list of canonical factors
            k = (read-only) length of list a

        Test that we at least permute the strands correctly
            >>> Braid([1, 2, -1, -2], 5).getPermutation()
            [2, 0, 1, 3, 4]
            >>> Braid([-3], 5).getPermutation()
            [0, 1, 3, 2, 4]

        Test identity elements
            >>> Braid([], 5).getPermutation()
            [0, 1, 2, 3, 4]
            >>> Braid().getPermutation()
            []

        Check that we get the right number of factors
            >>> Braid([1, 2, -1, -2], 5).k
            2
            >>> Braid([1, 2, -1, 2], 5).k
            4

        """
        super().__init__(*args, **kwargs)
        
        # Easy: Copy braid properties
        if isinstance(obj, Braid):
            self.n = obj.n
            self.p = obj.p
            self.a = list(obj.a)
        elif isinstance(obj, str):
            m = re.match('^[\s]*\[([0-9]+)\] D\^\((-?[0-9]+)\) \* \[([0-9,\[\]\s]*)\][\s]*$', obj)
            if not m:
                return NotImplemented
            m = m.groups()
            self.n = int(m[0])
            self.p = int(m[1])
            self.a = []
            if m[2]:
                m = [int(x.strip().strip('[]')) for x in m[2].split(',')]
                while m:
                    a, m = m[:self.n], m[self.n:]
                    self.a.append(self.CanonicalFactor([int(x) for x in a]))
        elif isinstance(obj, list) and n is not None:
            self.n = n
            # Quick exit for identity elements and powers of D
            if not obj:
                self.p = p or 0
                self.a = []
                return
            if isinstance(obj[0], Braid.CanonicalFactor):
                # A list of canonical factors? Copy so we can modify in place
                if p is not None:
                    self.p = p
                    self.a = list(obj)
                else:
                    raise NotImplementedError
            else:
                # Starting from a word in generators
                if isinstance(obj[0], list):
                    # Copy list of band generators so we can modify in place
                    if len(obj[0]) == 2:
                        bandgens = list(obj)
                    else:
                        raise NotImplementedError
                else:
                    # Convert Artin generators to band generators
                    if isinstance(obj[0], int):
                        bandgens = []
                        for x in obj:
                            if 0 < x < n:
                                bandgens.append([x + 1, x])
                            elif -n < x < 0:
                                bandgens.append([-x, 1 - x])
                    else:
                        raise NotImplementedError
                # Now, bandgens is guaranteed to be a list of band generators
                # We convert to canonical factors
                # First, self.p counts the occurrences of negative generators
                self.p = 0
                for i in range(len(bandgens) - 2, -1, -1):
                    if bandgens[i+1][0] < bandgens[i+1][1]:
                        self.p -= 1
                    Braid._tauband(bandgens[i], self.n, self.p)
                # Don't forget to check the first (leftmost) generator!
                if bandgens[0][0] < bandgens[0][1]:
                    self.p -= 1
                self.a = [Braid.CanonicalFactor.createFromPair(x, self.n) for x in bandgens]
            # Now, self.a is guaranteed to be a list of canonical factors
            # Sort (sort of) this list
            l = len(self.a)
            d = Braid.d(self.n)
            dcycles = [self.n - 1] * self.n
            leftmost = -1
            rightmost = l - 2
            while leftmost < rightmost:
                newleft = rightmost
                for j in range(rightmost, leftmost, -1):
                    # I know the paper by Cha et al says d * ~self.a[j]
                    # But I think our permutations mean different things
                    # And the paper without pseudocode does it this way.
                    b = (~self.a[j] * d).meet(self.a[j+1])
                    if b:
                        # Shift b one factor to the left
                        newleft = j
                        a_jplus = ~b * self.a[j+1]
                        # Keep track of identity elements
                        if a_jplus:
                            self.a[j+1] = a_jplus
                        else:
                            if rightmost == j:
                                rightmost -= 1
                            self.a[j+1] = Braid.CanonicalFactor()
                        self.a[j] = self.a[j] * b
                leftmost = newleft
            # Clean up the list of canonical factors
            # Cut out the identity elements from the right
            l = rightmost + 2
            del self.a[l:]
            while l > 0 and not self.a[-1]:
                del self.a[-1]
                l -= 1
            # Cut out copies of D from the left
            while l > 0 and self.a[0] == d:
                del self.a[0]
                l -= 1
                self.p += 1
        elif obj is 1 or not obj:
            self.n = n or 0
            self.p = 0
            self.a = []
        else:
            raise NotImplementedError

    ####################
    # Group Arithmetic #
    ####################

    def __mul__(self, other):
        """
        Multiplication of braids.

        >>> one = Braid([1], 5)
        >>> two = Braid([2], 5)
        >>> one * two == Braid([1, 2], 5)
        True

        """
        # Shortcut for identity elements
        if not self:
            try:
                return Braid(other)
            except NotImplementedError:
                return NotImplemented
        if other is 1 or not other:
            return Braid(self)
        # Ensure compatible braids
        if not isinstance(other, Braid) or self.n != other.n:
            return NotImplemented
        # Combine information and construct the product
        a = [x.tau(other.p) for x in self.a] + other.a
        return Braid(a, n = self.n, p = self.p + other.p)

    def __pow__(self, other):
        """Compute self^other."""
        try:
            exponent = int(other)
        except ValueError:
            return NotImplemented
        if exponent >= 0:
            return reduce(self.__class__.__mul__, [self] * exponent, self.__class__())
        else:
            return reduce(self.__class__.__mul__, [~self] * -exponent)

    def __invert__(self):
        """
        Inverse of a braid.

        >>> x = Braid([5, 1, -2, 4, 3, -1, -2], 6)
        >>> not x * ~x
        True

        """
        # Shortcut for identity elements
        if not self:
            return Braid(self)
        # Transform and reverse the list of canonical factors
        d = Braid.d(self.n)
        k = self.k
        a = [(~self.a[i] * d).tau(-self.p - i - 1) for i in range(k - 1, -1, -1)]
        # Construct the inverse
        return Braid(a, n = self.n, p = -self.p - k)

    def __eq__(self, other):
        """
        Equality test.
        A quirk: identity elements from different B_n are considered equal.

        Test the Artin relations
        >>> Braid([1, 3], 4) == Braid([3, 1], 4)
        True
        >>> Braid([1, 2, 1], 3) == Braid([2, 1, 2], 3)
        True

        More tests, just for good measure
        >>> Braid([1, 4, 4, 1], 7) == Braid([4, -5, 1, 1, 5, 4], 7)
        True
        >>> Braid([1, 4, 4, 1], 7) == Braid([4, -5, 1, 1, 6, 4], 7)
        False

        """
        if other is 1 or not other:
            return not self
        if not isinstance(other, Braid):
            return NotImplemented
        if self.n != other.n:
            return NotImplemented
        return self.n == other.n and self.p == other.p and self.a == other.a

    def __nonzero__(self):
        """Override the default boolean casting, since we have a fast way."""
        return self.p != 0 or self.k != 0
    __bool__ = __nonzero__


    ###########
    # Helpers #
    ###########

    @staticmethod
    def _tauband(generator, n, power=1):
        """
        Transformation: generator --> t^{power}(generator)

        Where:          generator * D = D * t(generator)
        Modifies a 2-element list in place. No error checks.

        >>> a = [1, 5]
        >>> Braid._tauband(a, 5)
        >>> a
        [1, 2]
        >>> Braid._tauband(a, 7, -3)
        >>> a
        [5, 6]
        """

        if power == 0:
            return
        generator[0] += power
        generator[1] += power
        positivity = (generator[0] > generator[1])
        generator[0] = (generator[0] - 1) % n + 1
        generator[1] = (generator[1] - 1) % n + 1
        if positivity != (generator[0] > generator[1]):
            generator.reverse()

    _d = {}
    @classmethod
    def d(cls, n):
        """The fundamental factor D (lowercase delta here) in B_n."""
        if not n in cls._d:
            cls._d[n] = cls.CanonicalFactor([n - 1] + list(range(0, n - 1)))
        return cls._d[n]

    def _get_k(self):
        """Return the length of the canonical factor list."""
        return len(self.a)
    k = property(_get_k)
    
    def twist(self, i):
        """
        Make a new braid with one twist added.

        >>> a = Braid([3, 2, -3, 1, 2], 5)
        >>> a.twist(2) == a
        False
        >>> a.twist(2).twist(-2) == a
        True
        >>> a.twist(2).twist(-3) == a * Braid([2, -3], 5)
        True

        This tests that it leaves the original unchanged
        >>> a == Braid([3, 2, -3, 1, 2], 5)
        True

        """
        if self.n == 0:
            return NotImplemented
        if 0 < i < self.n:
            p = self.p
            a = self.a + [Braid.CanonicalFactor.createFromPair([i + 1, i], self.n)]
        elif -self.n < i < 0:
            p = self.p - 1
            a = [x.tau(-1) for x in self.a] + [Braid.CanonicalFactor.createFromPair([-i, 1 - i], self.n)]
        return Braid(a, self.n, p)

    ###########
    # Display #
    ###########

    def getPermutation(self):
        """A diagnostic function."""
        if self.n == 0:
            return CanonicalBandPermutation()
        d = Braid.d(self.n)
        right = reduce(lambda x, y: x * y, self.a, Braid.CanonicalFactor(list(range(0, self.n))))
        return (d ** self.p) * right

    def numMixedTranspositions(self):
        """ Number of transpositions in mixed canonical form. """
        if self.p >= 0:
            return (self.n - 1) * self.p + sum(a.numTranspositions() for a in self.a)
        elif self.p <= -self.k:
            return (self.n - 1) * (-self.p) + sum(self.n - 1 - a.numTranspositions() for a in self.a)
        else:
            return sum(self.n - 1 - a.numTranspositions() for a in self.a[:self.p]) + sum(a.numTranspositions() for a in self.a[self.p:])

    def __str__(self):
        return '[%s] D^(%s) * %s' % (self.n, self.p, self.a)
    __repr__ = __str__


class _BraidConstructorIndex(object):
    def __getitem__(self, key):
        return lambda obj=None, p=None, *args, **kwargs: Braid(obj, *args, n=key, p=p, **kwargs)
B = _BraidConstructorIndex()


if __name__ == '__main__':
    import doctest
    doctest.testmod()
