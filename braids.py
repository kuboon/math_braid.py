#!/usr/bin/python

"""
Braid arithmetic.

This implementation follows the description in the paper:
J. Cha et al, "An Efficient Implementation of Braid Groups",
Advances in Cryptology: Proceedings of ASIACRYPT 2001,
Lecture Notes in Computer Science (2001), 144--156.

"""

import re

from group import GroupElement
from permutation import Permutation


class CanonicalBandPermutation(Permutation):
    """Canonical factor in the band-generator presentation."""

    @classmethod
    def createFromPair(cls, pair, n):
        """
        Given "pair"=[t,s], return the corresponding canonical factor.

        The list [t,s] should denote the generator a_{ts}.
        If a negative generator a^{-1} is represented,
            then return the permutation corresponding to Da^{-1}.

        >>> CanonicalBandPermutation.createFromPair([5, 2], 7)
        [0, 4, 2, 3, 1, 5, 6]
        >>> CanonicalBandPermutation.createFromPair([2, 5], 7)
        [6, 3, 1, 2, 0, 4, 5]

        """
        if pair[0] > pair[1]:
            # Positive: Just twist
            lst = range(0, n)
            lst[pair[0] - 1] = pair[1] - 1
            lst[pair[1] - 1] = pair[0] - 1
        elif pair[0] < pair[1]:
            # Negative: Break into an inner cycle and an outer cycle
            lst = [n - 1] + range(0, n - 1)
            lst[pair[0] - 1] = (pair[1] - 2) % n
            lst[pair[1] - 1] = (pair[0] - 2) % n
        return cls(lst)
    
    def tau(self, power=1):
        """
        Transformation A --> t(A) where AD = Dt(A).

        That is, t(A) = D^{-1} A D.
        This function computes t^{power}(A).

        >>> x = CanonicalBandPermutation([0, 4, 2, 3, 1, 5, 6])
        >>> d = CanonicalBandPermutation([6, 0, 1, 2, 3, 4, 5])
        >>> x.tau()
        [0, 1, 5, 3, 4, 2, 6]
        >>> x.tau(3) == d ** -3 * x * d ** 3
        True
        >>> x.tau(0) == x
        True

        """
        return self.__class__([
            (self.table[(i - power) % self.n] + power) % self.n
            for i in xrange(0, self.n)])

    def numTranspositions(self):
        """
        Count the band generators (transpositions) required to write this.

        >>> b = LeftBraid([1, 2, 1, 3, 3, 2, 3, 1, 2], 4)
        >>> all(len(x.getTranspositions()) == x.numTranspositions() for x in b.a)
        True

        """
        return len([None for i in xrange( 0, self.n ) if self.table[i] < i])

    def getTranspositions(self):
        """
        Write a canonical factor as band generators.

        >>> b = LeftBraid([1, 2, 1, 3, 3, 2, 3, 1, 2], 4)
        >>> c = [(LeftBraid(x.getTranspositions(), 4), x) for x in b.a]

        Test that we get canonical factors back.
        >>> all(x[0].p == 0 for x in c)
        True
        >>> all(x[0].k == 1 for x in c)
        True

        Test that we get the same canonical factor.
        >>> all(x[0].a[0] == x[1] for x in c)
        True

        """
        ans = []
        for i in range(self.n - 1, -1, -1):
            if self.table[i] < i:
                ans.append([i + 1, self.table[i] + 1])
        return ans

    def computeCycles(self):
        """
        Compute the maxima of a permutation's descending cycles.

        This method is a mutator.

        >>> one = CanonicalBandPermutation([0, 2, 1, 3, 4])
        >>> one.computeCycles()
        >>> one.cycles
        [0, 2, 2, 3, 4]
        >>> two = CanonicalBandPermutation([6, 3, 1, 2, 0, 4, 5])
        >>> two.computeCycles()
        >>> two.cycles
        [6, 3, 3, 3, 6, 6, 6]

        """
        # Break the abstraction barrier for some speed
        self.cycles = range(0, self.n)
        for i in xrange(self.n - 1, -1, -1):
            if self.table[i] < i:
                self.cycles[self.table[i]] = self.cycles[i]
        return

    def meet(self, other):
        """
        Compute the meet, self /\ other:
        The max of all canonical factors smaller than both self and other.

        >>> ident = CanonicalBandPermutation()
        >>> one = CanonicalBandPermutation([0, 4, 2, 3, 1, 6, 5])
        >>> two = CanonicalBandPermutation([0, 4, 3, 2, 1, 5, 6])
        >>> one.meet(two)
        [0, 4, 2, 3, 1, 5, 6]
        >>> two.meet(one)
        [0, 4, 2, 3, 1, 5, 6]
        >>> one.meet(one) == one
        True
        >>> two.meet(two) == two
        True
        >>> one.meet(ident) == ident
        True

        """
        # These safeguards might not be necessary
        # Time cost is about 1 part in 50
        try:
            other = self.__class__(other)
        except NotImplementedError:
            return NotImplemented

        # Shortcut for an identity element
        if self.n == 0 or other.n == 0:
            return self.__class__()
        if self.n != other.n:
            return NotImplemented

        # Compute descending-cycle maxima
        self.computeCycles()
        other.computeCycles()

        # This part isn't written exactly as given in the paper
        # That version seemed to have a few confusing redundancies
        # Specifically, switching between 1...n and n...1 unnecessarily
        order = range(0, self.n)
        order.sort(key = lambda x: (self.cycles[x], other.cycles[x]))
        order.reverse()

        j = order[0]
        cycles = [j] * self.n
        for x in order[1:]:
            if self.cycles[j] != self.cycles[x] or other.cycles[j] != other.cycles[x]:
                j = x
            cycles[x] = j

        # Convert the cycles back into a permutation
        prev = [-1] * self.n
        lst = [-1] * self.n
        for i in xrange(0, self.n):
            if prev[cycles[i]] < 0:
                lst[i] = cycles[i]
            else:
                lst[i] = prev[cycles[i]]
            prev[cycles[i]] = i

        return self.__class__(lst)


class LeftBraid(GroupElement):
    """
    Implements a braid in left-greedy normal form, e.g. D^p A_1 A_2 ... A_k

    Properties:
        n: braid width (number of strands)
        p: power of fundamental element D
        k: number of canonical factors (implemented via a getter)
        a: list of canonical factors

    Some tests that are a little ridiculous to write out:
        >>> LeftBraid([1, 2, 1, 3, 2, 1, 4, 3, 2, 1, 5, 4, 3, 2, 1], 6) ** 3 * LeftBraid([3, 1, -4, 1, -5, 2, 5, 3, 5], 6)**3 == LeftBraid([1, 2, 3, 2, 4, 3, 2, 1, 5, 4, 3, 2, 1, 1, 2, 1, 4, 3, 2, 1, 5, 1, 2, 3, 2, 4, 5, 4, 3, 2, 1, 3, 2, 4, 5, 4, 2, 1, 3, 4, 3, 5, 4, 3, 2, 1, 2, 3, 4, 3, 2, 5, 4, 3, 2, 1, 1, 2, 3, 5], 6)
        True
        >>> LeftBraid([1, 2, 1, 3, 2, 1, 4, 3, 2, 1, 5, 4, 3, 2, 1], 6) ** 2 *  LeftBraid([3, 1, -4, 1, -5], 6)**3 == LeftBraid([1, 2, 1, 3, 2, 1, 4, 3, 5, 4, 1, 2, 4, 3, 2, 1, 5, 4, 3, 2, 1, 1, 4, 5, 1, 5, 1, 5, 4, 3, 1, 1, 1], 6)
        True
        >>> LeftBraid([1, 2, 1, 3, 2, 1, 4, 3, 2, 1, 5, 4, 3, 2, 1], 6) ** 2 *  LeftBraid([3, 1, -4, 1, -5], 6)**3 == LeftBraid([1, 2, 1, 3, 2, 1, 4, 3, 5, 4, 1, 2, 4, 3, 2, 1, 5, 4, 3, 2, 1, 1, 4, 5, 1, 5, 1, 5, 4, 3, 1, 1, 2], 6)
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
            >>> LeftBraid([1, 2, -1, -2], 5).getPermutation()
            [2, 0, 1, 3, 4]
            >>> LeftBraid([-3], 5).getPermutation()
            [0, 1, 3, 2, 4]

        Test identity elements
            >>> LeftBraid([], 5).getPermutation()
            [0, 1, 2, 3, 4]
            >>> LeftBraid().getPermutation()
            []

        Check that we get the right number of factors
            >>> LeftBraid([1, 2, -1, -2], 5).k
            2
            >>> LeftBraid([1, 2, -1, 2], 5).k
            4

        """
        super(LeftBraid, self).__init__(*args, **kwargs)
        
        # Easy: Copy braid properties
        if isinstance(obj, LeftBraid):
            self.n = obj.n
            self.p = obj.p
            self.a = list(obj.a)
        elif isinstance(obj, basestring):
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
            if isinstance(obj[0], LeftBraid.CanonicalFactor):
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
                for i in xrange(len(bandgens) - 2, -1, -1):
                    if bandgens[i+1][0] < bandgens[i+1][1]:
                        self.p -= 1
                    LeftBraid._tauband(bandgens[i], self.n, self.p)
                # Don't forget to check the first (leftmost) generator!
                if bandgens[0][0] < bandgens[0][1]:
                    self.p -= 1
                self.a = [LeftBraid.CanonicalFactor.createFromPair(x, self.n) for x in bandgens]
            # Now, self.a is guaranteed to be a list of canonical factors
            # Sort (sort of) this list
            l = len(self.a)
            d = LeftBraid.d(self.n)
            dcycles = [self.n - 1] * self.n
            leftmost = -1
            rightmost = l - 2
            while leftmost < rightmost:
                newleft = rightmost
                for j in xrange(rightmost, leftmost, -1):
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
                            self.a[j+1] = LeftBraid.CanonicalFactor()
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

        >>> one = LeftBraid([1], 5)
        >>> two = LeftBraid([2], 5)
        >>> one * two == LeftBraid([1, 2], 5)
        True

        """
        # Shortcut for identity elements
        if not self:
            try:
                return LeftBraid(other)
            except NotImplementedError:
                return NotImplemented
        if other is 1 or not other:
            return LeftBraid(self)
        # Ensure compatible braids
        if not isinstance(other, LeftBraid) or self.n != other.n:
            return NotImplemented
        # Combine information and construct the product
        a = [x.tau(other.p) for x in self.a] + other.a
        return LeftBraid(a, n = self.n, p = self.p + other.p)

    def __invert__(self):
        """
        Inverse of a braid.

        >>> x = LeftBraid([5, 1, -2, 4, 3, -1, -2], 6)
        >>> not x * ~x
        True

        """
        # Shortcut for identity elements
        if not self:
            return LeftBraid(self)
        # Transform and reverse the list of canonical factors
        d = LeftBraid.d(self.n)
        k = self.k
        a = [(~self.a[i] * d).tau(-self.p - i - 1) for i in xrange(k - 1, -1, -1)]
        # Construct the inverse
        return LeftBraid(a, n = self.n, p = -self.p - k)

    def __eq__(self, other):
        """
        Equality test.
        A quirk: identity elements from different B_n are considered equal.

        Test the Artin relations
        >>> LeftBraid([1, 3], 4) == LeftBraid([3, 1], 4)
        True
        >>> LeftBraid([1, 2, 1], 3) == LeftBraid([2, 1, 2], 3)
        True

        More tests, just for good measure
        >>> LeftBraid([1, 4, 4, 1], 7) == LeftBraid([4, -5, 1, 1, 5, 4], 7)
        True
        >>> LeftBraid([1, 4, 4, 1], 7) == LeftBraid([4, -5, 1, 1, 6, 4], 7)
        False

        """
        if other is 1 or not other:
            return not self
        if not isinstance(other, LeftBraid):
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
        >>> LeftBraid._tauband(a, 5)
        >>> a
        [1, 2]
        >>> LeftBraid._tauband(a, 7, -3)
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
        if not cls._d.has_key(n):
            cls._d[n] = cls.CanonicalFactor([n - 1] + range(0, n - 1))
        return cls._d[n]

    def _get_k(self):
        """Return the length of the canonical factor list."""
        return len(self.a)
    k = property(_get_k)
    
    def twist(self, i):
        """
        Make a new braid with one twist added.

        >>> a = LeftBraid([3, 2, -3, 1, 2], 5)
        >>> a.twist(2) == a
        False
        >>> a.twist(2).twist(-2) == a
        True
        >>> a.twist(2).twist(-3) == a * LeftBraid([2, -3], 5)
        True

        This tests that it leaves the original unchanged
        >>> a == LeftBraid([3, 2, -3, 1, 2], 5)
        True

        """
        if self.n == 0:
            return NotImplemented
        if 0 < i < self.n:
            p = self.p
            a = self.a + [LeftBraid.CanonicalFactor.createFromPair([i + 1, i], self.n)]
        elif -self.n < i < 0:
            p = self.p - 1
            a = [x.tau(-1) for x in self.a] + [LeftBraid.CanonicalFactor.createFromPair([-i, 1 - i], self.n)]
        return LeftBraid(a, self.n, p)

    ###########
    # Display #
    ###########

    def getPermutation(self):
        """A diagnostic function."""
        if self.n == 0:
            return Permutation()
        d = LeftBraid.d(self.n)
        right = reduce(lambda x, y: x * y, self.a, LeftBraid.CanonicalFactor(range(0, self.n)))
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

Braid = LeftBraid
class _BraidConstructorIndex(object):
    def __getitem__(self, key):
        return lambda obj=None, p=None, *args, **kwargs: LeftBraid(obj, *args, n=key, p=p, **kwargs)
B = _BraidConstructorIndex()

if __name__ == '__main__':
    import doctest
    doctest.testmod()
