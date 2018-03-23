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
            lst = list(range(0, n))
            lst[pair[0] - 1] = pair[1] - 1
            lst[pair[1] - 1] = pair[0] - 1
        elif pair[0] < pair[1]:
            # Negative: Break into an inner cycle and an outer cycle
            lst = [n - 1] + list(range(0, n - 1))
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
            for i in range(0, self.n)])

    def numTranspositions(self):
        """
        Count the band generators (transpositions) required to write this.

        >>> from braid import Braid
        >>> b = Braid([1, 2, 1, 3, 3, 2, 3, 1, 2], 4)
        >>> all(len(x.getTranspositions()) == x.numTranspositions() for x in b.a)
        True

        """
        return len([None for i in range( 0, self.n ) if self.table[i] < i])

    def getTranspositions(self):
        """
        Write a canonical factor as band generators.

        >>> x = CanonicalBandPermutation([0, 4, 2, 3, 1, 5, 6])
        >>> x.getTranspositions()
        [[5, 2]]

        >>> from braid import Braid
        >>> b = Braid([1, 2, 1, 3, 3, 2, 3, 1, 2], 4)
        >>> c = [(Braid(x.getTranspositions(), 4), x) for x in b.a]

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
        for i, v in enumerate(self.table):
            if v < i:
                ans.append([i + 1, v + 1])
        ans.reverse()
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
        self.cycles = list(range(0, self.n))
        for i in range(self.n - 1, -1, -1):
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
        order = list(range(0, self.n))
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
        for i in range(0, self.n):
            if prev[cycles[i]] < 0:
                lst[i] = cycles[i]
            else:
                lst[i] = prev[cycles[i]]
            prev[cycles[i]] = i

        return self.__class__(lst)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
