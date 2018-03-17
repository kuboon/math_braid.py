#!/usr/bin/python

"""Example implementation of permutation group S_n."""

from __future__ import print_function

from group import GroupElement

class Permutation(GroupElement):
    """
    Permutation group S_n.

    >>> x = Permutation([0, 1, 3, 2, 4])
    >>> y = Permutation([2, 3, 4, 0, 1])
    >>> z = Permutation([1, 4, 0, 3, 2])

    >>> x == y
    False
    >>> y != z
    True
    >>> x == x
    True

    >>> x ** 2 == []
    True
    >>> y ** 5 == 1
    True

    >>> x * y
    [3, 2, 4, 0, 1]
    >>> y * x
    [2, 3, 0, 4, 1]
    >>> ~z
    [2, 0, 4, 3, 1]

    >>> (x * y) * z == x * (y * z)
    True
    >>> y * z == z * y
    False

    >>> Permutation([0, 2, 1]) * 1
    [0, 2, 1]
    >>> [1,2,0] * Permutation([0, 2, 1])
    [1, 0, 2]
    """
    
    def __init__(self, obj=None, *args, **kwargs):
        """Initialize with list, None, or other Permutation."""

        super(Permutation, self).__init__(*args, **kwargs)

        if isinstance(obj, Permutation):
            # Copy another permutation
            self.n = obj.n
            self.table = list(obj.table)
        elif obj is 1 or not obj:
            # Identity element
            self.n = 0
            self.table = []
        elif isinstance(obj, list):
            # Initialize from a list
            self.n = len(obj)
            self.table = list(obj)
        else:
            # Stumped
            raise NotImplementedError

    ##############################
    # Make it behave like a list #
    ##############################

    def __len__(self):
        """Asking for the length of an element of S_n returns n."""
        return self.n

    def __getitem__(self, key):
        """Define x[i] to be the number to which x takes i."""
        if self.n == 0:
            return key
        if not (0 <= int(key) < self.n):
            raise KeyError, 'Index out of range.'
        return self.table[key]

    def __setitem__(self, key, value):
        """Assign directly to the mapping table."""
        if self.n != 0:
            if not (0 <= int(key) < self.n):
                raise KeyError, 'Index out of range.'
        self.table[key] = value

    def __iter__(self):
        """Iterator just loops through the entries of the table."""
        return iter(self.table)
    iterkeys = __iter__

    # For friendly printing
    def __str__(self):
        """Print as a mapping."""
        return str(self.table)
    __repr__ = __str__

    ####################
    # Group Arithmetic #
    ####################
    
    def __eq__(self, other):
        """Equality test."""
        try:
            other = self.__class__(other)
        except NotImplementedError:
            return NotImplemented

        if other.n == 0:
            return self.n == 0 or self.table == range(0, self.n)
        else:
            return self.table == other.table
    
    def __invert__(self):
        """Inverse of an element."""

        # Shortcut for identity
        if self.n == 0:
            return self.__class__(self)

        # Initialize a list, and then permute it
        mapping = [0] * self.n
        for i in xrange(0, self.n):
            # Break the abstraction barrier for some speed
            mapping[self.table[i]] = i

        return self.__class__(mapping)
    
    def __mul__(self, other):
        """
        Multiply two permutations.

        It's possible to multiply a permutation by a list or by the integer 1.

        """

        # Produce a permutation
        try:
            ans = self.__class__(other)
        except NotImplementedError:
            return NotImplemented

        # Shortcut for identity elements
        if self.n == 0:
            return ans
        if ans.n == 0:
            return self.__class__(self)

        # Check compatibility and then multiply
        if self.n != ans.n:
            raise TypeError, 'Incompatible operands'
        # Break the abstraction barrier for a little speed
        ans.table = map(self.table.__getitem__, ans.table)

        return ans
    
    def __nonzero__(self):
        """Nonzero test. Overridden because we can do it faster."""
        return self.n != 0 and self.table != range(0, self.n)
    __bool__ = __nonzero__

if __name__ == '__main__':
    import doctest
    doctest.testmod()
