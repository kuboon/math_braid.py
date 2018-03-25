from functools import reduce


class Permutation:
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
    [2, 3, 0, 4, 1]
    >>> y * x
    [3, 2, 4, 0, 1]
    >>> ~z
    [2, 0, 4, 3, 1]

    >>> (x * y) * z == x * (y * z)
    True
    >>> y * z == z * y
    False

    >>> Permutation([0, 2, 1]) * 1
    [0, 2, 1]
    """

    def __init__(self, obj=None, *args, **kwargs):
        """Initialize with list, None, or other Permutation."""

        if isinstance(obj, Permutation):
            # Copy another permutation
            self.size = obj.size
            self.array_form = list(obj.array_form)
        elif obj is 1 or not obj:
            # Identity element
            self.size = 0
            self.array_form = list()
        elif isinstance(obj, list):
            # Initialize from a list
            self.size = len(obj)
            self.array_form = list(obj)
        else:
            # Stumped
            raise NotImplementedError

    ##############################
    # Make it behave like a list #
    ##############################

    def __len__(self):
        """Asking for the length of an element of S_n returns n."""
        return self.size

    def __getitem__(self, key):
        """Define x[i] to be the number to which x takes i."""
        if self.size == 0:
            return key
        if not (0 <= int(key) < self.size):
            raise KeyError('Index out of range.')
        return self.array_form[key]

    def __setitem__(self, key, value):
        """Assign directly to the mapping table."""
        if self.size != 0:
            if not (0 <= int(key) < self.size):
                raise KeyError('Index out of range.')
        self.array_form[key] = value

    def __iter__(self):
        """Iterator just loops through the entries of the table."""
        return iter(self.array_form)
    iterkeys = __iter__

    # For friendly printing
    def __str__(self):
        """Print as a mapping."""
        return str(self.array_form)
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

        if other.size == 0:
            return self.size == 0 or self.array_form == list(range(0, self.size))
        else:
            return self.array_form == other.array_form

    def __invert__(self):
        """Inverse of an element."""

        # Shortcut for identity
        if self.size == 0:
            return self.__class__(self)

        # Initialize a list, and then permute it
        mapping = [0] * self.size
        for i in range(0, self.size):
            # Break the abstraction barrier for some speed
            mapping[self.array_form[i]] = i

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
        if self.size == 0:
            return ans
        if ans.size == 0:
            return self.__class__(self)

        # Check compatibility and then multiply
        if self.size != ans.size:
            raise TypeError('Incompatible operands')
        # Break the abstraction barrier for a little speed
        ans.array_form = list(map(list(other.array_form).__getitem__, self.array_form))

        return ans

    def __pow__(self, exponent):
        if exponent >= 0:
            return reduce(
                self.__class__.__mul__,
                [self] * exponent,
                self.__class__())
        else:
            return reduce(self.__class__.__mul__, [~self] * -exponent)

    def __nonzero__(self):
        """Nonzero test. Overridden because we can do it faster."""
        return self.size != 0 and self.array_form != list(range(0, self.size))
    __bool__ = __nonzero__


if __name__ == '__main__':
    import doctest
    doctest.testmod()
