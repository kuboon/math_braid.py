import collections

from .group import GroupElement


class Dihedral(GroupElement):
    def __init__(self, obj=None, flip=0, n=None, *args, **kwargs):
        """Initialize element s^{flip} r_{obj} in D_{2n}."""
        super().__init__(*args, **kwargs)

        if isinstance(obj, Dihedral):
            # Copy another permutation
            self.n = obj.n
            self.i = obj.i
            self.flip = obj.flip
        elif n is not None:
            # Initialize from an integer
            self.n = n
            self.i = int(obj)
            self.flip = flip
        elif obj is 1 or not obj:
            # Identity element
            self.n = 0
            self.i = 0
            self.flip = 0
        else:
            # Stumped
            raise NotImplementedError
        # Ensure the numbers lie in range
        self.flip %= 2
        if self.n != 0:
            self.i %= self.n

    # For friendly printing
    def __str__(self):
        """Print as a mapping."""
        ans = 'r_%s' % self.i
        if self.flip:
            return 's' + ans
        return ans
    __repr__ = __str__

    ####################
    # Group Arithmetic #
    ####################

    def __nonzero__(self):
        return self.i != 0 or self.flip != 0
    __bool__ = __nonzero__

    def __eq__(self, other):
        """Equality test."""
        try:
            other = self.__class__(other)
        except NotImplementedError:
            return NotImplemented

        if other.n == 0:
            return self.n == 0 or self.i == self.flip == 0
        else:
            ans = (self.n, self.i, self.flip) == (other.n, other.i, other.flip)
            print(ans)
            return ans

    def __invert__(self):
        if self.flip:
            return Dihedral(self)
        else:
            return Dihedral(-self.i, n=self.n)

    def __mul__(self, other):
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
            raise TypeError('Incompatible operands')
        if ans.flip:
            ans.i -= self.i
        else:
            ans.i += self.i
        ans.flip += self.flip
        ans.i %= ans.n
        ans.flip %= 2

        return ans


class IterDihedral(collections.Iterator):
    """ Iterator over all elements of D_2n. """

    def __init__(self, n):
        self.n = n
        self._current = Dihedral(0, flip=0, n=n)
        self._current.i -= 1

    def next(self):
        self._current.i += 1
        if self._current.i == self.n:
            if self._current.flip == 0:
                self._current.flip = 1
            else:
                raise StopIteration
            self._current.i = 0
        return Dihedral(self._current)
