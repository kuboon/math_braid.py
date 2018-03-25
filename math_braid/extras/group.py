#!/usr/bin/python

"""Group element abstract base class."""

from functools import reduce
from abc import ABCMeta, abstractmethod


class GroupElement(metaclass=ABCMeta):
    ######################
    # Required functions #
    ######################

    @abstractmethod
    def __eq__(self, other):
        """Equality test."""
        return NotImplemented

    @abstractmethod
    def __invert__(self):
        """Invert an element."""
        return NotImplemented

    @abstractmethod
    def __mul__(self, other):
        """Multiply self * other."""
        return NotImplemented

    ################################################
    # Helper functions, to make arithmetic easier. #
    # Usually no need to override                  #
    ################################################

    def __nonzero__(self):
        """
        Test whether an element is the identity.

        May need to be overridden, e.g. if blank initialization is unsupported.

        """
        return self.__ne__(self.__class__())
    __bool__ = __nonzero__  # Preparing for python 3.0 compatibility

    def __ne__(self, other):
        """Not-equal test. Default: be consistent with __eq__."""
        ans = self.__eq__(other)
        if ans is NotImplemented:
            return ans
        else:
            return not ans

    def __rmul__(self, other):
        """Multiply other * self. Default: be consistent with __mul__."""
        try:
            return self.__class__(other).__mul__(self)
        except BaseException:
            return NotImplemented

    def __pow__(self, other):
        """Compute self^other. Override if you have a shortcut."""
        try:
            exponent = int(other)
        except ValueError:
            return NotImplemented
        if exponent >= 0:
            return reduce(
                self.__class__.__mul__,
                [self] * exponent,
                self.__class__())
        else:
            return reduce(self.__class__.__mul__, [~self] * -exponent)
