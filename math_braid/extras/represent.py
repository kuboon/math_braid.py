#!/usr/bin/python

"""Methods for manipulating representations in GL(n, F_p) for p prime."""

import itertools

from sympy import Matrix

from modular import GLFinite, modularInverse, modularized
from .dihedral import IterDihedral, Dihedral


def findRepresentations(generators, relations, n, p):
    """Find all representations of the given presentation in GL(n, F_p)."""
    ans = []
    eye = Matrix(n, n, lambda i, j: i == j and 1 or 0)
    indexmap = dict((x, i) for i, x in enumerate(generators))
    # Check all possible combinations
    # groups = [GLFinite(n, p) for x in generators]
    for images in itertools.product(GLFinite(n, p), repeat=len(generators)):
        # Check each relation
        violated = False
        for rel in relations:
            # Construct the product in this relation
            prod = Matrix(eye)
            for x in rel:
                if x < 0:
                    prod *= modularInverse(images[indexmap[-x]], p)
                else:
                    prod *= images[indexmap[x]]
            # Is the relation satisfied? (Is the product the identity?)
            # If not, stop checking for this set of generators.
            if modularized(prod, p) != eye:
                violated = True
                break
        # If any relation broke, we won't add this one.
        if not violated:
            ans.append(images)
    return ans


def findDihedral(generators, relations, n):
    """Find all homomorphisms to the dihedral group of order 2n."""
    ans = []
    eye = Dihedral()
    indexmap = dict((x, i) for i, x in enumerate(generators))
    # Check all possible combinations
    for images in itertools.product(IterDihedral(n), repeat=len(generators)):
        # Check each relation
        violated = False
        for rel in relations:
            # Construct the product in this relation
            prod = Dihedral(eye)
            for x in rel:
                if x < 0:
                    prod *= ~images[indexmap[-x]]
                else:
                    prod *= images[indexmap[x]]
            # Is the relation satisfied? (Is the product the identity?)
            # If not, stop checking for this set of generators.
            if prod != eye:
                violated = True
                break
        # If any relation broke, we won't add this one.
        if not violated:
            ans.append(images)
    return ans
