#!/usr/bin/python

"""Methods for simplifying factorizations in groups."""

from __future__ import print_function, division

import random
import collections
import operator
from .braidextras import complexity_mixed, lineout


def factorization_twist(factors, i):
    """
    Perform Hurwitz move <i> on the factorization <factors>.

    Caller is responsible for ensuring that abs(i) < len(factors).

    """
    if i > 0:  # positive twist
        factors[i - 1:i + 1] = [factors[i - 1] *
                                factors[i] * ~factors[i - 1], factors[i - 1]]
    else:  # negative twist
        factors[-i - 1:-i + 1] = [factors[-i],
                                  ~factors[-i] * factors[-i - 1] * factors[-i]]


class Search(collections.Iterator):
    """Base class for searches."""

    def __init__(
            self,
            factors,
            f_complexity=complexity_mixed,
            bias=2.0,
            *args,
            **kwargs):
        # Copy parameters
        self.factors = factors
        self.f_complexity = f_complexity
        self.bias = bias
        # Some properties for storing results
        self.complexity_map = {}


class WeightSearch(Search):
    def __init__(self, *args, **kwargs):
        super(WeightSearch, self).__init__(*args, **kwargs)
        self.n = len(self.factors)
        self.default_moves = set(range(1 - self.n, 0) + range(1, self.n))
        self.best = {
            'complexity': self.f_complexity(self.factors),
            'factors': list(self.factors),
            'moves_to_try': set(self.default_moves),
            'moves_to_get_here': [],
            'weight': 1.0,
            'key': str(self.factors),
        }
        # Collections of factorizations
        self.finished = {}
        self.unfinished = {
            str(self.factors): self.best,
        }

    def next(self):
        # Anything left to explore?
        if len(self.unfinished) == 0:
            raise StopIteration('Accessible factorizations exhausted.')

        # Select the top-weighted factorization; maybe would be more efficient
        # as a heap.
        curinfo = max(self.unfinished.itervalues(), key=lambda x: x['weight'])
        # Try all moves from this factorization
        for i in curinfo['moves_to_try']:
            newfactors = list(curinfo['factors'])
            factorization_twist(newfactors, i)
            # Compute complexity and weight.
            new_key = str(newfactors)
            if new_key not in self.finished:
                new_complexity = self.f_complexity(newfactors)
                new_weight = curinfo['weight'] * \
                    self.bias ** (curinfo['complexity'] - new_complexity)
                # If this is completely new to us, store some info.
                if new_key not in self.unfinished:
                    self.unfinished[new_key] = {
                        'complexity': new_complexity,
                        'factors': newfactors,
                        'moves_to_try': self.default_moves - set([-i]),
                        'moves_to_get_here': curinfo['moves_to_get_here'] + [i],
                        'weight': new_weight,
                        'key': new_key,
                    }

                    # Update our record of the best factorization.
                    if new_complexity < self.best['complexity']:
                        self.best = self.unfinished[new_key]
                # If it's already in our 'unfinished' pile, update weight and
                # moves.
                else:
                    if new_key != curinfo['key']:
                        self.unfinished[new_key]['weight'] += new_weight
                        self.unfinished[new_key]['moves_to_try'].discard(-i)
            else:
                # Indicates a bug; we should have already explored this one.
                lineout(
                    '%s, %s, %s\n' %
                    (self.finished[new_key]['moves_to_get_here'],
                     curinfo['moves_to_get_here'],
                        i))

        # We're done with this one; move it into the "finished" pile.
        self.finished[curinfo['key']] = curinfo
        del self.unfinished[curinfo['key']]

        return self.best['complexity']

    def run(self, update_interval=10, stop_at=None):
        counter = update_interval
        current_complexity = self.f_complexity(self.factors)
        # "for complexity in self" repeatedly sets complexity=self.next()
        # it loops forever or until we run out of states to examine
        try:
            for complexity in self:
                counter -= 1
                if counter == 0:
                    counter = update_interval
                    lineout('Explored %s factorizations (%s queued)' %
                            (len(self.finished), len(self.unfinished)))
                if complexity < current_complexity:
                    lineout(
                        'New best complexity: %s at %s\n    %s\n' %
                        (self.best['complexity'], len(
                            self.finished), self.best['moves_to_get_here']))
                    current_complexity = complexity
                if stop_at is not None and complexity <= stop_at:
                    break
        except KeyboardInterrupt:
            lineout('Interrupted.\n')
        lineout('Total of %s factorizations explored.\n' % len(self.finished))


class RandSearch(Search):
    def __init__(self, *args, **kwargs):
        super(RandSearch, self).__init__(*args, **kwargs)
        self.n = len(self.factors)
        self.default_moves = range(1 - self.n, 0) + range(1, self.n)
        self.positive_only = kwargs.get('positive_only', False) and True
        # Best factorization seen so far
        self.best = {
            'complexity': self.f_complexity(self.factors),
            'factors': list(self.factors),
            'moves_to_get_here': [],
        }
        # Current factorizations
        self.current = {
            'complexity': self.f_complexity(self.factors),
            'factors': list(self.factors),
            'moves_to_get_here': [],
        }

    def next(self):
        i = random.choice(self.default_moves)
        if self.positive_only:
            i = abs(i)
        newfactors = list(self.current['factors'])
        factorization_twist(newfactors, i)
        # Calculuate complexities
        new_complexity = self.f_complexity(newfactors)
        diff_complexity = self.current['complexity'] - new_complexity
        # Accept a transformation that decreases complexity
        # Or with probability bias**diff_complexity, one that increases it.
        if diff_complexity > 0 or random.random() < self.bias**diff_complexity:
            self.current['complexity'] = new_complexity
            self.current['factors'] = newfactors
            self.current['moves_to_get_here'].append(i)
            # Is this a new record?
            if new_complexity < self.best['complexity']:
                self.best['complexity'] = new_complexity
                self.best['factors'] = list(newfactors)
                self.best['moves_to_get_here'] = list(
                    self.current['moves_to_get_here'])
        return self.best['complexity']

    def run(self, update_interval=10, stop_at=None, limit=None):
        countdown = update_interval
        countup = 0
        current_complexity = self.f_complexity(self.factors)
        # "for complexity in self" repeatedly sets complexity=self.next()
        # it loops forever or until we run out of states to examine
        try:
            for complexity in self:
                countdown -= 1
                countup += 1
                if countdown == 0:
                    countdown = update_interval
                    lineout(
                        'Current complexity: %s at %s' %
                        (self.current['complexity'], countup))
                if complexity < current_complexity:
                    lineout(
                        'New best complexity: %s at %s\n    %s\n' %
                        (self.best['complexity'], countup, self.best['moves_to_get_here']))
                    current_complexity = complexity
                if stop_at is not None and complexity <= stop_at:
                    break
                if limit is not None and countup >= limit:
                    break
        except KeyboardInterrupt:
            lineout('Interrupted.\n')


def Bridge(search_type, one, two, f_complexity=complexity_mixed, **kwargs):
    if f_complexity(one) < f_complexity(two):
        smaller = one
        larger = two
        lineout('Attempting to transform second factorization into first.\n')
    else:
        larger = one
        smaller = two
        lineout('Attempting to transform first factorization into second.\n')
    smaller_inv = [~x for x in smaller]

    def f2(l):
        if len(l) == 1:
            return f_complexity(l)
        return f_complexity(map(operator.mul, smaller_inv, l))
    return search_type(larger, f_complexity=f2, **kwargs)
