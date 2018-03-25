#!/usr/bin/python

"""Methods for manipulating presentations of complement groups."""

from .braidextras import tally


def simplifyWord(word):
    """
    Perform cancellations in a word.

    Input word should be a list of integers;
        [1, 1, -2] means x_1^2 x_2^{-1}

    """
    not_stuck = True
    while not_stuck:
        not_stuck = False
        prev = None
        new_word = []
        for x in word:
            if prev == -x:
                prev = None
                not_stuck = True
            else:
                if prev is not None:
                    new_word.append(prev)
                prev = x
        # Make sure we get the last letter in too
        if prev is not None:
            new_word.append(prev)
        word = new_word
    return word


def trimPresentation(generators, relations):
    """Trim a presentation."""
    def _findIsolated(rel):
        """Find an isolated generator. Return None if not found."""
        t = tally(abs(x) for x in rel)
        for x, count in t.iteritems():
            if count == 1:
                return x
        return None

    def _simplifyRel(rel):
        """Simplify a relation."""
        ans = simplifyWord(rel)
        # Conjugates to relations are relations, so maybe we can chop off ends.
        while len(ans) > 1 and ans[0] == -ans[-1]:
            ans = ans[1:-1]
        return ans

    def _expandRel(rel, target, replacement):
        """
        Expand target (a generator) to its replacement.

        Input types should be:
            rel: a list of integers
            target: an integer
            replacement: a list of integers

        """
        inverse = [-x for x in reversed(replacement)]
        ans = []
        for x in rel:
            if x == target:
                ans.extend(replacement)
            elif x == -target:
                ans.extend(inverse)
            else:
                ans.append(x)
        return ans

    newgens = list(generators)
    oldrels = [list(y) for y in sorted(relations, key=lambda x: (-len(x), x))]
    newrels = []
    while len(oldrels) > 0:
        rel = _simplifyRel(oldrels.pop())
        # Look for a generator to replace
        g = _findIsolated(rel)
        if g is None:
            # None found? I guess we just have to take this relation.
            if len(rel) > 0:
                newrels.append(rel)
        else:
            # Determine the word to substitute
            if -g not in rel:
                rel = [-x for x in reversed(rel)]
            i = rel.index(-g)
            sub = rel[i + 1:] + rel[:i]
            # Make the substitution in all old *and* new relations
            newrels = [_expandRel(x, g, sub) for x in newrels]
            oldrels = [_expandRel(x, g, sub) for x in oldrels]
            # Drop this generator
            newgens.remove(g)
    return [newgens, [_simplifyRel(rel) for rel in newrels]]
