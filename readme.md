# Braid Group implementation on pure python

This project is based on http://www-math.mit.edu/~seidel/geng/#download

Changes from original:
- Original code is only for python2. Mine is only for python3
- Removed Permutation and use sympy.combinatorics.Permutation
- Built directory structures and packaged for pip.
- Removed LeftBraid. just use Braid.

# Usage

```sh
pip3 install math_braid
```

```python3
from math_braind import Braid

# Identity Element
Braid([], 4)

# Artin generation. num should be from 1 to 3
Braid([1,2,3], 4)

# Band generation. num should be from 1 to 4
b = Braid([[3,1], [1,2]], 4); b

# Length of Permutation
b.k

# inverse element
~b

```
