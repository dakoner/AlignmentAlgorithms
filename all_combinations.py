# The number of alignments can be expressed using binomials (see https://math.stackexchange.com/questions/1694210/number-of-paths-of-length-ni-in-n-x-n-grid-with-diagonal) which also leads to a bijective correspondence, crudely implemented below in Python.

import itertools

def path_to_align(s, t, path):
  print(path)
  si = iter(s)
  ti = iter(t)
  for item in path:
    if item == "↗":
      print(next(si), next(ti))
    elif item == "→":
      print('-', next(ti))
    elif item == "↑":
      print(next(si), '-')
  print()

def paths(u, r):
    for d in range(min(u, r) + 1):
        n = u + r - d
        for d_steps in itertools.combinations(range(n), d):
            for r_steps in itertools.combinations(set(range(n)) - set(d_steps), r - d):
                yield [
                    ("↗" if i in d_steps else "→" if i in r_steps else "↑")
                    for i in range(n)
                ]

s = "SMILES"
t = "ILEAC"

for path in paths(len(s), len(t)):
  path_to_align(s, t, path)
