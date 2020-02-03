def recurse(s, t, i, j, s1, t1):
    # print(i, s[i], j, t[j], s1, t1)
    if i == len(s) and j == len(t):
        print(''.join(s1))
        print(''.join(t1))
        print()
        return
    if i < len(s):
        recurse(s, t, i+1, j, s1 + [s[i]], t1 + ['-'])
    if j < len(t):
        recurse(s, t, i, j+1, s1 + ['-'], t1 + [t[j]])
    if i < len(s) and j < len(t):
        recurse(s, t, i+1, j+1, s1 + [s[i]], t1 + [t[j]])
    

s = "SMILES"
t = "ILEAC"

recurse(s, t, 0, 0, [], [])
