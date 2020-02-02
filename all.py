def recurse(s, t, i, j, s1, t1):
    # print(i, s[i], j, t[j], s1, t1)
    if i == len(s):
        print(''.join(s1))
        print(''.join(t1))
        return
    if j == len(t):
        print(''.join(s1))
        print(''.join(t1))
        return
    recurse(s, t, i+1, j, s1 + [s[i]], t1 + ['-'])
    recurse(s, t, i, j+1, s1 + ['-'], t1 + [t[j]])
    recurse(s, t, i+1, j+1, s1 + [s[i]], t1 + [t[j]])
    print()
    

s = "SMILES"
t = "ILEAC"

recurse(s, t, 0, 0, [], [])
