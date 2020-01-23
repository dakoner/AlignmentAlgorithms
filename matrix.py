alphabet = "A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","B","Z","X","*"

def makeIdentity(match, mismatch):
  matrix = {}
  for i in range(len(alphabet)):
    matrix[alphabet[i],alphabet[i]] = match
    for j in range(i+1,len(alphabet)):
      matrix[alphabet[i],alphabet[j]] = mismatch
      matrix[alphabet[j],alphabet[i]] = mismatch
  return matrix

def read_matrix(filename):
  f = open(filename)
  matrix = {}
  aa_row = None
  for line in f.readlines():
    if line.startswith("#"):
      continue
    elif line.startswith("   "):
      aa_row = line.split()
    else:
      aa_col = line[0]
      vals = line[1:].split()
      for i, val in enumerate(vals):
        matrix[aa_row[i], aa_col] = int(val)

  return matrix

blosum62 = read_matrix('blosum62.txt')
