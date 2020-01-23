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
        matrix[aa_row[i], aa_col] = val

  return matrix
