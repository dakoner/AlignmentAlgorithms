  

def makeIdentity(alphabet, match, mismatch):
    matrix = {}
    for i in range(len(alphabet)):
      matrix[alphabet[i],alphabet[i]] = match
      for j in range(i+1,len(alphabet)):
        matrix[alphabet[i],alphabet[j]] = mismatch
        matrix[alphabet[j],alphabet[i]] = mismatch
    return matrix

def read_matrix(filename):
  f = open(filename)
  aa_row = None
  matrix = {}
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


class SimilarityMatrix(object):
  def __init__(self):
    self.alphabet = None
    self.matrix = {}

class DNAMatrix(SimilarityMatrix):
  def __init__(self):
    super().__init__()
    self.alphabet = "A", "C", "G", "T"

class ProteinMatrix(SimilarityMatrix):
  def __init__(self):
    super().__init__()
    self.alphabet = "A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","B","Z","X","*"
    
class Blosum62(ProteinMatrix):
  def __init__(self):
    super().__init__()
    self.matrix = read_matrix('blosum62.txt')

blosum62 = Blosum62()
