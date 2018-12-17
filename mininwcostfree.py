## a miniature implementation of Needleman Wunsch algorithm using three
## matrices M, Ix, Iy, and a backtrace matrix

import sys, matrix, getopt, string

debug = 0
def readFASTA(file):
    f = open(file)
    sequence = []
    identifier = None
    while 1:
        line = f.readline()
        if line == "":
            break
        if line[0] == ">":
            if identifier:
                break
            identifier = line[1:]
        else:
            sequence.append(string.strip(line))

    return "".join(sequence)


## this argmax will return the highest index of the largest value
def argmax(args):
    best = None
    for i in range(len(args)):
        if best == None or args[i] > args[best]:
            best = i
    print "Argmax", args, "best index: ", i, "best value: ", args[best]
    return best

def printDPMatrix(s1, s2, m):
	print "DP matrix"
	sys.stdout.write("%4s %4s " % (" ", " "))
 	for i in range(1, len(m)):
            sys.stdout.write("%4s " % s1[i-1])
	sys.stdout.write("\n")
	
	for i in range(len(m[0])):
            if i == 0:
                sys.stdout.write("%4s " % " ")
            else:
                sys.stdout.write("%4s " % s2[i-1])

            for j in range(len(m)):
                sys.stdout.write("%4d " % m[j][i])
            sys.stdout.write("\n")

def build_align(sequence1, sequence2, M, Ix, Iy, Mbt, Ixbt, Iybt):
    aligned_seq1 = []
    aligned_seq2 = []

    ## instead of starting at largest matrix value in lower right corner,
    ## started with the largest matrix in left or bottom edge
    bestValue = M[len(sequence1)][len(sequence2)]
    bestIndex = len(sequence1), len(sequence2)
    matrix = Mbt
    
    for i in range(len(sequence1)+1):
        if M[i][len(sequence2)] > bestValue:
            bestValue = M[i][len(sequence2)]
            bestIndex = i, len(sequence2)
    for j in range(len(sequence2)+1):
        if M[len(sequence1)][j] > bestValue:
            bestValue = M[len(sequence1)][j]
            bestIndex = len(sequence1), j

    for i in range(len(sequence1)+1):
        if Ix[i][len(sequence2)] > bestValue:
            bestValue = Ix[i][len(sequence2)]
            bestIndex = i, len(sequence2)
            matrix = Ixbt
    for j in range(len(sequence2)+1):
        if Ix[len(sequence1)][j] > bestValue:
            bestValue = Ix[len(sequence1)][j]
            bestIndex = len(sequence1), j
            matrix = Ixbt

    for i in range(len(sequence1)+1):
        if Iy[i][len(sequence2)] > bestValue:
            bestValue = Iy[i][len(sequence2)]
            bestIndex = i, len(sequence2)
            matrix = Iybt
    for j in range(len(sequence2)+1):
        if Iy[len(sequence1)][j] > bestValue:
            bestValue = Iy[len(sequence1)][j]
            bestIndex = len(sequence1), j
            matrix = Iybt

    i, j = bestIndex

    reverse = lambda s: ''.join([s[i] for i in xrange(len(s)-1, -1, -1)])
    ## add any gapped suffix regions skipped by cost-free end-gap
    if i < len(sequence1):
        addseq = reverse(sequence1[i:])
        aligned_seq1.extend(addseq)
        aligned_seq2.extend(len(addseq) * "-")
    elif j < len(sequence2):
        addseq = reverse(sequence2[j:])
        aligned_seq1.extend(len(addseq) * "-")
        aligned_seq2.extend(addseq)
        
        
    ## insert suffix alignment
    
##     if debug:
##         print "Best initial matrix = ", nameset[m]
    ## nasty and ugly to deal with end-conditions
##     import pdb
##     pdb.set_trace()
    while i > 0 or j > 0:
        if debug:
            print i, j
        ## if current matrix is match, emit the relevant sequence
        if id(matrix) == id(Mbt):
            if debug:
                print "Match matrix",
            ## if either sequence has reached its left edge, gap it rather than using a negative index
            ## why would this be the case, ever?
            if i < 1:
                if debug:
                    print "emit s1 -",
                s1 = "-"
            else:
                if debug:
                    print "emit s1", sequence1[i-1],
                s1 = sequence1[i-1]
            if j < 1:
                if debug:
                    print "emit s2 -"
                s2 = "-"
            else:
                if debug:
                    print "emit s2", sequence2[j-1]
                s2 = sequence2[j-1]

            ## transition between matrices

            ## M->M
            if matrix[i][j] == 0:
                if i > 0: i -= 1
                if j > 0: j -= 1
            ## M->Ix
            elif matrix[i][j] == 1:
                if i > 0: i -= 1
                if j > 0: j -= 1
                if debug:
                    print "Transition to Ixbt"
                matrix = Ixbt
            ##M->Iy
            elif matrix[i][j] == 2:
                if i > 0: i -= 1
                if j > 0: j -= 1
                if debug:
                    print "Transition to Iybt"
                matrix = Iybt

        ## if current matrix is Ix (insert in second sequence) then emit a char from the first sequence and gap the second sequence
        elif id(matrix) == id(Ixbt):
            if debug:
                print "Ix matrix",
            s1 = sequence1[i-1]
            s2 = "-"
            if debug:
                print "emit s1", sequence1[i-1], "emit s2", "-"
            if matrix[i][j] == 0:
                if debug:
                    print "Transition to Mbt"
                matrix = Mbt
            i -= 1

        ## if current matrix is Iy (insert in first sequence) then gap the first sequence and emit a char from the second sequence
        elif id(matrix) == id(Iybt):
            if debug:
                print "Iy matrix"
            s1 = "-"
            s2 = sequence2[j-1]
            if debug:
                print "emit s1", "-", "emit s2", sequence2[j-1]
            if matrix[i][j] == 0:
                if debug:
                    print "Transition to Mbt"
                matrix = Mbt
            j -= 1

        aligned_seq1.append(s1)
        aligned_seq2.append(s2)

    aligned_seq1.reverse()
    aligned_seq2.reverse()

    return "".join(aligned_seq1), "".join(aligned_seq2)

def nw(seq1, seq2, matrix, d, e):
    l_1 = len(seq1) + 1
    l_2 = len(seq2) + 1

    M = []
    Match = []
    Ix = []
    Iy = []
    Mbt = []
    Ixbt = []
    Iybt = []
    for i in range(l_1):
        M.append( [0] * l_2) 
        Ix.append( [0] * l_2)
        Iy.append( [0] * l_2)
        Mbt.append( [0] * l_2)
        Ixbt.append( [0] * l_2)
        Iybt.append( [0] * l_2)

    for i in range(1, l_1):
##         Ix[i][0] = i * -e
##         Ix[i][0] = -d + (i-1) * -e
        Ixbt[i][0] = 1
    for i in range(1, l_2):
##         Iy[0][i] = i * -e
##         Iy[0][i] = -d + (i-1) * -e
        Iybt[0][i] = 1


    if debug:
        print "Initial"
        printDPMatrix(seq1, seq2, M)
        printDPMatrix(seq1, seq2, Ix)
        printDPMatrix(seq1, seq2, Iy)
        print ""

    for j in range(1, l_2):
        for i in range(1, l_1):
            if debug:
                print "Cell", i,j
            match = matrix[seq1[i-1],seq2[j-1]]

            set =  (  M[i-1][j-1] + match,
                     Ix[i-1][j-1] + match,
                     Iy[i-1][j-1] + match)
            setmax = argmax((set))
            M[i][j] = set[setmax]
            Mbt[i][j] = setmax

            if debug:
                print
                print "M M", seq1[i-1],seq2[j-1], "%d + %d = %d" % (M[i-1][j-1], match, M[i-1][j-1] + match)
                print "IxM", seq1[i-1], "-", "%d + %d = %d" % (Ix[i-1][j-1], match, Ix[i-1][j-1] + match)
                print "IyM", "-", seq2[j-1], "%d + %d = %d" % (Iy[i-1][j-1], match, Iy[i-1][j-1] + match )
                print "Best M", setmax, set[setmax]
                
            set =  ( M[i-1][j] - d,\
                    Ix[i-1][j] - e)
            setmax = argmax(set)
            Ix[i][j] = set[setmax]
            Ixbt[i][j] = setmax
            if debug:
                print
                print "MI", seq1[i-1],"-", "%d + %d = %d" % (M[i-1][j], -d, M[i-1][j] - d)
                print "II", seq1[i-1],"-", "%d + %d = %d" % (Ix[i-1][j], -e,Ix[i-1][j] - e )
                print "Best Ix", setmax, set[setmax]

            set =  ( M[i][j-1] - d,
                     Iy[i][j-1] - e)
            setmax = argmax((set))
            Iy[i][j] = set[setmax]
            Iybt[i][j] = setmax
            if debug:
                print
                print "M Ix", "-", seq2[j-1], "%d + %d = %d" % (M[i][j-1], -d, M[i][j-1] - d)
                print "IxIx", "-", seq2[j-1], "%d + %d = %d" % (Iy[i][j-1], -e, Iy[i][j-1] - e)
                print "Best Iy", setmax, set[setmax]
                print 
            
    if 1:
        print ""
        print "Final matrix"
        printDPMatrix(seq1, seq2, M)
        printDPMatrix(seq1, seq2, Ix)
        printDPMatrix(seq1, seq2, Iy)
        print ""
        print "Final Traceback"
        printDPMatrix(seq1, seq2, Mbt)
        printDPMatrix(seq1, seq2, Ixbt)
        printDPMatrix(seq1, seq2, Iybt)
        print ""


    if debug:
        print "Score:", M[l_1-1][l_2-1]

    if debug:
        print ""
    aseq1, aseq2 = build_align(seq1, seq2, M, Ix, Iy, Mbt, Ixbt, Iybt)
    return aseq1, aseq2

def printFASTA(queryFile, alignedQuery, subjectFile, alignedSubject, o=sys.stdout):
    print >>o, ">%s" % queryFile
    print >>o, alignedQuery
    print >>o, ">%s" % subjectFile
    print >>o, alignedSubject
    
def main():
    queryFile = None
    subjectFile = None
    openCost = 12
    extendCost = 2
    global debug
    try:
        opts, args = getopt.getopt(sys.argv[1:], "q:s:o:e:d",["query=", "subject=","openCost=","extendCost=","debug="])
    except getopt.GetoptError, what:
	raise RuntimeError, "usage"

    for o, a in opts:

        if o in ('-q', 'query'):
            queryFile = a
        if o in ('-s', 'subject'):
            subjectFile = a
        if o in ('-o', 'openCost'):
            openCost = int(a)
        if o in ('-e', 'extendCost'):
            extendCost = int(a)
        if o in ('-d', 'debug'):
            debug = 1

    if not queryFile or not subjectFile:
        raise RuntimeError, "No query or subject sequence supplied"
 
    alignedQuery, alignedSubject = nw(readFASTA(queryFile), readFASTA(subjectFile),  matrix.blosum50, openCost, extendCost)
    printFASTA(queryFile, alignedQuery, subjectFile, alignedSubject)

def test():
    global debug
    debug = 0
    openCost = 12
    extendCost = 2
    imatrix = matrix.makeIdentity(2, -1)
    # imatrix = matrix.blosum50

    seqpairs = [("".join(map(string.strip, open(sys.argv[1]).readlines()[1:])),
                "".join(map(string.strip, open(sys.argv[2]).readlines()[1:])))]
    print seqpairs

    for seq1, seq2 in seqpairs:
        alignedQuery, alignedSubject = nw(seq1, seq2, imatrix, openCost, extendCost)

        printFASTA("seq1", alignedQuery, "seq2", alignedSubject)

        if string.replace(alignedQuery, "-", "") != seq1:
            print "First sequence was mangled"
        if string.replace(alignedSubject, "-", "") != seq2:
            print "Second sequence was mangled"
        print

if __name__ == '__main__':
    test()
##     main()
