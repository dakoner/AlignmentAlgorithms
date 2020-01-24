## a miniature implementation of Needleman Wunsch algorithm using three
## matrices M, Ix, Iy, and a backtrace matrix

import sys
import argparse
import logging
import io
import matrix
from dp_matrix import printDPMatrix
from fasta_reader import readFASTA

logging.basicConfig()
logger = logging.getLogger('mininwcostfree')
logger.setLevel(logging.INFO)


def argmax(array):
    return array.index(max(array))

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

    reverse = lambda s: ''.join([s[i] for i in range(len(s)-1, -1, -1)])
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
        logger.debug("%d %d" % (i,j))
        ## if current matrix is match, emit the relevant sequence
        if id(matrix) == id(Mbt):
            logger.debug("Match matrix")
            ## if either sequence has reached its left edge, gap it rather than using a negative index
            ## why would this be the case, ever?
            if i < 1:
                logger.debug("emit s1 -")
                s1 = "-"
            else:
                logger.debug("emit s1: %s", sequence1[i-1])
                s1 = sequence1[i-1]
            if j < 1:
                logger.debug("emit s2 -")
                s2 = "-"
            else:
                logger.debug("emit s2: %s", sequence2[j-1])
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
                logger.debug("Transition to Ixbt")
                matrix = Ixbt
            ##M->Iy
            elif matrix[i][j] == 2:
                if i > 0: i -= 1
                if j > 0: j -= 1
                logger.debug("Transition to Iybt")
                matrix = Iybt

        ## if current matrix is Ix (insert in second sequence) then emit a char from the first sequence and gap the second sequence
        elif id(matrix) == id(Ixbt):
            logger.debug("Ix matrix")
            s1 = sequence1[i-1]
            s2 = "-"
            logger.debug("emit s1 %s emit s2 -" % sequence1[i-1])
            if matrix[i][j] == 0:
                logger.debug("Transition to Mbt")
                matrix = Mbt
            i -= 1

        ## if current matrix is Iy (insert in first sequence) then gap the first sequence and emit a char from the second sequence
        elif id(matrix) == id(Iybt):
            logger.debug("Iy matrix")
            s1 = "-"
            s2 = sequence2[j-1]
            logger.debug("emit s1 - emit s2 %s" % sequence2[j-1])
            if matrix[i][j] == 0:
                logger.debug("Transition to Mbt")
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


    logger.debug("Initial")
    printDPMatrix(seq1, seq2, M)
    printDPMatrix(seq1, seq2, Ix)
    printDPMatrix(seq1, seq2, Iy)

    for j in range(1, l_2):
        for i in range(1, l_1):
            logger.debug( "Cell %d %d" % (i,j))
            match = matrix[seq1[i-1],seq2[j-1]]

            set =  (  M[i-1][j-1] + match,
                     Ix[i-1][j-1] + match,
                     Iy[i-1][j-1] + match)
            setmax = argmax((set))
            M[i][j] = set[setmax]
            Mbt[i][j] = setmax

            logger.debug("M M %s %s %d + %d = %d" % (seq1[i-1],seq2[j-1], M[i-1][j-1], match, M[i-1][j-1] + match))
            logger.debug("IxM %s - %d + %d = %d" % (seq1[i-1], Ix[i-1][j-1], match, Ix[i-1][j-1] + match))
            logger.debug("IyM - %s %d + %d = %d" % (seq2[j-1], Iy[i-1][j-1], match, Iy[i-1][j-1] + match))
            logger.debug("Best M %s %s", setmax, set[setmax])
                
            set =  ( M[i-1][j] - d,
                    Ix[i-1][j] - e)
            setmax = argmax(set)
            Ix[i][j] = set[setmax]
            Ixbt[i][j] = setmax
            logger.debug("MI %s - %d + %d = %d" % (seq1[i-1], M[i-1][j], -d, M[i-1][j] - d))
            logger.debug("II %s - %d + %d = %d" % (seq1[i-1],Ix[i-1][j], -e,Ix[i-1][j] - e))
            logger.debug("Best Ix %s %s" % (setmax, set[setmax]))

            set =  ( M[i][j-1] - d,
                     Iy[i][j-1] - e)
            setmax = argmax((set))
            Iy[i][j] = set[setmax]
            Iybt[i][j] = setmax
            logger.debug("M Ix - %s %d + %d = %d" % (seq2[j-1], M[i][j-1], -d, M[i][j-1] - d))
            logger.debug("IxIx - %s %d + %d = %d" % (seq2[j-1], Iy[i][j-1], -e, Iy[i][j-1] - e))
            logger.debug("Best Iy %s %s" % (setmax, set[setmax]))
            
    if 1:
        logger.debug("Final matrix")
        printDPMatrix(seq1, seq2, M)
        printDPMatrix(seq1, seq2, Ix)
        printDPMatrix(seq1, seq2, Iy)
        logger.debug("Final Traceback")
        printDPMatrix(seq1, seq2, Mbt)
        printDPMatrix(seq1, seq2, Ixbt)
        printDPMatrix(seq1, seq2, Iybt)


    logger.debug("Score %d:", M[l_1-1][l_2-1])

    aseq1, aseq2 = build_align(seq1, seq2, M, Ix, Iy, Mbt, Ixbt, Iybt)
    return aseq1, aseq2

def printFASTA(queryFile, alignedQuery, subjectFile, alignedSubject, o=sys.stdout):
    print(">%s" % queryFile, file=o)
    print(alignedQuery, file=o)
    print(">%s" % subjectFile, file=o)
    print(alignedSubject, file=o)
    
def main():
    queryFile = None
    subjectFile = None
    parser = argparse.ArgumentParser("Parse alginment arguments")
    parser.add_argument('queryFile', type=str, default=None)
    parser.add_argument('subjectFile', type=str, default=None)
    parser.add_argument('-o', '--openCost', type=float, default=11)
    parser.add_argument('-e', '--extendCost', type=float, default=1)
    parser.add_argument('-d', dest='debug', action='store_true')
    args = parser.parse_args()
    if args.debug:
      logger.setLevel(logging.DEBUG)
    if not args.queryFile or not args.subjectFile:
        raise RuntimeError("No query or subject sequence supplied")
 
    alignedQuery, alignedSubject = nw(readFASTA(args.queryFile), readFASTA(args.subjectFile),  matrix.blosum62, args.openCost, args.extendCost)
    printFASTA(queryFile, alignedQuery, subjectFile, alignedSubject)

if __name__ == '__main__':
    main()
