import sys

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
            sequence.append(line.strip())

    return "".join(sequence)


def writeFASTA(queryFile, alignedQuery, subjectFile, alignedSubject, o=sys.stdout):
    print(">%s" % queryFile, file=o)
    print(alignedQuery, file=o)
    print(">%s" % subjectFile, file=o)
    print(alignedSubject, file=o)
