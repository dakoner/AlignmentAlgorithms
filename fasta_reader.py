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
