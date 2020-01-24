import logging
import io

logger = logging.getLogger('dp_matrix')
logger.setLevel(logging.INFO)

def printDPMatrix(t, s1, s2, m):
    s = io.StringIO()
    s.write("DP matrix %s\n" % t)
    s.write("%4s %4s " % (" ", " "))
    for i in range(1, len(m)):
        s.write("%4s " % s1[i-1])
    s.write("\n")

    for i in range(len(m[0])):
        if i == 0:
            s.write("%4s " % " ")
        else:
            s.write("%4s " % s2[i-1])

        for j in range(len(m)):
            s.write("%4d " % m[j][i])
        s.write("\n")
    logger.debug(s.getvalue())

