#!/usr/bin/env python

import sys

def main(file1, file2):
    A = {int(line.rstrip()) for line in open(file1)}
    B = {int(line.rstrip()) for line in open(file2)}
    sys.stdout.write("A='{}',B='{}'\n".format(file1, file2))
    sys.stdout.write("|A|={},|B|={}\n".format(len(A), len(B)))
    sys.stdout.write("|A union B| = {}\n".format(len(A.union(B))))
    sys.stdout.write("|A intersect B| = {}\n".format(len(A.intersection(B))))
    sys.stdout.write("|A \ B| = {}\n".format(len(A.difference(B))))
    sys.stdout.write("|B \ A| = {}\n".format(len(B.difference(A))))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("usage: {} <file1> <file2>\n".format(sys.argv[0]))
        sys.stderr.flush()
        sys.exit(-1)
    main(sys.argv[1], sys.argv[2])
    sys.exit(0)
