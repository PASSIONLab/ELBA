from decimal import *
getcontext().prec = 128

from math import factorial as fact
import sys

def binom(n, k):
    return fact(n)/(fact(k) * fact(n-k))

def P(m, d, e, k):
    p = Decimal(1-e)**k

    coef = binom(d,m)
    p1 = p ** m
    p2 = (1-p) **(d-m)
    return Decimal(coef) * p1 * p2

def main(d, e, k, minprob):

    mysum = 0
    m = 2
    while (mysum < minprob):
        prob = P(m, d, e, k)
        mysum += prob
        m += 1

    l = m - 1

    mysum = 0
    m = d
    while (mysum < minprob):
        prob = P(m, d, e, k)
        mysum += prob
        m -= 1

    u = m + 1

    print("lower = {}, upper = {}".format(l, u))

if __name__ == "__main__":
    main(int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]))
