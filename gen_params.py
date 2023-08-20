import itertools

L = [15]
U = [35,55]
k = [51]
A = [1]
B = [1]
f = [0.99]
s = [2000]
x = [15, 30]
c = [0.2, 0.5, 0.65]

parameters = list(itertools.product(L, U, k, A, B, f, s, x, c))

with open("parameters.txt", "w") as f:
    f.write("L,U,k,A,B,f,s,x,c,reads_fname,ref_fname,outdir\n")
    for i, p in enumerate(parameters):
        outdir = "out_L{}_U{}_k{}_A{}_B{}_f{}_s{}_x{}_c{}".format(p[0], p[1], p[2], p[3], p[4], int(p[5]*100), p[6], p[7], int(p[8]*100))
        f.write(",".join([str(item) for item in p]) + ",reads.fa,ref.fa,{}\n".format(outdir))
    
