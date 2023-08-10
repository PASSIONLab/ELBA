#!/usr/bin/env python

import sys
import subprocess as sp
import shutil
from pathlib import Path

command_id = 1

def logmsg(msg):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()

def run_command(cmdlist):
    global command_id
    step_id = 1
    logmsg("\n({}.*) executing: '{}'\n".format(command_id, " ".join(cmdlist)))
    proc = sp.Popen(cmdlist, stdout=sp.PIPE)
    while True:
        output = proc.stdout.readline().decode()
        if output == "" and proc.poll() is not None:
            break
        if output:
            logmsg("\t({}.{}). {}".format(command_id, step_id, output.rstrip()))
            step_id += 1
    command_id += 1

def get_experiments(params_fname):
    experiments = {}
    with open(params_fname, "r") as f:
        header = next(f).rstrip()
        assert header == "L,U,k,A,B,f,s,x,c,reads_fname,ref_fname,outdir"
        for line in f.readlines():
            L,U,k,A,B,fr,s,x,c,fasta_fname,ref_fname,outdir = (t(tok) for t,tok in zip((int,int,int,int,int,float,int,int,float,str,str,str), line.rstrip().split(",")))
            fasta_path = Path(fasta_fname).resolve()
            ref_path = Path(ref_fname).resolve()
            outdir_path = Path(outdir).resolve()
            assert fasta_path.is_file()
            assert ref_path.is_file()
            outdir_path.mkdir(exist_ok=False)
            if not L in experiments: experiments[L] = {}
            if not U in experiments[L]: experiments[L][U] = {}
            if not k in experiments[L][U]: experiments[L][U][k] = []
            experiments[L][U][k].append((A, B, fr, s, x, c, fasta_path, ref_path, outdir_path))
    return experiments

def create_slurm_script(experiment, exe_path):
    exe_path2 = experiment[-1].joinpath("elba").resolve()
    assert not exe_path2.exists()
    shutil.copyfile(str(exe_path), str(exe_path2))
    A, B, fr, s, x, c, fasta_path, ref_path, outdir_path = experiment
    script_path = outdir_path.joinpath("slurm.sh")
    assert not script_path.exists()
    with open(str(script_path), "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write("#SBATCH -q regular\n")
        f.write("#SBATCH -t 5\n")
        f.write("#SBATCH -C cpu\n")
        f.write("#SBATCH -N 1\n\n")
        f.write("chmod +x {}\n".format(str(exe_path2)))
        cmdlist = ["srun", "-n", "121", "-N", "1", "-c", "2", "--cpu_bind=cores", str(exe_path2)]
        cmdlist += ["-A", str(A), "-B", str(B), "-G", str(B), "-x", str(x), "-f", str(fr), "-s", str(s), "-c", str(c), str(fasta_path)]
        f.write("{}\n".format(" ".join(cmdlist)))

def compile_elba(L, U, k):
    elba_path = Path("/pscratch/sd/g/gabeh98/August2023/ELBA").resolve()
    assert elba_path.is_dir()

    proc = sp.Popen(["make", "-C", str(elba_path), "clean"], stdout=sp.PIPE)
    proc.wait()

    proc = sp.Popen(["make", "-C", str(elba_path), "L={}".format(L), "U={}".format(U), "K={}".format(k), "-j", "32"], sp.PIPE)
    proc.wait()

    exe_path = elba_path.joinpath("elba")
    assert exe_path.is_file()
    return exe_path

def main(argc, argv):
    experiments = get_experiments(argv[1])

    for L in experiments:
        for U in experiments[L]:
            for k in experiments[L][U]:
                exe_path = compile_elba(L, U, k)
                for exp in experiments[L][U][k]:
                    create_slurm_script(exp, exe_path)

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
