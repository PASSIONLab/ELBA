import subprocess as sp
from pathlib import Path
import sys

# Copy this script into the directory containing all your setup runs (obtained via run_elba.py -C)
# and then run it to submit every slurm job.

for p in Path.cwd().iterdir():
    if p.is_dir():
        for j in p.iterdir():
            if j.name.endswith("slurm.sh"):
                command = ["sbatch", "--chdir", str(p), str(j)]
                sys.stdout.write(" ".join(command) + "\n")       
                proc = sp.Popen(command)
                proc.wait()

