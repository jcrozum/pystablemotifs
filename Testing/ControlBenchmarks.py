import PyStableMotifs as sm
import PyBoolNet as pbn

import atexit
import os
import shutil
import subprocess
import sys
import tempfile

from timeit import default_timer

# The old Java version of StableMotifs relies on Nashorn, which is removed from recent
# Java versions. As a result, we need to specify a path to an old version of Java.
java_path = "C:\\Users\\jcroz\\Documents\\Work\\StableMotifsJavaVersion\\OLDJAVA\\bin\\java.exe"
jar_path = "C:\\Users\\jcroz\\Documents\\Work\\StableMotifsJavaVersion\\dist\\PyStableMotifs.jar"


def JavaMotifTimes(filepath,verbose=False):
    # prepare temporary working space
    wd = tempfile.mkdtemp(prefix="StableMotifs-")
    def cleanup():
        shutil.rmtree(wd)
    atexit.register(cleanup)

    shutil.copy(filepath, wd)
    base = os.path.basename(filepath)

    argv = [java_path,'-jar',jar_path, base]
    proc = subprocess.Popen(argv, cwd=wd,stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    result = proc.communicate()
    if verbose:
        print(result[0])

    atime = -1.0
    ctime = -1.0
    try:
        for line in result[0].decode('utf8').split('\r\n'):
            if line.startswith("Total time for finding quasi-attractors:"):
                atime = float(line.split(':')[1].strip(" s"))
            elif line.startswith("Total time for finding stable motif control sets:"):
                ctime = float(line.split(':')[1].strip(" s"))
    except:
        print("Error during time extraction.")
        pass
    return (atime,ctime)

def CompareTimes(filepath,display=True):
    p=sm.Format.import_primes(filepath,remove_constants=False)
    a2013,c2013 = JavaMotifTimes(filepath)

    start=default_timer()
    pbn.PrimeImplicants.percolate_and_remove_constants(p) # modifies p in-place
    ar = sm.AttractorRepertoire.from_primes(p,max_simulate_size=0) # max_simulate_size = 0 ensures equivalence with Java method
    end=default_timer()
    a2021 = end-start

    start=default_timer()
    for att in ar.attractors:
        fixed = att.logically_fixed_nodes
        reprogram = ar.succession_diagram.reprogram_to_trap_spaces(fixed,target_method='history',driver_method='internal')
    end=default_timer()
    c2021 = end-start
    if display:
        print("Trap Space ID (old,new):  ",a2013,", ",round(a2021,9),sep="")
        print("      Control (old,new):  ",c2013,", ",round(c2021,9),sep="")
        print("     Combined (old,new):  ",round(a2013+c2013,9),", ",round(a2021+c2021,9),sep="")
    return (a2013,c2013,a2021,c2021)

Ns = [10,20,30,40,50]
K=2
bias=0.5
N_ensemble=10
seed = 0
rbn_fnames = []
for N in Ns:
    rbns=sm.RandomBooleanNetworks.Random_Boolean_Network_Ensemble_Kauffman(N,K,bias,N_ensemble,seed=seed)
    for i,rbn in enumerate(rbns):
        rbn_bnet = sm.Format.booleannet2bnet(rbn)
        p = sm.Format.longbnet2primes(rbn_bnet,remove_constants=False)
        fname = "../models/Control Benchmarks/RBNs_seed="+str(seed)+"/"+str(N)+"_"+str(i)+".txt"
        rbn_fnames.append(fname)
        f = open(fname,"w")
        f.write(sm.Format.primes2booleannet(p,header="#BOOLEAN RULES\n"))
        f.close()

f = open("../models/Control Benchmarks/RBNs_seed="+str(seed)+"/timings.csv","w")
f.write("model,a2013,c2013,a2021,c2021\n")

for fname in rbn_fnames:
    print()
    print(fname)
    times = CompareTimes(fname)
    f.write(','.join([fname]+list(map(str,times)))+'\n')

f.close()

quit()
print()
print("SmallTest:")
filepath = "../models/Control Benchmarks/SmallTest.txt"
CompareTimes(filepath)
print()
print("PhaseSwitch:")
filepath = "../models/Control Benchmarks/PhaseSwitch.txt"
CompareTimes(filepath)
print()
print("ThNetwork:")
filepath = "../models/Control Benchmarks/ThNetwork.txt"
CompareTimes(filepath)
print()
print("TLGL_Small:")
filepath = "../models/Control Benchmarks/TLGL_Small.txt"
CompareTimes(filepath)
print()
print("TLGL_Large_Fixed_Inputs:")
filepath = "../models/Control Benchmarks/TLGL_Large_Fixed_Inputs.txt"
CompareTimes(filepath)
print()
