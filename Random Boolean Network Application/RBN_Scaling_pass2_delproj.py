import pyboolnet
from pyboolnet.external.bnet2primes import bnet_text2primes
import pystablemotifs as sm
import multiprocessing as mp
import time
from timeit import default_timer
import pandas as pd

from pathlib import Path

max_simulate_size=10

# the name of the file that the summary will be printed to
summary_output_fname = "network_summaries_delproj.csv"

# The task that will be passed to the workers in the parallelization loop:
# Given (i,primes), compute the attractor repertoire object. We keep track of i
# so that we can do follow-up on networks that timeout or have other difficulties.
def analyze_rules(x):
    i = x[0]
    primes,constants = sm.reduction.deletion_reduction(x[1])
    return i,sm.AttractorRepertoire.from_primes(primes,max_simulate_size=max_simulate_size)

# Function to turn ensemble into a dictionary of pyboolnet primes with index keys
def primify_rbn(rbn_ensemble_rules):
    prime_dict = {}
    for i,x in enumerate(rbn_ensemble_rules):
        rules = sm.format.booleannet2bnet(x)
        primes = bnet_text2primes(rules)
        pyboolnet.prime_implicants.percolation(primes,True)
        primes, constants = sm.reduction.remove_outdag(primes)
        prime_dict[i] = primes
    return prime_dict

def prune_primes(rbn_primes,prev_data,N):
    df = prev_data.loc[prev_data['N'] == N]

    new_primes = {}
    for i,p in rbn_primes.items():
        row = df[df['index']==i]
        if len(row) == 1:
            L=int(row['lbound'].iloc[0])
            U=int(row['ubound'].iloc[0])
            if L < U:
                new_primes[i] = p
        else:
            new_primes[i] = p

    return new_primes

# Do not remove this if statement; multiprocessing needs it
if __name__ == '__main__':
    # How many cores do we want to use?
    N_proc = mp.cpu_count()

    # Ensemble parameters
    K=2 # in-degree
    p_bias=sm.random_boolean_networks.get_criticality_p_Kauffman(K)[0] # bias
    N_ensemble=90 # ensemble size
    seed=1000

    # Timeout parameters
    # After SoftTimeCapSeconds seconds, we give up if we go more than
    # GetterTimeoutSeconds without obtaining a new result
    SoftTimeCapSeconds = 3600
    GetterTimeoutSeconds = 30

    pname = "delproj_p="+str(p_bias)+"_K="+str(K)+"_seed="+str(seed)+"_simsize="+\
        str(max_simulate_size)+"_timeouts="+str(GetterTimeoutSeconds)+","+\
        str(SoftTimeCapSeconds)+"/"
    Path("./"+pname).mkdir(exist_ok=True)
    fname_in = pname + "network_summaries_prev.csv"
    fname_out = pname + summary_output_fname

    prev_data = pd.read_csv(fname_in,sep=',',header=0)

    # clear file contents and write header:
    with open(fname_out,"w") as file:
        print("N,index,lbound,ubound",file=file)

    chunksize = 1 # How big should the list be that we pass to each worker?

    attbounds = {} # these are the results we're interested in

    # Main loop
    for N in [181]: # Number of nodes (before reduction)
        print("Generating ensemble Kauffman RBNs for N =",N,". . .")
        rbn_ensemble_rules=sm.random_boolean_networks.Random_Boolean_Network_Ensemble_Kauffman(N,K,p_bias,N_ensemble,seed=seed,write_Boolean_network=False)
        rbn_primes = primify_rbn(rbn_ensemble_rules)
        rbn_primes = prune_primes(rbn_primes,prev_data,N)
        print("Ensemble generated.")
        ars={} # will store calculated attractor repertoires

        start=default_timer() # Timer used for testing timeouts

        # The parallel loop
        with mp.Pool(processes=N_proc) as pool:
            res = pool.imap_unordered(analyze_rules,rbn_primes.items(),chunksize)
            while 1:
                try: # we calculated in time
                    i,r = res.next(timeout=GetterTimeoutSeconds)
                    ars[i]=r
                    G=sm.export.networkx_succession_diagram_reduced_network_based(r)
                    sm.export.save_to_graphml(G,pname+"N="+str(N)+",i="+str(i))
                except mp.TimeoutError: # we are too slow . . .
                    if default_timer() - start > SoftTimeCapSeconds: # try again?
                        break
                except StopIteration: # finished calculating everything
                    break

        with open(fname_out, "a") as file:
            for ind,ar in ars.items():
                lbound = ar.fewest_attractors
                # complex attractors in deletion projection can be "fake"
                for att in ar.attractors:
                    if att.n_unfixed > 0 and lbound > 1:
                        lbound -= 1
                ubound = ar.most_attractors
                print(N,ind,lbound,ubound,file=file,sep=",")
