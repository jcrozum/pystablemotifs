import PyBoolNet
import StableMotifs as sm
import multiprocessing as mp
import time
from timeit import default_timer

# The task that will be passed to the workers in the parallelization loop:
# Given (i,primes), compute the attractor repertoire object. We keep track of i
# so that we can do follow-up on networks that timeout or have other difficulties.
def analyze_rules(x):
    i = x[0]
    primes = x[1]
    return i,sm.AttractorRepertoire.from_primes(primes,max_simulate_size=10)

# Function to turn ensemble into a dictionary of PyBoolNet primes with index keys
def primify_rbn(rbn_ensemble_rules):
    prime_dict = {}
    for i,x in enumerate(rbn_ensemble_rules):
        rules = sm.Format.booleannet2bnet(x)
        prime_dict[i] = sm.Format.longbnet2primes(rules, remove_constants = True)
    return prime_dict

# Do not remove this if statement; multiprocessing needs it
if __name__ == '__main__':

    # Ensemble parameters
    K=2 # in-degree
    p_bias=sm.RandomBooleanNetworks.get_criticality_p_Kauffman(K)[0] # bias
    N_ensemble=20 # ensemble size
    seed=1000

    # Timeout parameters
    # After SoftTimeCapSeconds seconds, we give up if we go more than
    # GetterTimeoutSeconds without obtaining a new result
    SoftTimeCapSeconds = 30
    GetterTimeoutSeconds = 5

    followups_prime = {} # The primes that timed out
    followups_ar = {} # these are attractor repertoire objects that have
                      #different upper and lower bounds for the attractor number
    attbounds = {} # these are the results we're interested in

    # How many cores do we want to use?
    N_proc = mp.cpu_count()
    if N_proc < N_ensemble:
        N_proc = N_ensemble

    chunksize = 1 # How big should the list be that we pass to each worker?

    for N in [32,64,128,256,512]: # Number of nodes (before reduction)
        print("Generating ensemble Kauffman RBNs for N =",N,". . .")
        rbn_ensemble_rules=sm.RandomBooleanNetworks.Random_Boolean_Network_Ensemble_Kauffman(N,K,p_bias,N_ensemble,seed=seed,write_Boolean_network=False)
        rbn_primes = primify_rbn(rbn_ensemble_rules)
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
                except mp.TimeoutError: # we are too slow . . .
                    if default_timer() - start > SoftTimeCapSeconds: # try again?
                        print("Timeout")
                        break
                except StopIteration: # finished calculating everything
                    print("Done")
                    break

        # Collect and display the primes that timed out
        followups_prime[N] = [p for i,p in rbn_primes.items() if i not in ars]
        for p in followups_prime[N]:
            print()
            print("Timeout Primes")
            sm.Format.pretty_print_prime_rules(p)
            print()

        # Grab the bounds on attractor number, and keep track of when they differ
        followups_ar[N] = []
        attbounds[N] = []
        for ar in ars.values():
            attbounds[N].append((ar.fewest_attractors,ar.most_attractors))
            if ar.fewest_attractors != ar.most_attractors:
                followups_ar[N].append(ar)
            print(ar.fewest_attractors,ar.most_attractors)

    # Display average upper and lower bounds on attractor number
    for k in attbounds:
        L = attbounds[k]
        print(k,sum([x[0] for x in L]) / len(L),sum([x[1] for x in L]) / len(L))
