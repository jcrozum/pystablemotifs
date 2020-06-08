import PyBoolNet
import StableMotifs as sm
import multiprocessing as mp
import time
from timeit import default_timer

def analyze_rules(x):
    i = x[0]
    primes = x[1]
    return i,sm.AttractorRepertoire.from_primes(primes,max_simulate_size=10)

def primify_rbn(rbn_ensemble_rules):
    prime_dict = {}
    for i,x in enumerate(rbn_ensemble_rules):
        rules = sm.Format.booleannet2bnet(x)
        prime_dict[i] = sm.Format.longbnet2primes(rules, remove_constants = True)
    return prime_dict

if __name__ == '__main__':
    K=2 # in-degree
    p=sm.RandomBooleanNetworks.get_criticality_p_Kauffman(K)[0] # bias
    N_ensemble=20 # ensemble size
    seed=1000

    # After SoftTimeCapSeconds seconds, we give up if we go more than
    # GetterTimeoutSeconds without obtaining a new result
    SoftTimeCapSeconds = 30
    GetterTimeoutSeconds = 5

    followups_prime = {}
    followups_ar = {}
    attbounds = {}

    N_proc = mp.cpu_count()
    if N_proc < N_ensemble:
        N_proc = N_ensemble

    chunksize = 1#int(N_ensemble / N_proc)

    for N in [32,64]:#,128,256,512]: # Number of nodes
        start=default_timer()
        print("Generating ensemble Kauffman RBNs for N =",N,". . .")
        rbn_ensemble_rules=sm.RandomBooleanNetworks.Random_Boolean_Network_Ensemble_Kauffman(N,K,p,N_ensemble,seed=seed,write_Boolean_network=False)
        rbn_primes = primify_rbn(rbn_ensemble_rules)
        print("Ensemble generated.")

        ars={}
        with mp.Pool(processes=N_proc) as pool:
            res = pool.imap_unordered(analyze_rules,rbn_primes.items(),chunksize)
            while 1:
                try:
                    i,r = res.next(timeout=GetterTimeoutSeconds)
                    ars[i]=r
                except mp.TimeoutError:
                    if default_timer() - start > SoftTimeCapSeconds:
                        print("Timeout")
                        break
                except StopIteration:
                    print("Done")
                    break

        followups_prime[N] = [p for i,p in rbn_primes.items() if i not in ars]
        for p in followups_prime[N]:
            print()
            print("Timeout Primes")
            sm.Format.pretty_print_prime_rules(p)
            print()

        followups_ar[N] = []
        attbounds[N] = []
        for ar in ars.values():
            attbounds[N].append((ar.fewest_attractors,ar.most_attractors))
            if ar.fewest_attractors != ar.most_attractors:
                followups_ar[N].append(ar)
            print(ar.fewest_attractors,ar.most_attractors)

    for k in attbounds:
        L = attbounds[k]
        print(k,sum([x[0] for x in L]) / len(L),sum([x[1] for x in L]) / len(L))
