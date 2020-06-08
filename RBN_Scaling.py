import PyBoolNet
import StableMotifs as sm
import multiprocessing as mp


def analyze_rules(primes):
    return sm.AttractorRepertoire.from_primes(primes,max_simulate_size=10)

def primify_rbn(rbn_ensemble_rules):
    prime_list = []
    for x in rbn_ensemble_rules:
        rules = sm.Format.booleannet2bnet(x)
        prime_list.append(sm.Format.longbnet2primes(rules, remove_constants = True))
    return prime_list

if __name__ == '__main__':
    K=2 # in-degree
    p=sm.RandomBooleanNetworks.get_criticality_p_Kauffman(K)[0] # bias
    N_ensemble=30 # ensemble size
    seed=1000
    TimeoutSeconds = 5

    followups_prime = {}
    followups_ar = {}
    attbounds = {}

    N_proc = mp.cpu_count()
    if N_proc < N_ensemble:
        N_proc = N_ensemble

    chunksize = int(N_ensemble / N_proc)

    for N in [32,64,128,256,512]: # Number of nodes
        print("Generating ensemble Kauffman RBNs for N =",N,". . .")
        rbn_ensemble_rules=sm.RandomBooleanNetworks.Random_Boolean_Network_Ensemble_Kauffman(N,K,p,N_ensemble,seed=seed,write_Boolean_network=False)
        rbn_primes = primify_rbn(rbn_ensemble_rules)
        print("Ensemble generated.")

        # ars = []
        # for rules in rbn_ensemble_rules:
        #     ars.append(analyze_rules(rules))
        ars=[]
        with mp.Pool(processes=N_proc) as pool:
            res = pool.imap_unordered(analyze_rules,rbn_primes,chunksize)

            while 1:
                try:
                    ars.append(res.next(timeout=TimeoutSeconds))
                except mp.TimeoutError:
                    print("Timeout")
                except StopIteration:
                    break

        followups_prime[N] = []
        for x in rbn_primes:
            found = False
            for ar in ars:
                if PyBoolNet.PrimeImplicants.are_equal(x,ar.succession_diagram.unreduced_primes):
                    found = True
                    break
            if not found:
                followups_prime[N].append(x)
                print()
                print("Timed-out Rules:")
                sm.Format.pretty_print_prime_rules(followups_prime[N][-1])


        followups_ar[N] = []
        attbounds[N] = []
        for ar in ars:
            attbounds[N].append((ar.fewest_attractors,ar.most_attractors))
            if ar.fewest_attractors != ar.most_attractors:
                followups_ar[N].append(ar)
            print(ar.fewest_attractors,ar.most_attractors)
