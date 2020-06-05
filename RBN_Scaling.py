import PyBoolNet
import StableMotifs as sm
import multiprocessing as mp


def analyze_rules(rules):
    rules = sm.Format.booleannet2bnet(rules)
    primes = sm.Format.longbnet2primes(rules, remove_constants = True)

    return sm.AttractorRepertoire.from_primes(primes,max_simulate_size=10)



if __name__ == '__main__':
    K=2 # in-degree
    p=sm.RandomBooleanNetworks.get_criticality_p_Kauffman(K)[0] # bias
    N_ensemble=30 # ensemble size
    seed=1000
    TimeoutSeconds = None
    followups = {}
    attbounds = {}

    N_proc = mp.cpu_count()
    if N_proc < N_ensemble:
        N_proc = N_ensemble

    chunksize = int(N_ensemble / N_proc)

    for N in [32,64,128,256,512]: # Number of nodes
        print("Generating ensemble Kauffman RBNs for N =",N,". . .")
        rbn_ensemble_rules=sm.RandomBooleanNetworks.Random_Boolean_Network_Ensemble_Kauffman(N,K,p,N_ensemble,seed=seed,write_Boolean_network=False)
        print("Ensemble generated.")

        # ars = []
        # for rules in rbn_ensemble_rules:
        #     ars.append(analyze_rules(rules))

        with mp.Pool(processes=N_proc) as pool:
            res = pool.map_async(analyze_rules,rbn_ensemble_rules,chunksize)

            try:
                ars = res.get(timeout=TimeoutSeconds)
            except mp.TimeoutError:
                ars = []
                print("Timeout")

        followups[N] = []
        attbounds[N] = []
        for ar in ars:
            attbounds[N].append((ar.fewest_attractors,ar.most_attractors))
            if ar.fewest_attractors != ar.most_attractors:
                followups[N].append(ar)
            print(ar.fewest_attractors,ar.most_attractors)
