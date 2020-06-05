import PyBoolNet
import StableMotifs as sm


K=2 # in-degree
p=sm.RandomBooleanNetworks.get_criticality_p_Kauffman(K)[0] # bias
N_ensemble=5 # ensemble size
seed=1000
followups = {}
attbounds = {}
for N in [32,64,128,256,512]: # Number of nodes
    print("Generating ensemble Kauffman RBNs for N =",N,". . .")
    rbn_ensemble_rules=sm.RandomBooleanNetworks.Random_Boolean_Network_Ensemble_Kauffman(N,K,p,N_ensemble,seed=seed,write_Boolean_network=False)
    print("Ensemble generated.")
    followups[N] = []
    attbounds[N] = []
    for rules in rbn_ensemble_rules:
        rules = sm.Format.booleannet2bnet(rules)
        primes = sm.Format.longbnet2primes(rules, remove_constants = True)

        ar = sm.AttractorRepertoire.from_primes(primes,max_simulate_size=10)
        attbounds[N].append((ar.fewest_attractors,ar.most_attractors))

        if ar.fewest_attractors != ar.most_attractors:
            followups[N].append(ar)
        print(ar.fewest_attractors,ar.most_attractors)
