import PyBoolNet
import StableMotifs as sm

print("Generating ensemble Kauffman RBNs . . .")
N=200 # Number of nodes
K=2 # in-degree
p=sm.RandomBooleanNetworks.get_criticality_p_Kauffman(K)[0] # bias
N_ensemble=30 # ensemble size
seed=1000
rbn_ensemble_rules=sm.RandomBooleanNetworks.Random_Boolean_Network_Ensemble_Kauffman(N,K,p,N_ensemble,seed=seed,write_Boolean_network=True)
print("Ensemble generated.")

followups = []
for rules in rbn_ensemble_rules:
    rules = sm.Format.booleannet2bnet(rules)
    primes = sm.Format.longbnet2primes(rules)
    ar = sm.AttractorRepertoire.from_primes(primes,max_simulate_size=10)
    att_bounds = (ar.fewest_attractors,ar.most_attractors)
    if att_bounds[0] != att_bounds[1]:
        followups.append(ar)
    print(att_bounds)
