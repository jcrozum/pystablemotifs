import PyBoolNet
import PyStableMotifs as sm
from timeit import default_timer

print("Generating ensemble of 100 Kauffman RBN with K=2 and N=500 at criticality (p=0.5) . . .")
N=500
K=2
p=sm.RandomBooleanNetworks.get_criticality_p_Kauffman(K)[0]
N_ensemble=100
seed=1000
start=default_timer()
rbn_ensemble_rules_K2_pcrit=sm.RandomBooleanNetworks.Random_Boolean_Network_Ensemble_Kauffman(N,K,p,N_ensemble,seed=seed,write_Boolean_network=True)
end=default_timer()
print("Time (s) generating ensemble:",end-start)

print("\nGenerating ensemble of 100 Kauffman RBN with K=2 and N=500 not at criticality (p=0.4) . . .")
N=500
K=2
p=0.4
N_ensemble=100
seed=1000
start=default_timer()
rbn_ensemble_rules_K2_p04=sm.RandomBooleanNetworks.Random_Boolean_Network_Ensemble_Kauffman(N,K,p,N_ensemble,seed=seed,write_Boolean_network=True)
end=default_timer()
print("Time (s) generating ensemble:",end-start)

print("\nGenerating ensemble of 100 Kauffman RBN with K=3 and N=500 at criticality . . .")
N=500
K=3
p=sm.RandomBooleanNetworks.get_criticality_p_Kauffman(K)[0]
N_ensemble=100
seed=1000
start=default_timer()
rbn_ensemble_rules_K3_pcrit=sm.RandomBooleanNetworks.Random_Boolean_Network_Ensemble_Kauffman(N,K,p,N_ensemble,seed=seed,write_Boolean_network=True)
end=default_timer()
print("Time (s) generating ensemble:",end-start)

print("\nGenerating ensemble of 100 Kauffman RBN with K=3 and N=500 not at criticality . . .")
N=500
K=3
p=0.5
N_ensemble=100
seed=1000
start=default_timer()
rbn_ensemble_rules_K3_p05=sm.RandomBooleanNetworks.Random_Boolean_Network_Ensemble_Kauffman(N,K,p,N_ensemble,seed=seed,write_Boolean_network=True)
end=default_timer()
print("Time (s) generating ensemble:",end-start)

print("\nPrinting out original rules and reduced prime implicant rules for a RBN")
rules=rbn_ensemble_rules_K3_pcrit[0]
print("\nOriginal rules:")
print(rules)
rules = sm.Format.booleannet2bnet(rules)
primes = sm.Format.longbnet2primes(rules)
PyBoolNet.PrimeImplicants._percolation(primes,True)
print("\nReduced prime implicant rules:")
sm.Format.pretty_print_prime_rules(primes)
