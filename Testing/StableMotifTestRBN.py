import pyboolnet.prime_implicants
from pyboolnet.external.bnet2primes import bnet_text2primes
import pystablemotifs as sm
from timeit import default_timer

print("Generating ensemble of 100 Kauffman RBN with K=2 and N=500 at criticality (p=0.5) . . .")
N=500
K=2
p=sm.random_boolean_networks.get_criticality_p_Kauffman(K)[0]
N_ensemble=100
seed=1000
start=default_timer()
rbn_ensemble_rules_K2_pcrit=sm.random_boolean_networks.random_boolean_network_ensemble_kauffman(N,K,p,N_ensemble,seed=seed,write_boolean_network=True)
end=default_timer()
print("Time (s) generating ensemble:",end-start)

print("\nGenerating ensemble of 100 Kauffman RBN with K=2 and N=500 not at criticality (p=0.4) . . .")
N=500
K=2
p=0.4
N_ensemble=100
seed=1000
start=default_timer()
rbn_ensemble_rules_K2_p04=sm.random_boolean_networks.random_boolean_network_ensemble_kauffman(N,K,p,N_ensemble,seed=seed,write_boolean_network=True)
end=default_timer()
print("Time (s) generating ensemble:",end-start)

print("\nGenerating ensemble of 100 Kauffman RBN with K=3 and N=500 at criticality . . .")
N=500
K=3
p=sm.random_boolean_networks.get_criticality_p_Kauffman(K)[0]
N_ensemble=100
seed=1000
start=default_timer()
rbn_ensemble_rules_K3_pcrit=sm.random_boolean_networks.random_boolean_network_ensemble_kauffman(N,K,p,N_ensemble,seed=seed,write_boolean_network=True)
end=default_timer()
print("Time (s) generating ensemble:",end-start)

print("\nGenerating ensemble of 100 Kauffman RBN with K=3 and N=500 not at criticality . . .")
N=500
K=3
p=0.5
N_ensemble=100
seed=1000
start=default_timer()
rbn_ensemble_rules_K3_p05=sm.random_boolean_networks.random_boolean_network_ensemble_kauffman(N,K,p,N_ensemble,seed=seed,write_boolean_network=True)
end=default_timer()
print("Time (s) generating ensemble:",end-start)

print("\nPrinting out original rules and reduced prime implicant rules for a RBN")
rules=rbn_ensemble_rules_K3_pcrit[0]
print("\nOriginal rules:")
print(rules)
rules = sm.format.booleannet2bnet(rules)
primes = bnet_text2primes(rules)
pyboolnet.prime_implicants.percolation(primes,True)
print("\nReduced prime implicant rules:")
sm.format.pretty_print_prime_rules(primes)
