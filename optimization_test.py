import PyBoolNet
import StableMotifs as sm
import networkx as nx
from timeit import default_timer

N=1000
K=2
p=sm.RandomBooleanNetworks.get_criticality_p_Kauffman(K)[0]
N_ensemble=10
seed=1000
rbn_ensemble=sm.RandomBooleanNetworks.Random_Boolean_Network_Ensemble_Kauffman(N,K,p,N_ensemble,seed=seed,write_Boolean_network=True)

start=default_timer()
rules=rbn_ensemble[0]

rules = sm.Format.booleannet2bnet(rules)
primes = sm.Format.longbnet2primes(rules)
PyBoolNet.PrimeImplicants._percolation(primes,True)
end=default_timer()
print("Time (s) creating reduced networks:",end-start)
print("Reduced network size: ",str(len(primes)))
sm.Format.pretty_print_prime_rules(primes)


start=default_timer()
diag = sm.Succession.build_succession_diagram(primes,max_simulate_size=20)
end=default_timer()
print("Time (s) finding atttactors:",end-start)
print("Number of attractors: ",str(len(diag.attractor_fixed_nodes_list)))
diag.attractor_candidate_summary(show_reduced_rules=False)
