import PyBoolNet
import StableMotifs as sm
import networkx as nx
from timeit import default_timer
#
# rules = """xA* = xB
# xB* = xC
# xC* = xA & !xB
# """
# rules_pbn = sm.Format.booleannet2bnet(rules)
# primes = PyBoolNet.FileExchange.bnet2primes(rules_pbn)
# sm.Format.pretty_print_prime_rules(primes)
#
# print("========================")
# print("Rules:")
# rp,rc = sm.Reduction.deletion_reduction(primes)
# sm.Format.pretty_print_prime_rules(rp)
# print("Constants:")
# print(rc)
#
# G = PyBoolNet.InteractionGraphs.primes2igraph(primes)
# GC = nx.condensation(G)
# print(nx.number_weakly_connected_components(G))
# print(nx.number_strongly_connected_components(G))
# print(list(GC.nodes))
# od = PyBoolNet.InteractionGraphs.find_outdag(G)
# print(len(od))
# print(len(primes))
# print()
#print("Building diagram . . .")
#diag = sm.Succession.build_succession_diagram(primes)
# print("ATTRACTOR SUMMARY")
# print()
# diag.attractor_candidate_summary()

N=1500
K=2
p=sm.randomBooleanNetworks.get_criticality_p_Kauffman(K)[0]
N_ensemble=1
seed=1000
rbn_ensemble=sm.randomBooleanNetworks.Random_Boolean_Network_Ensemble_Kauffman(N,K,p,N_ensemble,seed=seed,write_Boolean_network=True)

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
