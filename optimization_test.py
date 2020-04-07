import PyBoolNet
import StableMotifs as sm
import networkx as nx
from timeit import default_timer


# primes = sm.Format.import_primes("./models/test9.txt")
#
# #new_primes = sm.Reduction.delete_node(primes,'Z2')
#
# red_primes = sm.Reduction.deletion_reduction(primes)
# sm.Format.pretty_print_prime_rules(primes)
# print()
# #sm.Format.pretty_print_prime_rules(new_primes)
# print()
# print('reduction order',sorted(primes))
# sm.Format.pretty_print_prime_rules(red_primes)
#
# print()
# print()
# print()
# print()
# print()
# primes = sm.Format.import_primes("./models/test10.txt")
#
# #new_primes = sm.Reduction.delete_node(primes,'Z2')
#
# red_primes = sm.Reduction.mediator_reduction(primes)
# sm.Format.pretty_print_prime_rules(primes)
# print()
# #sm.Format.pretty_print_prime_rules(new_primes)
# print()
# sm.Format.pretty_print_prime_rules(red_primes)
# rules='''
# n100* = n68
# n111* = !n69
# n112* = !n50
# n121* = !n146 & n305
# n127* = !n417
# n128* = !n22
# n13* = n457 | !n338
# n137* = !n386 | !n217
# n141* = !n298
# n146* = !n477
# n163* = n316
# n170* = n127
# n171* = n13
# n19* = !n284
# n191* = !n432
# n201* = n308
# n203* = n296
# n208* = !n32
# n21* = !n260
# n211* = !n484
# n217* = !n146
# n22* = !n203
# n228* = !n295
# n230* = n335
# n238* = !n248
# n248* = !n361
# n254* = !n323
# n258* = n350
# n26* = !n13
# n260* = n476
# n267* = !n404
# n269* = !n409
# n274* = !n35
# n284* = n260
# n295* = n141
# n296* = n208
# n298* = n350
# n304* = !n376
# n305* = n432
# n308* = n295 & n416
# n314* = !n392
# n316* = n111
# n32* = !n324
# n322* = !n476
# n323* = n208
# n324* = n485 & n98
# n327* = n141
# n334* = n304
# n335* = !n228 | !n171
# n338* = n83
# n346* = !n203 & n485
# n35* = n68
# n350* = !n170
# n355* = !n203 & !n83
# n358* = !n457
# n361* = !n141
# n364* = !n22 & !n26 | n22 & n26
# n376* = !n324
# n378* = n346
# n386* = !n305 | n267
# n392* = n296
# n394* = !n203
# n404* = !n211
# n406* = !n170
# n409* = n248 & !n443
# n416* = n493
# n417* = n304
# n423* = n83
# n432* = n98 | !n112
# n437* = n81
# n44* = n35 & n474
# n442* = !n203
# n443* = !n69
# n454* = !n121
# n457* = !n21
# n458* = !n146
# n474* = !n324
# n475* = n493
# n476* = !n308
# n477* = !n203
# n478* = n394
# n484* = n35
# n485* = !n322 & n35
# n487* = !n361
# n492* = !n127
# n493* = !n137
# n50* = n274
# n52* = n35
# n55* = !n146
# n68* = !n128
# n69* = !n474
# n73* = !n392
# n78* = !n55
# n8* = n485
# n81* = !n386
# n83* = n98
# n98* = !n121
# '''
#
#
#
# rules_pbn = sm.Format.booleannet2bnet(rules)
# primes = PyBoolNet.FileExchange.bnet2primes(rules_pbn)
# print("RBN Example:")
# sm.Format.pretty_print_prime_rules(primes)
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

N=500
K=2
p=sm.randomBooleanNetworks.get_criticality_p_Kauffman(K)[0]
N_ensemble=1
seed=1000
rbn_ensemble=sm.randomBooleanNetworks.Random_Boolean_Network_Ensemble_Kauffman(N,K,p,N_ensemble,seed=seed,write_Boolean_network=True)

start=default_timer()
rules=rbn_ensemble[0]

rules = sm.Format.booleannet2bnet(rules)
# print(rules)
primes = PyBoolNet.FileExchange.bnet2primes(rules)
PyBoolNet.PrimeImplicants._percolation(primes,True)
primes = sm.Reduction.mediator_reduction(primes)
primes = sm.Reduction.deletion_reduction(primes)
end=default_timer()
print("Time (s) creating reduced networks:",end-start)
print("Reduced network size: ",str(len(primes)))
sm.Format.pretty_print_prime_rules(primes)

G = PyBoolNet.InteractionGraphs.primes2igraph(primes)
C = nx.condensation(G)
for x in nx.topological_sort(C):
    print("===========================")
    print(len(C.nodes[x]['members'])," nodes in this SCC")
    sm.Format.pretty_print_prime_rules({k:v for k,v in primes.items() if k in C.nodes[x]['members']})
    print("===========================")
print(nx.number_weakly_connected_components(G))
print(nx.number_strongly_connected_components(G))
print(len(list(C.edges)))
od = PyBoolNet.InteractionGraphs.find_outdag(G)
print(len(od))
print(len(primes))

# quit()

start=default_timer()
diag = sm.Succession.build_succession_diagram(primes,search_partial_STGs=True)
end=default_timer()
print("Time (s) finding atttactors:",end-start)
print("Number of attractors: ",str(len(diag.attractor_fixed_nodes_list)))
diag.attractor_candidate_summary()
