import PyBoolNet
import PyStableMotifs as sm
import networkx as nx
print("Loading network . . .")
primes = sm.Format.import_primes('./models/test2.txt')
print("Network loaded.")
print()
print("RULES")
sm.Format.pretty_print_prime_rules({k:primes[k] for k in sorted(primes)})
print()
print("Analyzing network . . .")
ar = sm.AttractorRepertoire.from_primes(primes,max_simulate_size=7)

print("Analysis complete.")
print()
ar.summary()


print("Time reversal:")
tr = sm.TimeReversal.time_reverse_primes(primes)

sm.Format.pretty_print_prime_rules(tr)

print("\nTR Stable Notifs:")

for x in ar.succession_diagram.motif_reduction_dict[0].stable_motifs:
    print(x)

print("Producing and exporting STG . . .")
G=PyBoolNet.StateTransitionGraphs.primes2stg(primes,'asynchronous')
H=nx.DiGraph()
for u,v in G.edges():
    H.add_edge(u,v)
for n in H.nodes():
    H.nodes[n]["label"]=n
nx.write_graphml(H,"Fig3STG.graphml")
print("done.")

# print("Producing and exporting expanded network . . .")
# K = sm.Export.expanded_network(primes)
# KTR = sm.Export.expanded_network(sm.TimeReversal.time_reverse_primes(primes))
# nx.write_graphml(K,"EN_PhaseSwitch.graphml")
# nx.write_graphml(KTR,"EN_PhaseSwitch_TR.graphml")
# print("done.")
