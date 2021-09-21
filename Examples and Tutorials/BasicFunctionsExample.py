import pyboolnet.state_transition_graphs
import pystablemotifs as sm
import networkx as nx
print("Loading network . . .")
primes = sm.format.import_primes('./models/test2.txt')
print("Network loaded.")
print()
print("RULES")
sm.format.pretty_print_prime_rules({k:primes[k] for k in sorted(primes)})
print()
print("Analyzing network . . .")
ar = sm.AttractorRepertoire.from_primes(primes,max_simulate_size=7)

print("Analysis complete.")
print()
ar.summary()


print("Time reversal:")
tr = sm.time_reversal.time_reverse_primes(primes)

sm.format.pretty_print_prime_rules(tr)

print("\nTR Stable Notifs:")

for x in ar.succession_diagram.motif_reduction_dict[0].stable_motifs:
    print(x)

print("Producing and exporting STG . . .")
G=pyboolnet.state_transition_graphs.primes2stg(primes,'asynchronous')
H=nx.DiGraph()
for u,v in G.edges():
    H.add_edge(u,v)
for n in H.nodes():
    H.nodes[n]["label"]=n
nx.write_graphml(H,"Fig3STG.graphml")
print("done.")
