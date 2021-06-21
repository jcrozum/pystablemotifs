import PyBoolNet
import PyStableMotifs as sm
import networkx as nx
print("Loading network . . .")
primes = sm.Format.import_primes('./models/TLGL_Small.txt',remove_constants=True)
print("Network loaded.")
print()
print("RULES")
sm.Format.pretty_print_prime_rules({k:primes[k] for k in sorted(primes)})
print()
print("Analyzing network . . .")
ar = sm.AttractorRepertoire.from_primes(primes,max_simulate_size=7)
ar.summary()
root = ar.succession_diagram.motif_reduction_dict[0]

print("Producing and exporting STG . . .")
G=PyBoolNet.StateTransitionGraphs.primes2stg(primes,'asynchronous')
H=nx.DiGraph()
for u,v in G.edges():
    H.add_edge(u,v)
for n in H.nodes():
    H.nodes[n]["label"]=n
nx.write_graphml(H,"STG_TLGL_Small.graphml")
print("done.")

print("Producing and exporting expanded network . . .")
K = sm.Export.expanded_network(primes)
KTR = sm.Export.expanded_network(sm.TimeReversal.time_reverse_primes(primes))
nx.write_graphml(K,"EN_TLGL_Small.graphml")
nx.write_graphml(KTR,"EN_TR_TLGL_Small.graphml")
print("done.")

print("\nTime-reversal analysis")
print("\nGarden Spaces:")
print(root.time_reverse_stable_motifs)
# outputs: [{'cer': 1, 's1p': 1}, {'apo': 0}]

print("Thus, the state-space is partitioned into four regions:")
print("1) cer & s1p = 1, apo = 0")
print("2) cer & s1p = 1, apo = 1")
print("3) cer & s1p = 0, apo = 0")
print("4) cer & s1p = 0, apo = 1")
print("and no attractor lies in more than one region.")
print("Note that these are NOT the Garden spaces, which are instead regions 1, 1+2, and 1+3.")

print("\nConsider regions 1 and 2 and the LDOI of {(cer,1),(s1p,1)}:")
imp,con = sm.DomainOfInfluence.logical_domain_of_influence({'cer': 1, 's1p': 1},primes)
print("The contradiction boundary is non-empty:",con)
print("Therefore, neither region 1 nor 2 contains any attractor.")

print("\nNext, consider region 4 and the LDOI of {(apo,1)}:")
impss1,con = sm.DomainOfInfluence.logical_domain_of_influence({'apo': 1},primes)
print("It contains",impss1)
print("This fixes all values, so region 4 contains this fixed point and no other attractors.")

print("\nAll that remains is to consider region 3, in which apo = 0.")
print("The condition cer & s1p = 0 can be satisfied for s1p = 0 or for cer = 0 (or both).")
imp,con = sm.DomainOfInfluence.logical_domain_of_influence({'s1p': 0,'apo':0},primes)
print("However, {(s1p,0),(apo,0)} has non-empty contradiction boundary:",con)
impss2,con = sm.DomainOfInfluence.logical_domain_of_influence({'cer': 0,'apo':0},primes)
print("Meanwhile, {(cer,0),(apo,0)} implies a steady state:",impss2)

print("\nThus, the only two attractors are in regions 3 and 4 and are as follows:")
print({k:impss1[k] for k in sorted(impss1)})
print({k:impss2[k] for k in sorted(impss2)})
