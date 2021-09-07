import pyboolnet
import pystablemotifs as sm
import networkx as nx

print("Loading network . . .")
model = 'test5'
primes = sm.format.import_primes('./models/'+model+'.txt')
print("Network loaded.")
print()
print("RULES")
sm.format.pretty_print_prime_rules({k:primes[k] for k in sorted(primes)})

stable_motifs = pyboolnet.trap_spaces.trap_spaces(primes, "min",MaxOutput=10000)
print(stable_motifs)

import time
start_time = time.time()

print()
print("Analyzing network . . .")
ar = sm.AttractorRepertoire.from_primes(primes,max_simulate_size=20)
print("Analysis complete.")
print()
ar.summary()

for i in range(0,len(ar.attractors)):
    if ar.attractors[i].stg is not None:
        print("stg format information")
        print("1. attractor")
        print(ar.attractors[i].attractor_dict)
        print("2. type of the stg object")
        print(type(ar.attractors[i].stg))
        print("3. attribute's name, type, and value.")
        for attribute in ar.attractors[i].stg.__dict__.keys():
            print(str(attribute)+" : "+str(type(ar.attractors[i].stg.__dict__[attribute])))
            print(ar.attractors[i].stg.__dict__[attribute])

        nx.write_graphml(ar.attractors[i].stg, model+"_"+str(i)+"_GAU.graphml")

print("--- %s seconds for non-MBPN ---" % (time.time() - start_time))
start_time = time.time()

print()
print("Analyzing network (MPBN) . . .")
ax = sm.AttractorRepertoire.from_primes(primes,max_simulate_size=20,MPBN_update=True)
print("Analysis complete.")
print()
ax.summary()

for i in range(0,len(ax.attractors)):
    if ax.attractors[i].stg is not None:
        print("stg format information")
        print("1. attractor")
        print(ax.attractors[i].attractor_dict)
        print("2. type of the stg object")
        print(type(ax.attractors[i].stg))
        print("3. attribute's name, type, and value.")
        for attribute in ax.attractors[i].stg.__dict__.keys():
            print(str(attribute)+" : "+str(type(ax.attractors[i].stg.__dict__[attribute])))
            print(ax.attractors[i].stg.__dict__[attribute])

        nx.write_graphml(ax.attractors[i].stg, model+"_"+str(i)+"_MPBN.graphml")

print("--- %s seconds for MPBN ---" % (time.time() - start_time))

#
# print("Time reversal:")
# tr = sm.time_reversal.time_reverse_primes(primes)
#
# sm.format.pretty_print_prime_rules(tr)
#
# print("\nTR Stable Notifs:")
#
# for x in ar.succession_diagram.motif_reduction_dict[0].stable_motifs:
#     print(x)
#
# print("Producing and exporting STG . . .")
# G=pyboolnet.state_transition_graphs.primes2stg(primes,'asynchronous')
# H=nx.DiGraph()
# for u,v in G.edges():
#     H.add_edge(u,v)
# for n in H.nodes():
#     H.nodes[n]["label"]=n
# nx.write_graphml(H,"Fig3STG.graphml")
# print("done.")
#
# print("Producing and exporting expanded network . . .")
# K = sm.export.expanded_network(primes)
# KTR = sm.export.expanded_network(sm.time_reversal.time_reverse_primes(primes))
# nx.write_graphml(K,"EN_PhaseSwitch.graphml")
# nx.write_graphml(KTR,"EN_PhaseSwitch_TR.graphml")
# print("done.")
