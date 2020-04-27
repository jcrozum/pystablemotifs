import PyBoolNet
import StableMotifs as sm
import sys
from timeit import default_timer

# import and format rules
primes = sm.Format.import_primes(sys.argv[1],remove_constants=True)

# print rules after reduction
sm.Format.pretty_print_prime_rules({k:primes[k] for k in sorted(primes)})

diag = sm.Succession.build_succession_diagram(primes)
diag.summary()

diag.attractor_candidate_summary()

# Uncomment to write attractors and succession diagrams to file. Requires networkx and pandas libraries
print("Writing succession diagram network and attractors to file . . .")
import networkx as nx
diag.process_networkx_succession_diagram(include_attractors_in_diagram=True)
nx.write_graphml(diag.G_reduced_network_based_labeled, "SuccessionDiagram_reduced_network_based.graphml")
nx.write_gml(diag.G_reduced_network_based_labeled, "SuccessionDiagram_reduced_network_based.gml")
nx.write_graphml(diag.G_motif_based_labeled, "SuccessionDiagram_motif_based.graphml")
nx.write_gml(diag.G_motif_based_labeled, "SuccessionDiagram_motif_based.gml")

import pandas as pd
df_attractors=pd.DataFrame(columns=[node for node in diag.unreduced_primes])
df_attractors=df_attractors.merge(pd.DataFrame(diag.attractor_fixed_nodes_list),how="outer")
df_attractors.index=["Attractor_"+str(i) for i,value in enumerate(diag.attractor_fixed_nodes_list)]
df_attractors=df_attractors.fillna("X").astype(str).replace({"0.0": "0", "1.0": "1"}).sort_index(axis=1)
df_attractors.to_csv("Attractors.csv")

reprogramming_target = None
for reduction in diag.motif_reduction_list:
    if reduction.terminal == "yes" and len(reduction.logically_fixed_nodes) > 0:
        reprogramming_target = reduction.logically_fixed_nodes
        break
print()
print()
print()
if reprogramming_target is None:
    print("Could not find a terminal network reduction.")
else:
    print("We will control the system towards a state in which the following node states are fixed:")
    print(reprogramming_target)

print("Brute-force search for knockout/knockins that achieve",reprogramming_target,". . .")
start=default_timer()
koki = sm.DomainOfInfluence.knock_to_partial_state(reprogramming_target,primes,max_drivers=2)
end=default_timer()
print("Time running brute-force search method:",end-start)
print("Sets found:")
for x in koki: print(x)
print()
print("Computing driver sets (in multiple ways) . . .")
start=default_timer()
reprogram_sets_minimal = diag.reprogram_to_trap_spaces(reprogramming_target,target_method='merge',driver_method='minimal',max_drivers=4)
end=default_timer()
print("Time running minimal method:",end-start)
print("Sets found:")
for x in reprogram_sets_minimal: print(x)
print()
start=default_timer()
reprogram_sets_merge = diag.reprogram_to_trap_spaces(reprogramming_target,target_method='merge',driver_method='internal',max_drivers=4)
end=default_timer()
print("Time running merge method:",end-start)
print("Sets found:")
for x in reprogram_sets_merge: print(x)
print()
start=default_timer()
reprogram_sets_history = diag.reprogram_to_trap_spaces(reprogramming_target,target_method='history',driver_method='internal',max_drivers=4)
end=default_timer()
print("Time running history method:",end-start)
print("Sets found:")
for x in reprogram_sets_history: print(x)
#
start=default_timer()
reprogram_sets_minimal_history = diag.reprogram_to_trap_spaces(reprogramming_target,target_method='history',driver_method='minimal',max_drivers=None)
end=default_timer()
print()
print("Time running minimal_history method:",end-start)
print("-----")
print("Sets found:")
for x in reprogram_sets_minimal_history:
    print("One temporary intervention from each list, in order.")
    for y in x: print(y,"\n")
    print("---")
