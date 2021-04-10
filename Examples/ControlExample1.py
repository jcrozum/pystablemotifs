import PyBoolNet
import StableMotifs as sm
import networkx as nx
from timeit import default_timer
print("Loading network . . .")
primes = sm.Format.import_primes('./models/controlExample1.txt')
print("Network loaded.")
print()
print("RULES")
sm.Format.pretty_print_prime_rules({k:primes[k] for k in sorted(primes)})

print("Producing and exporting expanded network . . .")
EN = sm.Export.expanded_network(primes)
nx.write_graphml(EN,"EN_ControlExample1.graphml")
print("done.")

print()
print("Analyzing network . . .")
ar = sm.AttractorRepertoire.from_primes(primes)

print("Analysis complete.")
print()
ar.summary()


print()
print("Exporting Succession Diagram")
SD=sm.Export.networkx_succession_diagram_reduced_network_based(ar,include_attractors_in_diagram=True)
sm.Export.save_to_graphml(SD,"SD_ControlExample1")
print()

target = {'A': 1, 'B': 1, 'C': 1, 'D': 1, 'E': 1, 'F': 0, 'G': 0}

print("Brute-force search for knockout/knockins that achieve",target,". . .")
start=default_timer()
koki = sm.DomainOfInfluence.knock_to_partial_state(target,primes,max_drivers=2)
end=default_timer()
print("Time running brute-force search method:",end-start)
print("Sets found:")
for x in koki: print(x)

print()
print("GRASP search for knockout/knockins that achieve",target,". . .")
start=default_timer()
sols = sm.DomainOfInfluence.GRASP(target,primes,2000)
end=default_timer()
print("Time running GRASP search method:",end-start)
print("Control sets that fix",target)
for x in sorted(sols,key=lambda x: len(x)):
    print(x)

print()
print("Building succession diagram . . .")
diag = ar.succession_diagram

print("Computing driver sets (in multiple ways) that reprogram to an attractor with ",target,". . .")
start=default_timer()
reprogram_sets_minimal = diag.reprogram_to_trap_spaces(target,target_method='merge',driver_method='minimal',max_drivers=None)
end=default_timer()
print()
print("Time running minimal merge method:",end-start)
print("Sets found:")
for x in reprogram_sets_minimal: print(x)

start=default_timer()
reprogram_sets_merge = diag.reprogram_to_trap_spaces(target,target_method='merge',driver_method='internal',max_drivers=None)
end=default_timer()
print()
print("Time running internal merge method:",end-start)
print("Sets found:")
for x in reprogram_sets_merge: print(x)

start=default_timer()
reprogram_sets_GRASP_merge = diag.reprogram_to_trap_spaces(target,target_method='merge',driver_method='GRASP',GRASP_iterations=None)
end=default_timer()
print()
print("Time running GRASP merge method:",end-start)
print("Sets found:")
for x in sorted(reprogram_sets_GRASP_merge,key=lambda x: len(x)): print(x)

start=default_timer()
reprogram_sets_history = diag.reprogram_to_trap_spaces(target,target_method='history',driver_method='internal',max_drivers=None)
end=default_timer()
print()
print("Time running internal history method:",end-start)
print("Sets found:")
for x in reprogram_sets_history: print(x)

start=default_timer()
reprogram_sets_minimal_history = diag.reprogram_to_trap_spaces(target,target_method='history',driver_method='minimal',max_drivers=None)
end=default_timer()
print()
print("Time running minimal history method:",end-start)
print("Sets found:")
for x in reprogram_sets_minimal_history:
    print("---")
    print("One temporary intervention from each list, in order.")
    for y in x: print(y,"\n")

start=default_timer()
reprogram_sets_GRASP_history = diag.reprogram_to_trap_spaces(target,target_method='history',driver_method='GRASP',GRASP_iterations=None)
end=default_timer()
print()
print("Time running GRASP history method:",end-start)
print("Sets found:")
for x in reprogram_sets_GRASP_history:
    print("---")
    print("One temporary intervention from each list, in order.")
    for y in x: print(y,"\n")
