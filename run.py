import PyBoolNet
import StableMotifs as sm
import sys
from timeit import default_timer

# import and format rules
primes = sm.Format.import_primes(sys.argv[1],remove_constants=True)

# print rules after reduction
sm.Format.pretty_print_prime_rules(primes)

diag = sm.Succession.build_succession_diagram(primes)
diag.summary()

diag.attractor_candidate_summary()

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
