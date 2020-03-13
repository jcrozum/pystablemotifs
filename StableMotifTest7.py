import PyBoolNet
import StableMotifs as sm
from timeit import default_timer
import pickle

print("Loading EMT network . . .")
primes = sm.Format.import_primes("models/EMT.txt",remove_constants=True)
reprogramming_target = {'EMT':0}


print("Brute-force search for knockout/knockins that achieve",reprogramming_target,". . .")
start=default_timer()
koki = sm.DomainOfInfluence.knock_to_partial_state(reprogramming_target,primes,max_drivers=2)
end=default_timer()
print("Time running brute-force search method:",end-start)
print("Sets found:")
for x in koki: print(x)


print("Building succession diagram . . .")
# We do not need the complex attractors for this example, and ruling them
# out in the EMT case is extremely slow
diag = sm.Succession.build_succession_diagram(primes)

print("Computing driver sets (in multiple ways) that reprogram to an attractor with ",reprogramming_target,". . .")
start=default_timer()
reprogram_sets_minimal = diag.reprogram_to_trap_spaces(reprogramming_target,method='minimal',max_drivers=4)
end=default_timer()
print("Time running minimal method:",end-start)
print("Sets found:")
for x in reprogram_sets_minimal: print(x)

start=default_timer()
reprogram_sets_merge = diag.reprogram_to_trap_spaces(reprogramming_target,method='merge',max_drivers=4)
end=default_timer()
print("Time running merge method:",end-start)
print("Sets found:")
for x in reprogram_sets_merge: print(x)

start=default_timer()
reprogram_sets_history = diag.reprogram_to_trap_spaces(reprogramming_target,method='history',max_drivers=4)
end=default_timer()
print("Time running history method:",end-start)
print("Sets found:")
for x in reprogram_sets_history: print(x)
#
start=default_timer()
reprogram_sets_minimal_history = diag.reprogram_to_trap_spaces(reprogramming_target,method='minimal_history',max_drivers=2)
end=default_timer()
print("Time running minimal_history method:",end-start)
print("Sets found:")
for x in reprogram_sets_minimal_history: print(x)
