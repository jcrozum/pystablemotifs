import PyBoolNet
import StableMotifs as sm
from timeit import default_timer

print("Loading TLGL network and creating succssion diagram . . .")
primes = sm.Format.import_primes("models/test3.txt",remove_constants=True)
reprogramming_target = {'Apoptosis':1}
diag = sm.Succession.build_succession_diagram(primes)

print("Computing driver sets (in multiple ways) that reprogram to an attractor with ",reprogramming_target,". . .")
start=default_timer()
reprogram_sets_minimal = diag.reprogram_to_trap_spaces(reprogramming_target,target_method='history',driver_method='internal',max_drivers=None)
end=default_timer()
print()
print("Time running minimal method:",end-start)
print("Sets found:")
for x in reprogram_sets_minimal: print(x)

start=default_timer()
reprogram_sets_merge = diag.reprogram_to_trap_spaces(reprogramming_target,target_method='merge',driver_method='internal',max_drivers=None)
end=default_timer()
print()
print("Time running merge method:",end-start)
print("Sets found:")
for x in reprogram_sets_merge: print(x)

start=default_timer()
reprogram_sets_history = diag.reprogram_to_trap_spaces(reprogramming_target,target_method='history',driver_method='internal',max_drivers=None)
end=default_timer()
print()
print("Time running history method:",end-start)
print("Sets found:")
for x in reprogram_sets_history: print(x)
#
start=default_timer()
reprogram_sets_minimal_history = diag.reprogram_to_trap_spaces(reprogramming_target,target_method='history',driver_method='minimal',max_drivers=None)
end=default_timer()
print()
print("Time running minimal_history method:",end-start)
print("Sets found:")
for x in reprogram_sets_minimal_history:
    print("---")
    print("One temporary intervention from each list, in order.")
    for y in x: print(y,"\n")
