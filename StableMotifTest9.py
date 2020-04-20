import PyBoolNet
import StableMotifs as sm
from timeit import default_timer
import pickle

print("Loading EMT network . . .")
primes = sm.Format.import_primes("models/EMT.txt",remove_constants=True)
reprogramming_target = {'EMT':0}

print("Building succession diagram . . .")
diag = sm.Succession.build_succession_diagram(primes)

print("Computing attractor control sets using internal stable motif drivers.")
start=default_timer()
reprogram_sets_history = diag.reprogram_to_trap_spaces(reprogramming_target,target_method='history',driver_method='internal',max_drivers=4)
end=default_timer()
print("Time running history method:",end-start)
print("Sets found:")
for x in reprogram_sets_history: print(x)


print("Looking for knockout/knockins that achieve the controls . . .")
allowed_permanent = []
for x in reprogram_sets_history:
    print("Searching for smaller driver set for",x)
    koki = sm.DomainOfInfluence.knock_to_partial_state(x,primes,max_drivers=len(x)-1)


    print(len(koki),"smaller driver sets found.")
    for y in koki:
        imp,con = sm.DomainOfInfluence.logical_domain_of_influence(y,primes)
        if imp.items() >= reprogramming_target.items():
            print("Permanent intervention:",y)
        else:
            print("Temporary intervention:",y)
