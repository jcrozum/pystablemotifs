import PyBoolNet
import StableMotifs as sm
from timeit import default_timer

print("Loading EMT network . . .")
primes = sm.Format.import_primes("models/EMT.txt",remove_constants=True)
target = {'EMT':0}


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
sols = sm.DomainOfInfluence.GRASP(target,primes,5000)
end=default_timer()
print("Time running GRASP search method:",end-start)
print("GRASP control sets found:")
for x in sorted(sols,key=lambda x: len(x)):
    print(x)
print()
print("Building succession diagram . . .")
diag = sm.Succession.build_succession_diagram(primes)

print("Computing driver sets (in multiple ways) that reprogram to an attractor with ",target,". . .")


start=default_timer()
reprogram_sets_GRASP_merge = diag.reprogram_to_trap_spaces(target,target_method='merge',driver_method='GRASP',max_iterations=20000)
end=default_timer()
print()
print("Time running GRASP merge method:",end-start)
print("Sets found:")
for x in sorted(reprogram_sets_GRASP_merge,key=lambda x: len(x)): print(x)
print()
start=default_timer()
reprogram_sets_GRASP_history = diag.reprogram_to_trap_spaces(target,target_method='history',driver_method='GRASP',max_iterations=20000)
end=default_timer()
print()
print("Time running GRASP history method:",end-start)
print("Sets found:")
for x in reprogram_sets_GRASP_history:
    print("---")
    print("One temporary intervention from each list, in order.")
    for y in x: print(y,"\n")
