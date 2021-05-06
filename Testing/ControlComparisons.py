# This tests 6 attractor control methods and 2 target control methods.
# GRASP methods do not give determinstic results.

import PyBoolNet
import StableMotifs as sm
from timeit import default_timer


def compare_control_methods(ar,target):
    diag = ar.succession_diagram

    print("Brute-force search for knockout/knockins that achieve",target,". . .")
    start=default_timer()
    koki = sm.DomainOfInfluence.knock_to_partial_state(target,ar.primes,max_drivers=2)
    end=default_timer()
    print("Time running brute-force search method:",end-start)
    print("Sets found:")
    for x in koki: print(x)

    print()
    print("GRASP search for knockout/knockins that achieve",target,". . .")
    start=default_timer()
    sols = sm.DomainOfInfluence.GRASP(target,ar.primes,2000)
    end=default_timer()
    print("Time running GRASP search method:",end-start)
    print("Control sets that fix",target)
    for x in sorted(sols,key=lambda x: len(x)):
        print(x)
    print()
    print("=-"*3)
    print("\nComputing driver sets (in multiple ways) that reprogram to an attractor with ",target,". . .")
    start=default_timer()
    reprogram_sets_minimal = diag.reprogram_to_trap_spaces(target,
        target_method='merge',driver_method='minimal',max_drivers=4)
    end=default_timer()
    print()
    print("Time running minimal merge method:",end-start)
    print("Sets found:")
    for x in reprogram_sets_minimal: print(x)

    start=default_timer()
    reprogram_sets_merge = diag.reprogram_to_trap_spaces(target,
        target_method='merge',driver_method='internal',max_drivers=4)
    end=default_timer()
    print()
    print("Time running internal merge method:",end-start)
    print("Sets found:")
    for x in reprogram_sets_merge: print(x)

    start=default_timer()
    reprogram_sets_GRASP_merge = diag.reprogram_to_trap_spaces(target,
        target_method='merge',driver_method='GRASP',GRASP_iterations=2000)
    end=default_timer()
    print()
    print("Time running GRASP merge method:",end-start)
    print("Sets found:")
    for x in sorted(reprogram_sets_GRASP_merge,key=lambda x: len(x)): print(x)

    start=default_timer()
    reprogram_sets_history = diag.reprogram_to_trap_spaces(target,
        target_method='history',driver_method='internal',max_drivers=4)
    end=default_timer()
    print()
    print("Time running internal history method:",end-start)
    print("Sets found:")
    for x in reprogram_sets_history: print(x)

    start=default_timer()
    reprogram_sets_minimal_history = diag.reprogram_to_trap_spaces(target,
        target_method='history',driver_method='minimal',max_drivers=4)
    end=default_timer()
    print()
    print("Time running minimal history method:",end-start)
    print("Sets found:")
    for x in reprogram_sets_minimal_history:
        print("---")
        print("One temporary intervention from each list, in order.")
        print("("+str(len(x))+"interventions in total)")
        for y in x: print(y,"\n")

    start=default_timer()
    reprogram_sets_GRASP_history = diag.reprogram_to_trap_spaces(target,
        target_method='history',driver_method='GRASP',GRASP_iterations=500)
    end=default_timer()
    print()
    print("Time running GRASP history method:",end-start)
    print("Sets found:")
    for x in reprogram_sets_GRASP_history:
        print("---")
        print("One temporary intervention from each list, in order.")
        print("("+str(len(x))+"interventions in total)")
        for y in x: print(y,"\n")


print("Control comparisions: driving the T-LGL network to apoptosis.")
print('. . .')
primes = sm.Format.import_primes("../models/Control Benchmarks/TLGL_Large_Fixed_Inputs.txt",remove_constants=True)
ar = sm.AttractorRepertoire.from_primes(primes)
print("Attractor repertoire constructed.")
print('-'*10)
target = {'Apoptosis':1}
compare_control_methods(ar,target)
print('='*20)


print("Control comparisions: driving the EMT network to the epithelial state.")
print('. . .')
primes = sm.Format.import_primes("../models/EMT.txt",remove_constants=True)
ar = sm.AttractorRepertoire.from_primes(primes)
print("Attractor repertoire constructed.")
print('-'*10)
target = {'EMT':0}
compare_control_methods(ar,target)
print('='*20)
