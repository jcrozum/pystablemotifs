import PyBoolNet
import StableMotifs as sm
from timeit import default_timer

print("Loading TLGL network and creating succssion diagram . . .")
primes = sm.Format.import_primes("models/test3.txt",remove_constants=True)

diag = sm.Succession.build_succession_diagram(primes)

print("Computing driver sets (in multiple ways) that reprogram to an attractor with Apoptosis=0 . . .")
start=default_timer()
reprogram_sets_minimal = diag.reprogram_to_trap_spaces({'Apoptosis':0},method='minimal')
end=default_timer()
print("Time running minimal method:",end-start)

start=default_timer()
reprogram_sets_merge = diag.reprogram_to_trap_spaces({'Apoptosis':0},method='merge')
end=default_timer()
print("Time running merge method:",end-start)

start=default_timer()
reprogram_sets_history = diag.reprogram_to_trap_spaces({'Apoptosis':0},method='history')
end=default_timer()
print("Time running history method:",end-start)

start=default_timer()
reprogram_sets_minimal_history = diag.reprogram_to_trap_spaces({'Apoptosis':0},method='minimal_history')
end=default_timer()
print("Time running minimal_history method:",end-start)


if len(reprogram_sets_merge) == len(reprogram_sets_history) and all([x in reprogram_sets_merge for x in reprogram_sets_history]):
    reprogram_sets = reprogram_sets_merge
else:
    print("Something went horribly wrong!!")
    print("MERGE",reprogram_sets_merge)
    print("HISTORY",reprogram_sets_history)
    print("MINIMAL",reprogram_sets_minimal)
    quit()
print()
print("Internal drivers found:")
print(reprogram_sets)
for f in reprogram_sets:
    imp,con = sm.DomainOfInfluence.logical_domain_of_influence(f,primes)
    set_is_good = False
    if 'Apoptosis' in imp:
        set_is_good = imp['Apoptosis'] == 0

    if set_is_good:
        print("Set", f, "verified as driver set of Apoptosis=0.")
    else:
        print("Set", f, "failed to drive Apoptosis=0 as permanent intervention (cannot rule out success though).")
print()
print("Minimal drivers found:")
print(reprogram_sets_minimal)
for f in reprogram_sets_minimal:
    imp,con = sm.DomainOfInfluence.logical_domain_of_influence(f,primes)
    set_is_good = False
    if 'Apoptosis' in imp:
        set_is_good = imp['Apoptosis'] == 0

    if set_is_good:
        print("Set", f, "verified as driver set of Apoptosis=0.")
    else:
        print("Set", f, "failed to drive Apoptosis=0 as permanent intervention (cannot rule out success though).")

print()
print("Minimal hitory drivers found:")
for x in reprogram_sets_minimal_history: print(x)
