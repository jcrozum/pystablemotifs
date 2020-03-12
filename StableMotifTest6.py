import PyBoolNet
import StableMotifs as sm

def parse(stream):
    "Ingores commented out lines"
    lines = list(stream)
    lines = filter(lambda x: not x.startswith("#"), lines)
    rules = "\n".join(lines)
    return rules


rules='''
xA*= not xA and not xB or xC
xB*= not xA and not xB or xC
xC*= xA and xB
'''

rules_pbn = sm.Format.booleannet2bnet(rules)
primes = PyBoolNet.FileExchange.bnet2primes(rules_pbn)
print("\"Pathological\" Example:")
sm.Format.pretty_print_prime_rules(primes)
print()
motif = {'xA':1,'xB':1,'xC':1}
drivers = sm.DomainOfInfluence.find_internal_motif_drivers(motif,primes)
print("The following sets are internal driver sets of the stable motif", motif)
print(drivers)
diag = sm.Succession.build_succession_diagram(primes)

print("Driver sets that reprogram to an attractor with xA=1:")
print(diag.reprogram_to_trap_spaces({'xA':1}))

# now load in TLGL network
print("Loading TLGL Network")
rules = parse(open("models/test3.txt"))

# Reformat rules and get rid of constants
rules_pbn = sm.Format.booleannet2bnet(rules)
primes = PyBoolNet.FileExchange.bnet2primes(rules_pbn)
PyBoolNet.PrimeImplicants._percolation(primes,True)

diag = sm.Succession.build_succession_diagram(primes)
print("Driver sets that reprogram to an attractor with Apoptosis=0:")
reprogram_sets = diag.reprogram_to_trap_spaces({'Apoptosis':0})
print(reprogram_sets)
for set_size in range(len(reprogram_sets[-1])+1):
    n = len([x for x in reprogram_sets if len(x)==set_size])
    print("There are",n,"driver sets of size",set_size)

print("1-node driver sets for Apoptosis=0:")
print([x for x in reprogram_sets if len(x)==1])

print("2-node driver sets for Apoptosis=0:")
print([x for x in reprogram_sets if len(x)==2])
