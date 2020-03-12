import PyBoolNet
import StableMotifs as sm
# import networkx as nx

# rules='''
# xA*= not xA and not xB or xC
# xB*= not xA and not xB or xC
# xC*= xA and xB
# '''
#
# rules_pbn = sm.Format.booleannet2bnet(rules)
# primes = PyBoolNet.FileExchange.bnet2primes(rules_pbn)
# print("\"Pathological\" Example:")
# sm.Format.pretty_print_prime_rules(primes)
# print()
# motif = {'xA':1,'xB':1,'xC':1}
# drivers = sm.DomainOfInfluence.find_internal_motif_drivers(motif,primes)
# print("The following sets are internal driver sets of the stable motif", motif)
# print(drivers)
# diag = sm.Succession.build_succession_diagram(primes)
#
# print("Driver sets that reprogram to an attractor with xA=1:")
# print(diag.reprogram_to_trap_spaces({'xA':1}))

# now load in TLGL network
# print("Loading TLGL network and creating succssion diagram . . .")
# primes = sm.Format.import_primes("models/test3.txt",remove_constants=True)
#
# diag = sm.Succession.build_succession_diagram(primes)
# print(diag.edge_list)

# print("Computing driver sets that reprogram to an attractor with Apoptosis=0 . . .")
# reprogram_sets = diag.reprogram_to_trap_spaces({'Apoptosis':0})
# print("Drivers found:")
# print(reprogram_sets)
# for set_size in range(len(reprogram_sets[-1])+1):
#     n = len([x for x in reprogram_sets if len(x)==set_size])
#     print("There are",n,"driver sets of size",set_size)
#
# print("1-node driver sets for Apoptosis=0:")
# print([x for x in reprogram_sets if len(x)==1])
#
# print("2-node driver sets for Apoptosis=0:")
# print([x for x in reprogram_sets if len(x)==2])

# Showing off new mixed format
rules='''
xA*= xB | !xA & !xB
xB*= xA | !xA & !xB
xC*= xD | !xC & !xD
xD*= xC | !xC & !xD
xE*= xE & xA & xC or !(xA & xC) & ~xE
xF*= xF & xA & xC or !(xA & xC) & ~xF
'''

rules_pbn = sm.Format.booleannet2bnet(rules)
primes = PyBoolNet.FileExchange.bnet2primes(rules_pbn)
sm.Format.pretty_print_prime_rules(primes)
diag = sm.Succession.build_succession_diagram(primes)
diag.summary()
print(list(diag.digraph.edges()))
