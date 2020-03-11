import PyBoolNet
import StableMotifs as sm
import StableMotifs.DomainOfInfluence as di

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
drivers = di.find_internal_motif_drivers(motif,primes)
print("The following sets are internal driver sets of the stable motif", motif)
print(drivers)
