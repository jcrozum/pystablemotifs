import PyBoolNet
import StableMotifs as sm

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

diag = sm.Succession.build_succession_diagram(primes)
diag.summary()
