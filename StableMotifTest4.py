import PyBoolNet
import StableMotifs as sm

rules='''
xA*= xA
xB*= xB
xC*= xA and xB and xC or xD
xD*= not xC
'''



rules_pbn = sm.Format.booleannet2bnet(rules)
primes = PyBoolNet.FileExchange.bnet2primes(rules_pbn)
print("Source Node Example:")
sm.Format.pretty_print_prime_rules(primes)
print()
print("DIAGRAM SUMMARY")
diag = sm.Succession.build_succession_diagram(primes)
diag.summary()
print()
print("ATTRACTOR SUMMARY")
print()
diag.attractor_candidate_summary()
