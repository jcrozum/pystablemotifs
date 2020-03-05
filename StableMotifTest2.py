import PyBoolNet
import StableMotifs as sm

rules='''
xA*= xB
xB*= xA
xC*= xA or not xD
xD*= xC
xE*= xB and xF
xF*= xE
'''

rules_pbn = sm.Format.booleannet2bnet(rules)
primes = PyBoolNet.FileExchange.bnet2primes(rules_pbn)
print("Succession Diagram Example from CSM Paper:")
sm.Format.pretty_print_prime_rules(primes)
print()
print("DIAGRAM SUMMARY")
diag = sm.Succession.build_succession_diagram(primes)
diag.summary()
print()
print("ATTRACTOR SUMMARY")
print()
diag.attractor_candidate_summary()
