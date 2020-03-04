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

print()
print()
print()
print("============")
print("Potential \"Ghost Branch\":")
diag.MotifReductionList[0].summary()

print("============")
diag.MotifReductionList[0].build_partial_STG()
print("Partial STG edge-list in", *sorted(diag.MotifReductionList[0].reduced_primes), "order:")
print(diag.MotifReductionList[0].partial_STG.edges())
