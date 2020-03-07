import PyBoolNet
import StableMotifs as sm
import sys

rules = open(sys.argv[1]).read()

rules_pbn = sm.Format.booleannet2bnet(rules)
primes = PyBoolNet.FileExchange.bnet2primes(rules_pbn)
PyBoolNet.PrimeImplicants._percolation(primes,True) # Get rid of constants
print("TLGL Network (Fixed Inputs):")
sm.Format.pretty_print_prime_rules(primes)
print()

diag = sm.Succession.build_succession_diagram(primes)
diag.attractor_candidate_summary()
