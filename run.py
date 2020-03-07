import PyBoolNet
import StableMotifs as sm
import sys

def parse(stream):
    "Ingores commented out lines"
    lines = list(stream)
    lines = filter(lambda x: not x.startswith("#"), lines)
    rules = "\n".join(lines)
    return rules

rules = parse(open(sys.argv[1]))

# Reformat rules.
rules_pbn = sm.Format.booleannet2bnet(rules)

#
primes = PyBoolNet.FileExchange.bnet2primes(rules_pbn)

# Get rid of constants
PyBoolNet.PrimeImplicants._percolation(primes,True)

sm.Format.pretty_print_prime_rules(primes)

diag = sm.Succession.build_succession_diagram(primes)

diag.attractor_candidate_summary()
