import pystablemotifs as sm
import sys

"""
Example Usage:
python FindAttractors.py ../models/PhaseSwitch.txt

The argument should be a file (path) to to a model file. The rules can be in
either BooleanNet or BNet format.

This will print the rules and a summary of the system's attractor repertoire.
"""

print("Loading network . . .")
primes = sm.format.import_primes(sys.argv[1],remove_constants=True)
print("Network loaded.")
print()
print("RULES")
sm.format.pretty_print_prime_rules({k:primes[k] for k in sorted(primes)})
print()
print("Analyzing network . . .")
ar = sm.AttractorRepertoire.from_primes(primes,max_simulate_size=7)

print("Analysis complete.")
print()
ar.summary()
