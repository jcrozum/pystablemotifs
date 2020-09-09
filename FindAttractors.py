import PyBoolNet
import StableMotifs as sm
import sys

print("Loading network . . .")
primes = sm.Format.import_primes(sys.argv[1],remove_constants=True)
print("Network loaded.")
print()
print("RULES")
sm.Format.pretty_print_prime_rules({k:primes[k] for k in sorted(primes)})
print()
print("Analyzing network . . .")
ar = sm.AttractorRepertoire.from_primes(primes,max_simulate_size=7)

print("Analysis complete.")
print()
ar.summary()
