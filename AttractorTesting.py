import PyBoolNet
import StableMotifs as sm

print("Loading network . . .")
primes = sm.Format.import_primes("models/test2.txt",remove_constants=True)

print("Analyzing network . . .")
ar = sm.AttractorRepertoire.AttractorRepertoire(primes)

print("Analysis complete.")
ar.summary()
