import PyBoolNet
import StableMotifs as sm

primes = sm.Format.import_primes("models/test3.txt",remove_constants=True)
sm.Format.pretty_print_prime_rules(primes)

sols = sm.DomainOfInfluence.GRASP({'Apoptosis':0},primes,100)
print(sols)
