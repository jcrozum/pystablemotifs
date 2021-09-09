#this is testing that the source nodes a properly grouped together.
#The children of the root node in the SD should be the value combinations of the input nodes
#There should be no more SMs in the SD that involve A or B

import pyboolnet.file_exchange
import pystablemotifs as sm

rules='''
xA*= xA
xB*= xB
xC*= xA and xB and xC or xD
xD*= not xC
'''



rules_pbn = sm.format.booleannet2bnet(rules)
primes = pyboolnet.file_exchange.bnet2primes(rules_pbn)
print("Source Node Example:")
sm.format.pretty_print_prime_rules(primes)
print()
print("DIAGRAM SUMMARY")
diag = sm.succession.build_succession_diagram(primes)
diag.summary()
print()
print("ATTRACTOR SUMMARY")
print()
diag.attractor_candidate_summary()
