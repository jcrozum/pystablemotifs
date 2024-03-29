#check if the network reduction numbers are correct.

import pyboolnet.file_exchange
import pystablemotifs as sm

rules='''
xA* = xB
xB* = xA
xC* = xD
xD* = xC
xE* = xE or xC
'''



rules_pbn = sm.format.booleannet2bnet(rules)
primes = pyboolnet.file_exchange.bnet2primes(rules_pbn)
print("Source Node Example:")
sm.format.pretty_print_prime_rules(primes)
print()
print("DIAGRAM SUMMARY (UNMERGED)")
diag_unmerge = sm.succession.build_succession_diagram(primes, merge_equivalent_motifs=False)
diag_unmerge.summary()
print()
print()
print()
print("DIAGRAM SUMMARY (MERGED)")
diag = sm.succession.build_succession_diagram(primes)
diag = sm.succession.build_succession_diagram(primes)
diag.summary()

# The attractor summaries should both contain six attractors:
# AB motif can take two values, CD motif can take two values, if CD=00, then E motif can take two values for 2*(1+2)= 6 attractors

print()
print("ATTRACTOR SUMMARY (UNMERGED)")
print()
diag_unmerge.attractor_candidate_summary()

print()
print("ATTRACTOR SUMMARY (MERGED)")
print()
diag.attractor_candidate_summary()
print()
print("There are", len(diag_unmerge.motif_reduction_dict), "network reductions in the unmerged diagram.")
print("There are", len(diag.motif_reduction_dict), "network reductions in the merged diagram.")
print("For example . . .")
ss = {'xA': 0, 'xB': 0, 'xC': 0, 'xD': 0, 'xE': 0}
i = 0
for x in diag_unmerge.motif_reduction_dict.values():
    if ss == x.logically_fixed_nodes:
        i += 1
print(ss,"appears",i,"times in the unmerged diagram.")
i = 0
for x in diag.motif_reduction_dict.values():
    if ss == x.logically_fixed_nodes:
        i += 1
print(ss,"appears",i,"times in the merged diagram.")
