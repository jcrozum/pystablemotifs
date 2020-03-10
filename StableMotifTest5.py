import PyBoolNet
import StableMotifs as sm

rules='''
xA* = xB
xB* = xA
xC* = xD
xD* = xC
xE* = xE or xC
'''



rules_pbn = sm.Format.booleannet2bnet(rules)
primes = PyBoolNet.FileExchange.bnet2primes(rules_pbn)
print("Source Node Example:")
sm.Format.pretty_print_prime_rules(primes)
print()
print("DIAGRAM SUMMARY (UNMERGED)")
diag_unmerge = sm.Succession.build_succession_diagram(primes, merge_equivalent_motifs=False)
diag_unmerge.summary()
print()
print()
print()
print("DIAGRAM SUMMARY (MERGED)")
diag = sm.Succession.build_succession_diagram(primes)
diag = sm.Succession.build_succession_diagram(primes)
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
print("There are", len(diag_unmerge.MotifReductionList), "network reductions in the unmerged diagram.")
print("There are", len(diag.MotifReductionList), "network reductions in the merged diagram.")
print("For example . . .")
ss = {'xA': 0, 'xB': 0, 'xC': 0, 'xD': 0, 'xE': 0}
i = 0
for x in diag_unmerge.MotifReductionList:
    if ss == x.logically_fixed_nodes:
        i += 1
print(ss,"appears",i,"times in the unmerged diagram.")
i = 0
for x in diag.MotifReductionList:
    if ss == x.logically_fixed_nodes:
        i += 1
print(ss,"appears",i,"times in the merged diagram.")