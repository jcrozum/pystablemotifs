import PyBoolNet
import PyStableMotifs as sm
from PyStableMotifs.DomainOfInfluence import logical_domain_of_influence, domain_of_influence, fixed_implies_implicant

import itertools as it

def source_sets(primes, min_set_size=None, max_set_size=None, forbidden=None, fixed=None):
    """Short summary.
    Parameters
    ----------
    primes : PyBoolNet primes dictionary
        Update rules.
    min_set_size : int
        Minimum size of source set to consider.
        If None, is set to the size of the fixed set + 1.
    max_drivers : int
        Maximum size of source set to consider.
        If None, is set to the number of all the nodes in the prime that are not in the forbidden - 1.
    forbidden : set of str variable names
        Variables to be considered uncontrollable (the default is None).
    fixed : partial state dictionaries
        Node set that should be included in all the source sets (the default is None).
    Returns
    -------
    source_sets : list of partial state dictionaries
        Each state dictionary in the list is a possible source set.
    """
    if forbidden is None:
        forbidden = set()
    if fixed is None:
        fixed = {}

    # consistency check
    for node in forbidden:
        if node in fixed.keys():
            print("The forbidden nodes are inconsistent with the fixed nodes.")
            return
        if node not in primes.keys():
            print("The forbidden nodes are inconsistent with the primes.")
            return
    for node in fixed.keys():
        if node not in primes.keys():
            print("The fixed nodes are inconsistent with the primes.")
            return

    if min_set_size is None:
        min_set_size = len(fixed) + 1
    if max_set_size is None:
        max_set_size = len(primes) - len(forbidden) - 1

    # consistency check
    if min_set_size < len(fixed):
        print("min set size is inconsistent with the fixed nodes")
        return
    if max_set_size < min_set_size:
        print("min set size is inconsistent with the max set size")
        return


    nodes = []
    source_sets = []

    for node in primes.keys():
        if node not in forbidden:
            if node not in fixed.keys():
                nodes.append(node)

    for source_set_size in range(min_set_size-len(fixed), max_set_size+1-len(fixed)):
        for add_nodes in it.combinations(nodes,source_set_size):
            for s in it.product([0,1],repeat=len(add_nodes)):
                source_set = {k:s for k,s in zip(add_nodes,s)}
                source_set.update(fixed)
                source_sets.append(source_set)

    return source_sets

if __name__ == "__main__":
    print("Loading network . . .")
    primes = sm.Format.import_primes('./models/PhaseSwitch.txt')
    print("Network loaded.")
    print()
    print("RULES")
    sm.Format.pretty_print_prime_rules({k:primes[k] for k in sorted(primes)})
    print()

    fixed_sets = source_sets(primes,min_set_size=4, max_set_size=4)
    print("there are "+str(len(fixed_sets))+" source sets")
    print(fixed_sets)

    for source_set in fixed_sets:
        print("- - - - - - - - - -")
        print("fixed: ",source_set)
        GAU_DOI, GAU_DOI_contra, GAU_unknown, GAU_unknown_contra, GAU_ar = domain_of_influence(source_set,primes,MPBN_update=False)
        GAU_ar.summary()
