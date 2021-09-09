import pystablemotifs as sm
from pystablemotifs.drivers import logical_domain_of_influence, domain_of_influence, fixed_implies_implicant
import itertools as it
import time

def source_sets(primes, min_set_size=None, max_set_size=None, forbidden=None, fixed=None):
    """Short summary.
    Parameters
    ----------
    primes : pyboolnet primes dictionary
        Update rules.
    min_set_size : int
        Minimum size of source set to consider.
        If None, is set to the size of the fixed set + 1.
    max_drivers : int
        Maximum size of source set to consider.
        If None, is set to the number of all the nodes in the prime that are not in the forbidden - 1.
    forbidden : list of str variable names
        Variables to be considered uncontrollable (the default is None).
    fixed : partial state dictionaries
        Node set that should be included in all the source sets (the default is None).
    Returns
    -------
    source_sets : list of partial state dictionaries
        Each state dictionary in the list is a possible source set.
    """
    if forbidden is None:
        forbidden = []
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
                source_set = dict(sorted(source_set.items()))
                source_sets.append(source_set)

    return source_sets

if __name__ == "__main__":
    print("Loading network . . .")
    primes = sm.format.import_primes('./models/ABA_49_reduced.txt')
    print("Network loaded.")
    print()
    print("RULES")
    sm.format.pretty_print_prime_rules({k:primes[k] for k in sorted(primes)})
    print()

    fixed_sets = source_sets(primes,min_set_size=2,max_set_size=2,forbidden=['Closure'], fixed = {'ABA':0})
    # fixed_sets = [{'ABA':0}, {'ABA':1}]
    print("there are "+str(len(fixed_sets))+" source sets")
    print(fixed_sets)

    for source_set in fixed_sets:
        start_time = time.time()
        print("- - - - - - - - - -")
        print("fixed: ",source_set)
        LDOI, LDOI_contra = logical_domain_of_influence(source_set,primes)
        GAU_DOI, GAU_DOI_contra, GAU_unknown, GAU_unknown_contra, GAU_ar = domain_of_influence(source_set,primes,max_simulate_size=20,MPBN_update=False)
        MPBN_DOI, MPBN_DOI_contra, MPBN_unknown, MPBN_unknown_contra, MPBN_ar = domain_of_influence(source_set,primes,MPBN_update=True)

        set1 = set(LDOI.items())
        set2 = set(GAU_DOI.items())
        set3 = set(MPBN_DOI.items())

        # check the difference between DOI and LDOI in GAU
        print('LDOI: ',LDOI)
        print('LDOI_contra: ',LDOI_contra)

        print()
        print('GAU_DOI: ',GAU_DOI)
        print('GAU_DOI_contra: ',GAU_DOI_contra)
        print('GAU_unknown: ',GAU_unknown)
        print('GAU_unknown_contra: ',GAU_unknown_contra)
        if not len(set1 ^ set2) == 0:
            print('difference between LDOI and GAU_DOI:', sorted(set1 ^ set2))
        elif len(set1 ^ set2) == 0:
            print("GAU_DOI is the same as LDOI")

        print()
        print('MPBN_DOI: ',MPBN_DOI)
        print('MPBN_DOI_contra: ',MPBN_DOI_contra)
        print('MPBN_unknown: ',MPBN_unknown)
        print('MPBN_unknown_contra: ',MPBN_unknown_contra)
        if not len(set1 ^ set3) == 0:
            print('difference between LDOI and MPBN_DOI:', sorted(set1 ^ set3))
        elif len(set1 ^ set3) == 0:
            print("MPBN_DOI is the same as LDOI")
        # check the difference between GAU and MPBN
        if not len(set2 ^ set3) == 0:
            print()
            print('difference between GAU and MPBN:', sorted(set2 ^ set3))

            check1 = False
            if len(set(GAU_unknown.items())-set3) == 0:
                check1 = True
            print('is GAU_unknown included in MPBN_DOI?: ', check1)
            if check1 == False:
                print('Only in the GAU_unknown:', sorted(set(GAU_unknown.items())-set3))

            check2 = False
            if len(set2-set3) == 0:
                check2 = True
            print('is GAU_DOI included in MPBN_DOI?: ', check2)
            if check2 == False:
                print('Only in the GAU_DOI:', sorted(set2-set3))

            print('Only in the MPBN_DOI:', sorted(set3-set2-set(GAU_unknown.items())))
        elif len(set2 ^ set3) == 0:
            print("There is no difference between GAU and MPBN")
        print()
        print("Attractors for GAU")
        GAU_ar.summary()
        print("Attractors for MPBN")
        if len(MPBN_ar) == 1:
            print("There is 1 attractor.")
        else:
            print("There are " + str(len(MPBN_ar)) + " attractors.")
        for attractor in MPBN_ar:
            print(attractor)
            print()
        print("--- It took %s seconds to analyze ---" % (time.time() - start_time))
