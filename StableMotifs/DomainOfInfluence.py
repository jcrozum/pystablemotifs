import itertools as it

def fixed_implies_implicant(fixed,implicant):
    """
    Returns True iff the partial state "fixed" implies the implicant
    """
    rval = True
    for k,v in implicant.items():
        if not k in fixed:
            rval = False
            break
        elif fixed[k] != v:
            rval = False
            break
    return rval


def logical_domain_of_influence(state,primes):
    """
    Computes the logical domain of influence (LDOI) (see Yang et al. 2018)

    Inputs:
    state - a dict in the PyBoolNet implicant form that define fixed nodes
    primes - a PyBoolNet primes dictionary that define the update rules

    Outputs:
    implied - node states in the LDOI of state
    contradicted - node states that are implied by a subset of the LDOI,
                   but contradict the node states specified by state
    Note: implied and contradicted are dictionaries in the same format as state.
    """
    fixed = state.copy()
    implied = {}
    contradicted = {}
    primes_to_search = primes.copy()

    while True:
        states_added = False
        deletion_list = []
        for k,v in primes_to_search.items():
            for i in [0,1]:
                for p in v[i]:
                    if fixed_implies_implicant(fixed,p):
                        deletion_list.append(k)
                        states_added = True
                        if k in fixed:
                            if fixed[k] == i: implied[k] = i
                            else: contradicted[k] = i
                        else:
                            implied[k] = i
                            fixed[k] = i
        for k in set(deletion_list): del primes_to_search[k]
        if not states_added or len(primes_to_search) == 0: break
    return implied, contradicted

def single_drivers(ts,primes):
    """
    Finds all 1-node (logical) drivers of ts under the rules given by primes
    """
    drivers = []
    for k in primes:
        for val in [0,1]:
            ds = {k:val}
            ldoi,contra = logical_domain_of_influence(ds,primes)
            if all([kk in ldoi for kk in ts]):
                drivers.append(ds)
    return drivers

def find_internal_motif_drivers(motif,primes,max_drivers=None):
    """
    Find internal driver nodes of motif through brute-force
    """
    if max_drivers is None:
        max_drivers = len(motif)

    driver_sets = []

    for driver_set_size in range(max_drivers):
        for driver_vars in it.combinations(motif.keys(),driver_set_size):
            driver_dict = {k:motif[k] for k in driver_vars}

            is_candidate=True
            # The driver set we're looking at has too many nodes if we already
            # found a driver set that is a subset of it.
            for known_driver_set in driver_sets:
                if driver_dict.items() >= known_driver_set.items():
                    is_candidate = False
                    break
            if not is_candidate:
                continue

            implied,contradicted = logical_domain_of_influence(driver_dict,primes)
            fixed = {**implied, **driver_dict}

            motif_stabilized = True
            for k,v in motif.items():
                if not k in fixed or not fixed[k] == v:
                    motif_stabilized = False
                    break

            if motif_stabilized:
                driver_sets.append(driver_dict)

    if len(driver_sets) == 0:
        driver_sets.append(motif)

    return sorted(driver_sets, key = lambda x: len(x))
