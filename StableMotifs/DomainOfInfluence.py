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

def single_drivers(partial_state,primes):
    """
    Finds all 1-node (logical) drivers of partial_state under the rules given by primes

    Returns a list of length-1 dictionaries
    """
    drivers = []
    for k in primes:
        for val in [0,1]:
            ds = {k:val}
            ldoi,contra = logical_domain_of_influence(ds,primes)
            if all([kk in ldoi for kk in partial_state]):
                drivers.append(ds)
    return drivers

def all_drivers_of_size(driver_set_size,partial_state, primes, external_search_vars=None,internal_search_vars=None):
    if internal_search_vars is None:
        internal_search_vars = set(partial_state.keys())
    if external_search_vars is None:
        external_search_vars = set(x for x in primes if not x in partial_state)

    driver_sets = []

    search_vars = internal_search_vars | external_search_vars
    for driver_vars in it.combinations(search_vars,driver_set_size):
        # Any internal nodes must be in their partial_state state to be drivers
        internal_driver_dict = {k:partial_state[k] for k in driver_vars if k in partial_state}
        out_keys = [k for k in driver_vars if not k in partial_state]
        for s in it.product([0,1],repeat=len(out_keys)):
            external_driver_dict = {k:s for k,s in zip(out_keys,s)}
            driver_dict = {**internal_driver_dict, **external_driver_dict}
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

            partial_state_stabilized = True
            for k,v in partial_state.items():
                if not k in fixed or not fixed[k] == v:
                    partial_state_stabilized = False
                    break

            if partial_state_stabilized:
                driver_sets.append(driver_dict)
    return driver_sets

def internal_drivers(partial_state,primes,max_drivers=None):
    """
    Find internal (logical) driver nodes of partial_state through brute-force
    """
    if max_drivers is None:
        max_drivers = len(partial_state) # The partial_state itself is always its own driver

    driver_sets = []

    for driver_set_size in range(1,max_drivers+1):
        driver_sets += all_drivers_of_size(driver_set_size,partial_state,primes,external_search_vars=set())

    if len(driver_sets) == 0:
        driver_sets.append(partial_state)

    return sorted(driver_sets, key = lambda x: len(x))


def minimal_drivers(partial_state,primes,max_drivers=None):
    """
    Finds smallest set(s) of (logical) driver nodes of partial_state through brute-force
    Not limited to internal drivers nodes
    """
    if max_drivers is None:
        max_drivers = len(partial_state) # The partial_state itself is always its own driver

    driver_sets = []
    out_vars = set(x for x in primes if not x in partial_state)
    for driver_set_size in range(1,max_drivers+1):
        if len(driver_sets) > 0:
            break
        driver_sets += all_drivers_of_size(driver_set_size,partial_state,primes,external_search_vars=out_vars)

    if len(driver_sets) == 0:
        driver_sets.append(partial_state)

    return sorted(driver_sets, key = lambda x: len(x))

def knock_to_partial_state(target,primes,min_drivers=1,max_drivers=None):
    if max_drivers is None:
        max_drivers = len(primes) - len(target)
    knocked_nodes = []

    for n in range(min_drivers,max_drivers+1):
        found = all_drivers_of_size(n,target, primes)
        for x in found:
            add_x = True
            for y in knocked_nodes:
                if x.items() >= y.items():
                    add_x = False
                    break
            if add_x:
                knocked_nodes.append(x)

        n += 1
    return knocked_nodes
