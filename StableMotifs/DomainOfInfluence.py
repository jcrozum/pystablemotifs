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

import random

def initial_GRASP_candidates(target,primes):
    candidate_vars = [k for k in primes if not k in target]
    candidates = []
    for st in [0,1]:
        candidates += [{k:st} for k in candidate_vars]
    return candidates

def GRASP_scores(target,primes,candidates):
    scores = []
    m = 0 # will be the size of largest LDOI
    for candidate in candidates:
        imp,con = logical_domain_of_influence(candidate,primes)
        sc = len(imp)
        if sc > m:
            m = sc

        if any([k in target and target[k] == con[k] for k in con]):
            scores.append(-1*sc - 1)
        else:
            scores.append(sc)

    # score will be len(imp) if no contradictions in LDOI
    # but will be len(imp) - m otherwise
    scores = [x - int(x < 0)*(m+2*x+1) for x in scores]

    return scores

def construct_GRASP_solution(target,primes):
    solution = {}
    alpha = random.random()

    candidates = initial_GRASP_candidates(target,primes)
    scores = GRASP_scores(target,primes,candidates)
    pass_score = alpha * (max(scores) - min(scores)) + min(scores)
    RCL = [x for i,x in enumerate(candidates) if scores[i] >= pass_score]

    while len(RCL) > 0:
        s = random.choice(RCL)
        new_solution = solution.copy()
        new_solution.update(s)

        imp,con = logical_domain_of_influence(new_solution,primes)

        if any([k in target and target[k] != con[k] for k in con]):
            RCL = [x for x in RCL if x != s]
            continue

        if target.items() <= imp.items():
            return new_solution

        solution = new_solution
        RCL = [x for x in RCL if not next(iter(x.keys())) in solution]

    return {}

def local_GRASP_reduction(solution,target,primes):
    old_solution = solution.copy()
    for k in solution:
        new_solution = {x:v for x,v in old_solution.items() if x != k}
        if len(new_solution) == 0:
            return old_solution

        imp,con = logical_domain_of_influence(new_solution, primes)
        if target.items() <= imp.items():
            old_solution = new_solution

    return old_solution

def GRASP(target, primes, max_iterations):
    solutions = []
    for iter in range(max_iterations):
        solution = construct_GRASP_solution(target,primes)
        solution = local_GRASP_reduction(solution,target,primes)
        if not solution is None and len(solution) > 0 and not solution in solutions:
            solutions.append(solution)

    return solutions
