import itertools as it
import random

def fixed_implies_implicant(fixed,implicant):
    """
    Returns True if and only if the (possibly partial) state "fixed" implies the implicant.
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

def fixed_excludes_implicant(fixed,implicant):
    """
    Returns True if and only if the (possibly partial) state "fixed" contradicts the implicant.
    """
    rval = True
    for k,v in implicant.items():
        if not k in fixed:
            continue
        elif fixed[k] != v:
            rval = False
            break
    return not rval

def logical_domain_of_influence(partial_state,primes,implied_hint=None,contradicted_hint=None):
    """
    Computes the logical domain of influence (LDOI) (see Yang et al. 2018)

    INPUTS:
    partial_state - a dict in the PyBoolNet implicant form that define fixed nodes
    primes - a PyBoolNet primes dictionary that define the update rules
    implied_hint - set of states known to be a subset of the LDOI
    contradicted_hint - set of states known to be a subset of the contradicted states

    WARNING: hint states are NOT tested for consistency, as this would defeat
    the purpose of speedup and redundancy reduction.

    OUTPUTS:
    implied - node states in the LDOI of partial_state
    contradicted - node states that are implied by a subset of the LDOI,
                   but contradict the node states specified by partial_state
    Note: implied and contradicted are dictionaries in the same format as partial_state.
    """

    if implied_hint is None:
        implied_hint = {}
    if contradicted_hint is None:
        contradicted_hint = {}

    fixed = partial_state.copy() # fixed will be partial_state + its LDOI
    fixed.update(implied_hint)
    implied = implied_hint # the LDOI
    contradicted = contradicted_hint # states implied by partial_state that contradict partial_state
    primes_to_search = primes.copy()

    for k in implied_hint: del primes_to_search[k]
    for k in contradicted_hint: del primes_to_search[k]

    while True:
        states_added = False
        deletion_list = []
        for k,v in primes_to_search.items():
            kfixed = False
            for i in [0,1]:
                if kfixed: break
                for p in v[i]:
                    if kfixed: break
                    if fixed_implies_implicant(fixed,p):
                        kfixed = True
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
    """
    Finds all logical driver sets of the specified size for the target partial state

    INPUTS:
    driver_set_size - the number of driver nodes to try to find
    partial_state - dictionary specifying the target partial state
    primes - PyBoolNet state dictionary specifying the update rules
    external_search_vars - node set not in partial_state to consider as potential
                           drivers. Default: all nodes not fixed in partial_state
    internal_search_vars - node set in partial_state to consider as potential
                           drivers. Default: all nodes in partial state

    OUTPUT:
    driver_sets - a list of state dicts, each of which drives partial_state
    """
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

def knock_to_partial_state(target,primes,min_drivers=1,max_drivers=None,forbidden=None):
    """
    Find all partial states in primes that drive the target. Do not consider nodes in
    the forbidden list.
    """
    if max_drivers is None:
        max_drivers = len(primes) - len(target)
    if forbidden is None:
        forbidden = []

    knocked_nodes = []

    internal_search = set([k for k in target if not k in forbidden])
    external_search = set([k for k in primes if not k in forbidden and not k in target])

    for n in range(min_drivers,max_drivers+1):
        found = all_drivers_of_size(n,target,primes,
            internal_search_vars=internal_search,external_search_vars=external_search)
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

def initial_GRASP_candidates(target,primes,forbidden):
    """
    Helper function for GRASP driver search.
    Constructs initial candidates for driver nodes.
    """
    if forbidden is None:
        candidate_vars = list(primes.keys())
    else:
        candidate_vars = [k for k in primes if not k in forbidden]
    candidates = []
    for st in [0,1]:
        candidates += [{k:st} for k in candidate_vars]
    return candidates

def GRASP_default_scores(target,primes,candidates):
    """
    Helper function for GRASP driver search.
    Constructs scores for candidate driver nodes.
    """
    scores = []
    m = 0 # will be the size of largest LDOI
    for candidate in candidates:
        imp,con = logical_domain_of_influence(candidate,primes)

        sc = (1+len([x for x in imp if x in target and imp[x]==target[x]])) * len(imp)
        if sc > m:
            m = sc

        if any(k in target and target[k] == con[k] for k in con):
            scores.append(-1*sc - 1)
        else:
            scores.append(sc)

    # score will be len(imp) if no contradictions in LDOI
    # but will be len(imp) - m otherwise
    scores = [x - int(x < 0)*(m+2*x+1) for x in scores]

    return scores

def construct_GRASP_solution(target,primes,candidates,scores):
    """
    A helper funciton for GRASP driver search.
    Constructs an individual driver set using the GRASP search method.
    """
    solution = {}
    alpha = random.random()

    #pass_score = alpha * (max(scores) - min(scores)) + min(scores)
    scores_pos = [x for x in scores if x > 0] + [0] # add 0 to avoid min/max errors
    pass_score = alpha * (max(scores_pos) - min(scores_pos)) + min(scores_pos)

    RCL = [x for i,x in enumerate(candidates) if scores[i] >= pass_score]

    imp_old = None
    con_old = None

    while len(RCL) > 0:
        s = random.choice(RCL)
        new_solution = solution.copy()
        new_solution.update(s)

        imp,con = logical_domain_of_influence(new_solution,primes,implied_hint=imp_old,contradicted_hint=con_old)

        if any(k in target and target[k] != con[k] for k in con):
            RCL = [x for x in RCL if x != s]
            continue

        if target.items() <= imp.items():
            return new_solution

        imp_old = imp
        con_old = con
        solution = new_solution
        RCL = [x for x in RCL if not next(iter(x.keys())) in solution]

    return {}

def local_GRASP_reduction(solution,target,primes):
    """
    A helper funciton for GRASP driver search.
    Reduces valid solutions to attempt to remove redundancies.
    """
    keylist = list(solution.keys())
    random.shuffle(keylist)
    old_solution = {x:solution[x] for x in keylist}

    for k in keylist:
        new_solution = {x:v for x,v in old_solution.items() if x != k}
        if len(new_solution) == 0:
            return old_solution

        imp,con = logical_domain_of_influence(new_solution, primes)
        if target.items() <= imp.items():
            old_solution = new_solution.copy()

    return old_solution

def GRASP(target, primes, GRASP_iterations, forbidden=None, GRASP_scores=GRASP_default_scores):
    """
    Search for drivers of target in primes.

    INPUTS:
    target - dict representing the target partial state
    primes - PyBoolNet primes dict describing the update rules
    GRASP_iterations - the number of times to run the GRASP method
    forbidden - a list of variables that should not be considered as potential drivers
    GRASP_scores - a function for calculating heuristic scores of candidates;
                   takes arguments target,primes,candidates (in that order).
                   Must return a list of neumeric scores in the same order as
                   candidates. Candidates with negative scores will be removed
                   from consideration.

    OUTPUT:
    solutions - A list of driver sets in dictionary form
    """
    solutions = []
    candidates = initial_GRASP_candidates(target,primes,forbidden)
    scores = GRASP_scores(target,primes,candidates)
    for iter in range(GRASP_iterations):
        solution_big = construct_GRASP_solution(target,primes,candidates,scores)
        solution = local_GRASP_reduction(solution_big,target,primes)
        if not solution is None and len(solution) > 0 and not solution in solutions:
            solutions.append(solution)

    return solutions
