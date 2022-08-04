import pyboolnet.prime_implicants
import pystablemotifs as sm
import copy
import itertools as it
import random
from collections import namedtuple

def fixed_implies_implicant(fixed,implicant):
    """Returns True if and only if the (possibly partial) state "fixed" implies
    the implicant.

    Parameters
    ----------
    fixed : partial state dictionary
        State (or partial state) representing fixed variable states.
    implicant : partial state dictionary
        State (or partial state) representing the target implicant.

    Returns
    -------
    bool
        True if and only if the implicant is in the logical domain of influence
        of the fixed (partial) state.

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
    """Returns True if and only if the (possibly partial) state "fixed"
    contradicts the implicant.

    Parameters
    ----------
    fixed : partial state dictionary
        State (or partial state) representing fixed variable states.
    implicant : partial state dictionary
        State (or partial state) representing the target implicant.

    Returns
    -------
    bool
        True if and only if the implicant contradicts the logical domain of
        influence of the fixed (partial) state.

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
    """Computes the logical domain of influence (LDOI) (see Yang et al. 2018).
    In general, the LDOI is a subset of the full domain of influence (DOI), but it
    is much more easily (and quickly) computed.

    Parameters
    ----------
    partial_state : partial state dictionary
        pyboolnet implicant that defines fixed nodes.
    primes : pyboolnet primes dictionary
        Update rules.
    implied_hint : partial state dictionary
         Known subset of the LDOI; used during optimization.
    contradicted_hint : partial state dictionary
        Known subset of the contradiction boundary; used during optimization.

    Returns
    -------
    implied : partial state dictionary
        The logical domain of influence.
    contradicted : partial state dictionary
        The contradiction boundary.

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

def domain_of_influence(partial_state,primes,implied_hint=None,contradicted_hint=None,
    max_simulate_size=20,max_simulate_size_vc=None,max_stable_motifs=10000,max_in_degree=float('inf'),MPBN_update=False):
    """Computes the domain of influence (DOI) of the seed set. (see Yang et al. 2018)

    Parameters
    ----------
    partial_state : partial state dictionary
        pyboolnet implicant that defines fixed nodes (seed set).
    primes : pyboolnet primes dictionary
        Update rules.
    implied_hint : partial state dictionary
        Known subset of the DOI; used during optimization.
    contradicted_hint : partial state dictionary
        Known subset of the contradiction boundary; used during optimization.
    max_simulate_size : int
        Maximum number of variables for which to brute-force build a state
        transition graph (the default is 20).
    max_simulate_size_vc : int
        Maximum number of variables for which to brute-force build a state
        transition graph for the vc-reduced space (the default is the same as max_simulate_size).
    max_stable_motifs : int
        Maximum number of output lines for pyboolnet to process from the
        AspSolver (the default is 10000).
    max_in_degree : int or float
        Will not try to delete nodes that will result an increase in the
        in-degree of the downstream node so that it has in-degree larger than this.
        Deleting nodes with large in-degree can be computationally expensive (the default
        is float('inf')).
    MPBN_update : bool
        Whether MBPN update is used instead of general asynchronous update
        (see Pauleve et al. 2020)(the default is False).
    Returns
    -------
    data : namedtuple
        A namedtuple that contains the entries listed below.
    implied : partial state dictionary
        Nodes that are certain to be in the domain of influence.
    contradicted : partial state dictionary
        The contradiction boundary.
    possibly_implied : partial state dictionary
        Nodes that are possibly in the domain of influence.
    possibly_contradicted : partial state dictionary
        Nodes that are possibly in the contradiction boundary.
    attractor_repertoire : AttractorRepertoire
        The class that stores information about attractors.
    """
    # optional optimization values
    if implied_hint is None:
        implied_hint = {}
    if contradicted_hint is None:
        contradicted_hint = {}

    implied = implied_hint # the DOI
    contradicted = contradicted_hint # states implied by partial_state that contradict partial_state
    unknown = {}
    unknown_contra = {}
    primes_to_search = copy.deepcopy(primes)
    for node in implied_hint: del primes_to_search[node]
    for node in contradicted_hint: del primes_to_search[node]

    # fixed will be partial_state + its LDOI
    fixed = partial_state.copy()
    LDOI,LDOI_contra = logical_domain_of_influence(partial_state,primes)
    fixed.update(LDOI)
    implied.update(LDOI)
    contradicted.update(LDOI_contra)

    # reducing the primes by the LDOI
    primes_to_search, ps = sm.reduction.reduce_primes(fixed,primes_to_search)

    # Adding the seed nodes that are not in the LDOI as sinks in the reduced network
    sink_nodes = {}
    for node,value in fixed.items():
        if node not in LDOI:
            sink_nodes[node] = value

    # Adding back the simplifed rules for the sink_nodes
    rules_to_add = {}
    for node in sink_nodes:
        rules_to_add[node] = copy.deepcopy(primes[node])
        pyboolnet.prime_implicants.update_primes(rules_to_add,node,fixed)
    primes_to_search.update(rules_to_add)
    primes_to_search = sm.reduction.simplify_primes(primes_to_search)

    # Finding the attractors of this modified network
    if MPBN_update == True:
        attractor_dict_list = pyboolnet.trap_spaces.compute_trap_spaces(primes_to_search, "min", max_output=10000)
        sorted_attractor_dict_list = []
        for attractor in attractor_dict_list:
            for node in primes_to_search:
                if node not in attractor:
                    attractor[node] = 'X'
            sorted_attractor_dict_list.append(dict(sorted(attractor.items())))
        ar = sorted_attractor_dict_list
    else:
        ar = sm.AttractorRepertoire.from_primes(primes_to_search,
            max_simulate_size=max_simulate_size,max_simulate_size_vc=max_simulate_size_vc,
            max_stable_motifs=max_stable_motifs,max_in_degree=max_in_degree,MPBN_update=MPBN_update)
        attractor_dict_list = []
        for attractor in ar.attractors:
            attractor_dict_list.append(attractor.attractor_dict)

    # Determining which node values are shared by all attractors in the reduction
    if len(attractor_dict_list) == 0:   # there is no attractor
        print("Unable to properly count attractors.")
        return

    att = attractor_dict_list[0]  # there is at least one attractor.
    for node in att:
        value = None
        value_implied = None
        value_unknown = None
        if node in LDOI:    # if the node is in the LDOI, it is in the DOI.
            implied[node] = att[node]
            continue    # this node is in the DOI. Move on to the next.
        elif att[node] == 'X':
            continue    # this node is not in the DOI. Move on to the next.
        elif att[node] == '?' or att[node] == '!':
            value_unknown = True
        else:   # in case the node has a Boolean value
            value = att[node]

        for other in attractor_dict_list: # check the other attractors
            if node in other:
                if other[node] == 'X':
                    value_implied = False
                    break
                elif other[node] == '?' or other[node] == '!':
                    value_unknown = True
                else:   # the node in the other attractor has a Boolean value.
                    if value == None:
                        value = other[node]
                    elif int(value) == int(other[node]):
                        continue
                    else:
                        value_implied = False
                        break
            else:   # the node is not in the other attractor.
                value_implied = False
                break

        if value_implied == False:
            continue    # this node is not in the DOI, move on to the next node.
        elif value_unknown == True:
            if value == None:
                unknown[node] = 'Unknown'
            else:
                unknown[node] = int(value)
        else:
            implied[node] = int(value)

    for node in sink_nodes:
        if node in implied:
            if sink_nodes[node] != implied[node]:
                contradicted[node] = implied[node]
                del implied[node]
        elif node in unknown:
            if sink_nodes[node] != unknown[node]:
                unknown_contra[node] = unknown[node]

    doi_data = namedtuple('doi_data',
        'implied contradicted possibly_implied possibly_contradicted attractor_repertoire')
    data = doi_data(dict(sorted(implied.items())), dict(sorted(contradicted.items())),
        dict(sorted(unknown.items())), dict(sorted(unknown_contra.items())), ar)

    return data

def single_drivers(target,primes):
    """Finds all 1-node (logical) drivers of target under the rules given
    by primes.

    Parameters
    ----------
    target : partial state dictionary
        pyboolnet implicant that defines target fixed node states.
    primes : pyboolnet primes dictionary
        Update rules.

    Returns
    -------
    list of length-1 dictionaries
        Each dictionary describes a single node state that contains the target
        in its logical domain of influence.

    """
    drivers = []
    for k in primes:
        for val in [0,1]:
            ds = {k:val}
            ldoi,contra = logical_domain_of_influence(ds,primes)
            if all([kk in ldoi.items() for kk in target.items()]):
                drivers.append(ds)
    return drivers

def all_drivers_of_size(driver_set_size,target, primes, external_search_vars=None,internal_search_vars=None):
    """Finds all (logical) driver sets up to a specified size that drive the target.

    Parameters
    ----------
    driver_set_size : int
        The number of driver nodes to try to find.
    target : partial state dictionary
        pyboolnet implicant that defines target fixed node states.
    primes : pyboolnet primes dictionary
        Update rules.
    external_search_vars : set of str variable names
        Node set not in target to consider as potential drivers. If None,
        then all nodes not fixed in target (the default is None).
    internal_search_vars : set of str variable names
        Node set in target to consider as potential drivers. If None, all
        nodes in partial state (the default is None).

    Returns
    -------
    driver_sets : list of partial state dictionaries
        Each state dictionary in the list drives target.

    """

    if internal_search_vars is None:
        internal_search_vars = set(target.keys())
    if external_search_vars is None:
        external_search_vars = set(x for x in primes if not x in target)

    driver_sets = []

    search_vars = internal_search_vars | external_search_vars
    for driver_vars in it.combinations(search_vars,driver_set_size):
        # Any internal nodes must be in their target state to be drivers
        internal_driver_dict = {k:target[k] for k in driver_vars if k in target}
        out_keys = [k for k in driver_vars if not k in target]
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

            target_stabilized = True
            for k,v in target.items():
                if not k in fixed or not fixed[k] == v:
                    target_stabilized = False
                    break

            if target_stabilized:
                driver_sets.append(driver_dict)
    return driver_sets

def internal_drivers(target,primes,max_drivers=None):
    """Find internal (logical) driver nodes of target through brute-force.

    Parameters
    ----------
    target : partial state dictionary
        pyboolnet implicant that defines target fixed node states.
    primes : pyboolnet primes dictionary
        Update rules.
    max_drivers : int
        Maximum size of driver set to consider. If None, is set to the size of
        the parital state (as this is always sufficient to achieve the target)
        (the default is None).

    Returns
    -------
    driver_sets : list of partial state dictionaries
        Each state dictionary in the list drives target. These are sorted
        by length (smallest first).

    """
    if max_drivers is None:
        max_drivers = len(target) # The target itself is always its own driver

    driver_sets = []

    for driver_set_size in range(1,max_drivers+1):
        driver_sets += all_drivers_of_size(driver_set_size,target,primes,external_search_vars=set())

    if len(driver_sets) == 0:
        driver_sets.append(target)

    return sorted(driver_sets, key = lambda x: len(x))


def minimal_drivers(target,primes,max_drivers=None):
    """Finds smallest set(s) of (logical) driver nodes of target through
    brute-force. Unlike minimal_drivers, we are not limited to internal drivers
    nodes.

    Parameters
    ----------
    target : partial state dictionary
        pyboolnet implicant that defines target fixed node states.
    primes : pyboolnet primes dictionary
        Update rules.
    max_drivers : int
        Maximum size of driver set to consider. If None, is set to the size of
        the parital state (as this is always sufficient to achieve the target)
        (the default is None).

    Returns
    -------
    driver_sets : list of partial state dictionaries
        Each state dictionary in the list drives target. These are sorted
        by length (smallest first).

    """

    if max_drivers is None:
        max_drivers = len(target) # The target itself is always its own driver

    driver_sets = []
    out_vars = set(x for x in primes if not x in target)
    for driver_set_size in range(1,max_drivers+1):
        if len(driver_sets) > 0:
            break
        driver_sets += all_drivers_of_size(driver_set_size,target,primes,external_search_vars=out_vars)

    if len(driver_sets) == 0:
        driver_sets.append(target)

    return sorted(driver_sets, key = lambda x: len(x))

def knock_to_partial_state(target,primes,min_drivers=1,max_drivers=None,forbidden=None):
    """Find all partial states in primes that drive the target. Do not consider
    nodes in the forbidden list.

    Parameters
    ----------
    target : partial state dictionary
        pyboolnet implicant that defines target fixed node states.
    primes : pyboolnet primes dictionary
        Update rules.
    min_drivers : int
        Minimum size of driver set to consider. (the default is 1).
    max_drivers : int
        Maximum size of driver set to consider. If None, is set to the size of
        the parital state (as this is always sufficient to achieve the target)
        (the default is None).
    forbidden : set of str variable names
        Variables to be considered uncontrollable (the default is None).

    Returns
    -------
    knocked_nodes : list of partial state dictionaries
        Each state dictionary in the list drives target. Supersets of
        previously considered dictionaries are excluded.

    """
    if max_drivers is None:
        max_drivers = len(primes) - len(target)
    if forbidden is None:
        forbidden = set()

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

################################################################################
# Below are Greedy Random Adaptive Search Program (GRASP) methods.
# These were developed in Yang et al. 2018.
################################################################################

def _initial_GRASP_candidates(target,primes,forbidden):
    """Helper function for GRASP driver search. Constructs initial candidates
    for driver nodes.

    Parameters
    ----------
    target : partial state dictionary
        pyboolnet implicant that defines target fixed node states.
    primes : pyboolnet primes dictionary
        Update rules.
    forbidden : set of str variable names
        Variables to be considered uncontrollable (the default is None).

    Returns
    -------
    candidates : list of partial state dictionaries
        List of variable states that can potentially lead to the target.

    """
    if forbidden is None:
        candidate_vars = list(primes.keys())
    else:
        candidate_vars = [k for k in primes if not k in forbidden]
    candidates = []
    for st in [0,1]:
        candidates += [{k:st} for k in candidate_vars]
    return candidates

def _GRASP_default_scores(target,primes,candidates):
    """Helper function for GRASP driver search. Scores candidate driver nodes.

    Parameters
    ----------
    target : partial state dictionary
        pyboolnet implicant that defines target fixed node states.
    primes : pyboolnet primes dictionary
        Update rules.
    candidates : list of partial state dictionaries
        List of variable states that can potentially lead to the target.

    Returns
    -------
    scores : list of ints
        Logical domain of influence sizes for individual node states. If the
        node leads to a contradiction, the score will become the negative of the
        largest LDOI size. Scores are ordered in the same order as the
        candidates list.

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

def _construct_GRASP_solution(target,primes,candidates,scores):
    """Helper funciton for GRASP driver search. Constructs individual driver set
    using the GRASP search method.

    Parameters
    ----------
    target : partial state dictionary
        pyboolnet implicant that defines target fixed node states.
    primes : pyboolnet primes dictionary
        Update rules.
    candidates : list of partial state dictionaries
        List of variable states that can potentially lead to the target.
    scores : list of ints
        Logical domain of influence sizes for individual node states. If the
        node leads to a contradiction, the score will become the negative of the
        largest LDOI size. Scores are ordered in the same order as the
        candidates list.

    Returns
    -------
    partial state dictionary
        A partial state that contains the target in its LDOI. If no such partial
        state is found, returns an empty dictionary instead.

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

def _local_GRASP_reduction(solution,target,primes):
    """A helper funciton for GRASP driver search. Reduces valid solutions to
    attempt to remove redundancies.

    Parameters
    ----------
    solution : partial state dictionary
        Solution to be reduced; must contain the target in its LDOI.
    target : partial state dictionary
        pyboolnet implicant that defines target fixed node states.
    primes : pyboolnet primes dictionary
        Update rules.

    Returns
    -------
    partial state dictionary
        Reduced solution that also contains target in its LDOI.

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

def GRASP(target, primes, GRASP_iterations, forbidden=None, GRASP_scores=_GRASP_default_scores):
    """Search for drivers of target in primes using the method of Yang et al. 2018.

    Parameters
    ----------
    target : partial state dictionary
        pyboolnet implicant that defines target fixed node states.
    primes : pyboolnet primes dictionary
        Update rules.
    GRASP_iterations : int
        The number of times to run the GRASP method.
    forbidden : set of str variable names
        Variables to be considered uncontrollable (the default is None).
    GRASP_scores : function
        Function to score candiates (the default is _GRASP_default_scores; see
        that function for required inputs and outputs of the scoring function).

    Returns
    -------
    solutions : list of partial state dictionaries
        Each partial set dictionary represents a driver set whose LDOI contains
        the target.

    """
    solutions = []
    candidates = _initial_GRASP_candidates(target,primes,forbidden)
    scores = GRASP_scores(target,primes,candidates)
    for iter in range(GRASP_iterations):
        solution_big = _construct_GRASP_solution(target,primes,candidates,scores)
        solution = _local_GRASP_reduction(solution_big,target,primes)
        if not solution is None and len(solution) > 0 and not solution in solutions:
            solutions.append(solution)

    return solutions
