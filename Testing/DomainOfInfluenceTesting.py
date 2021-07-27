import PyBoolNet
import PyStableMotifs as sm
import copy

def domain_of_influence(partial_state,primes,implied_hint=None,contradicted_hint=None):
    """
    Computes the domain of influence (DOI) of the seed set.

    Parameters
    ----------
    partial_state : partial state dictionary
        PyBoolNet implicant that defines fixed nodes (seed set).
                                                                                {"node1": Boolean value, ...}
    primes : PyBoolNet primes dictionary
        Update rules.

                                                                                dictionary of lists of lists of dictionaries
                                                                                PyBoolNet Update rules.
                                                                                e.g following rules,
                                                                                A, A|B
                                                                                B, B
                                                                                are described as
                                                                                {'A':[[{'A':0,'B':0}],[{'A':1},{'B':1}]],'B':[[{'B':0}],[{'B':1}]]}

    (OPTIONAL)
    implied_hint : partial state dictionary
        Known subset of the DOI; used during optimization.
    contradicted_hint : partial state dictionary
        Known subset of the contradiction boundary; used during optimization.

    Returns
    -------
    implied : partial state dictionary
        Nodes that are certain to be in the domain of influence.
    contradicted : partial state dictionary
        The contradiction boundary.
    unknown : partial state dictionary
        Nodes that are possibly in the domain of influence.

    Calculating the LDOI and reducing the primes by this LDOI
    Adding the seed nodes that are not in the LDOI as sinks in the reduced network
    Finding the attractor repertoire of this modified network
    Determining which node values are shared by all attractors in the reduction
    Adding these common node values to the LDOI to get the DOI
    """
    # optional optimization values
    if implied_hint is None:
        implied_hint = {}
    if contradicted_hint is None:
        contradicted_hint = {}

    implied = implied_hint # the DOI
    contradicted = contradicted_hint # states implied by partial_state that contradict partial_state
    unknown = {}
    primes_to_search = copy.deepcopy(primes)
    for node in implied_hint: del primes_to_search[node]
    for node in contradicted_hint: del primes_to_search[node]

    # fixed will be partial_state + its LDOI
    fixed = partial_state.copy()
    LDOI,LDOI_contra = logical_domain_of_influence(partial_state,primes)
    fixed.update(LDOI)
    implied.update(LDOI)
    contradicted.update(LDOI_contra)

    print("updated fixed: ",fixed)

    # reducing the primes by the LDOI
    primes_to_search, ps = sm.Reduction.reduce_primes(fixed,primes_to_search)

    # Adding the seed nodes that are not in the LDOI as sinks in the reduced network
    sink_nodes = {}
    for node,value in fixed.items():
        if node not in LDOI:
            sink_nodes[node] = value

    print('sink_nodes: ',sink_nodes)
    print()

    # Adding back the simplifed rules for the sink_nodes
    rules_to_add = {}
    for node in sink_nodes:
        rules_to_add[node] = copy.deepcopy(primes[node])
    print("rules to add: ",rules_to_add)
    for node in sink_nodes:
        PyBoolNet.PrimeImplicants._substitute(rules_to_add,node,fixed)
    print("rules to add: ",rules_to_add)
    primes_to_search.update(rules_to_add)

    primes_to_search = sm.Reduction.simplify_primes(primes_to_search)

    print("MODIFIED RULES")
    sm.Format.pretty_print_prime_rules({k:primes_to_search[k] for k in sorted(primes_to_search)})
    print()

    # Finding the attractor repertoire of this modified network
    ar = sm.AttractorRepertoire.from_primes(primes_to_search)

    ar.summary()

    # Determining which node values are shared by all attractors in the reduction
    if ar.fewest_attractors == 0:   # there is no attractor
        implied.update(LDOI)
        contradicted.update(LDOI_contra)
        unknown = {}
        print("Unable to properly count attractors.")
        return implied, contradicted, unknown

    att = ar.attractors[0]  # there is at least one attractor.

    print("attractor: ", att.attractor_dict)

    for node in att.attractor_dict:
        value = None
        value_implied = None
        value_unknown = None
        if node in LDOI:    # if the node is in the LDOI, it is in the DOI.
            implied[node] = att.attractor_dict[node]
            continue    # this node is in the DOI. Move on to the next.
        elif att.attractor_dict[node] == 'X':
            continue    # this node is not in the DOI. Move on to the next.
        elif att.attractor_dict[node] == '?' or att.attractor_dict[node] == '!':
            value_unknown = True
        else:   # in case the node has a Boolean value
            value = att.attractor_dict[node]

        for other in ar.attractors: # check the other attractors
            if node in other.attractor_dict:
                if other.attractor_dict[node] == 'X':
                    value_implied = False
                    break
                elif other.attractor_dict[node] == '?' or other.attractor_dict[node] == '!':
                    value_unknown = True
                else:   # the node in the other attractor has a Boolean value.
                    if value == None:
                        value = other.attractor_dict[node]
                    elif int(value) == int(other.attractor_dict[node]):
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

    return dict(sorted(implied.items())), dict(sorted(contradicted.items())), dict(sorted(unknown.items()))

def logical_domain_of_influence(partial_state,primes,implied_hint=None,contradicted_hint=None):
    """Computes the logical domain of influence (LDOI) (see Yang et al. 2018)

    Parameters
    ----------
    partial_state : partial state dictionary
        PyBoolNet implicant that defines fixed nodes.
    primes : PyBoolNet primes dictionary
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


def main():
    print("Loading network . . .")
    primes = sm.Format.import_primes('./models/test5.txt')
    print("Network loaded.")
    print()
    print("RULES")
    sm.Format.pretty_print_prime_rules({k:primes[k] for k in sorted(primes)})
    print()

    fixed = {'ABA':1}
    print('fixed: ',fixed)

    LDOI, LDOI_contra = logical_domain_of_influence(fixed,primes)

    print('LDOI: ',LDOI)
    print('LDOI_contra: ',LDOI_contra)

    DOI, DOI_contra, unknown = domain_of_influence(fixed,primes)

    print('DOI: ',DOI)
    print('DOI_contra: ',DOI_contra)
    print('unknown: ',unknown)

    print()
    print("RULES")
    sm.Format.pretty_print_prime_rules({k:primes[k] for k in sorted(primes)})
    print()

main()
