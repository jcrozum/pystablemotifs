import pyboolnet.trap_spaces
import pyboolnet.prime_implicants
import pyboolnet.interaction_graphs
import pyboolnet.digraphs
from pyboolnet.external.bnet2primes import bnet_text2primes
import itertools as it
import networkx as nx
import re
import sympy

import pystablemotifs.time_reversal as sm_time
import pystablemotifs.restrict_space as sm_rspace
import pystablemotifs.format as sm_format
import pystablemotifs.drivers as sm_doi

def simplify_primes(primes):
    """Simplifies pyboolnet primes (e.g., A | A & B becomes A)

    Parameters
    ----------
    primes : pyboolnet primes dictionary
        Rules to simplify.

    Returns
    -------
    pyboolnet primes dictionary
        Simplified rules.

    """

    # reimport to force simplification
    if len(primes) > 0:
        return bnet_text2primes(sm_format.primes2bnet(primes))
    else:
        return primes

def reduce_primes(fixed,primes):
    """Simplifies boolean rules when some nodes are held fixed

    Parameters
    ----------
    fixed : partial state dictionary
        Node states to be held fixed.
    primes : pyboolnet primes dictionary
        Update rules.

    Returns
    -------
    reduced_primes : pyboolnet primes dictionary
        Simplified update rules
    percolated_states : partial state dictionary
        Fixed node states (including inputs) that were simplified and removed.

    """
    reduced_primes = pyboolnet.prime_implicants.percolate(primes, add_constants=fixed, remove_constants=False, copy=True)
    percolated_states = pyboolnet.prime_implicants.find_constants(reduced_primes)
    pyboolnet.prime_implicants.remove_all_constants(reduced_primes, copy=False) # remove constants in place

    return simplify_primes(reduced_primes), percolated_states

def delete_node(primes, node):
    """Reduces Boolean rules given by primes by deleting the variable specified by
    node. The deleted node may not appear in its own update function. Any update
    rules depending on the deleted node will have that dependence replaced by
    the update function of the deleted node. The rules are simplified after node
    deletion.

    Parameters
    ----------
    primes : pyboolnet primes dictionary
        Update rules.
    node : str
        Name of the node to delete.

    Returns
    -------
    new_primes : pyboolnet primes dictionary
        The reduced primes.
    constants : partial state dictionary
        Node states that became logically fixed during simplification.


    """
    succ = set()
    for p in primes:
        p_is_succ = False
        for v in [0,1]:
            for term in primes[p][v]:
                if node in term:
                    succ.add(p)
                    p_is_succ = True
                    break
            if p_is_succ:
                break


    # uncomment for additional (and hopefully unnecessary), but slow safety checks.
    #G = pyboolnet.interaction_graphs.primes2igraph(primes)
    #assert not G.has_edge(node,node), ' '.join(["Node",str(node),"has a self-loop and cannot be deleted."])

    assert not node in succ, ' '.join(["Node",str(node),"has a self-loop and cannot be deleted."])
    new_primes = {k:v for k,v in primes.items() if not k == node}
    constants = {}

    neg = "!"+node

    expr0 = sm_format.rule2bnet(primes[node][0])
    expr1 = sm_format.rule2bnet(primes[node][1])

    for child in succ:#G.successors(node):
        # If we have already simplified this child node, skip it
        if child in constants:
            continue
        crule1 = sm_format.rule2bnet(primes[child][1])
        crule1 = simplify_using_expression_and_negation(node,expr0,expr1,crule1)
        crule0 = sm_format.rule2bnet(primes[child][0])
        crule0 = simplify_using_expression_and_negation(node,expr0,expr1,crule0)
        new_primes[child] = sm_format._build_rule_using_bnet_dnfs(crule0,crule1)

        pyboolnet.prime_implicants.percolate(new_primes, remove_constants=False, copy=False)
        constants.update(pyboolnet.prime_implicants.find_constants(new_primes))
        #pyboolnet.prime_implicants.remove_all_constants(new_primes, copy=False)

    pyboolnet.prime_implicants.percolate(new_primes, remove_constants=False, copy=False)
    constants.update(pyboolnet.prime_implicants.find_constants(new_primes))
    pyboolnet.prime_implicants.remove_all_constants(new_primes, copy=False)
    new_primes = simplify_primes(new_primes)
    return new_primes, constants

def simplify_using_expression_and_negation(node,expr0,expr1,bnet):
    """Simplify the expression bnet by substituting the value for node given by
    node = expr1 = !expr0 (does not check that expr1=!expr0)

    Parameters
    ----------
    node : str
        Name of node to substitute
    expr0 : str
        Expression to substitute for !node
    expr1 : str
        Expression to substitute for node
    bnet : str
        BNET expression in which to perform the substitutions.

    Returns
    -------
    str
        Simplified BNET expression after substitutions are performed.

    """

    neg = "!"+node
    crule = re.sub(rf'\b{neg}\b',"("+expr0+")",bnet)
    crule = re.sub(rf'\b{node}\b',"("+expr1+")",crule)
    crule = sm_format.bnet2sympy(crule)
    crule = str(sympy.to_dnf(sympy.simplify(sympy.sympify(crule))))
    crule = sm_format.sympy2bnet(crule)
    return crule

def remove_outdag(primes):
    """Removes the terminal directed acyclic part of the regulatory network. This
    part of the network does not influence the attractor repertoire.

    Parameters
    ----------
    primes : pyboolnet primes dictionary
        Update rules.

    Returns
    -------
    reduced : pyboolnet primes dictionary
        The reduced primes.
    constants : partial state dictionary
        Node states that became logically fixed during reduction.

    """
    if len(primes) == 0:
        return primes, {}
    G = pyboolnet.interaction_graphs.primes2igraph(primes)
    od = pyboolnet.digraphs.find_outdag(G)
    reduced = primes.copy()
    constants = {}
    for node in od:
        if node in reduced:
            reduced, nc = delete_node(reduced, node)
            constants.update(nc)
    return reduced, constants

def deletion_reduction(primes, max_in_degree = float('inf')):
    """Implements the reduction method of Veliz-Cuba (2011).
    Deletion order is such that nodes with low in-degree are prioritized for
    removal. Deletion proceeds until all remaining nodes have self-loops.

    Parameters
    ----------
    primes : pyboolnet primes dictionary
        Update rules.
    max_in_degree : int or float
        Will not try to delete nodes that will result an increase in the
        in-degree of the downstream node so that it has in-degree larger than this.
        Deleting nodes with large in-degree can be computationally expensive (the default
        is float('inf')).

    Returns
    -------
    reduced : pyboolnet primes dictionary
        The reduced primes.
    constants : partial state dictionary
        Node states that became logically fixed during reduction.

    """

    reduced, constants = remove_outdag(primes)
    G = pyboolnet.interaction_graphs.primes2igraph(reduced)
    cur_order = sorted(reduced,key=lambda x: G.in_degree(x))

    change = True
    while change and len(reduced) > 0:
        change = False
        for node in cur_order:
            retry_node = True
            if not node in reduced:
                continue
            stop = False
            if not max_in_degree == float('inf') and G.in_degree(node) > 1:
                ups_nodes = set(G.predecessors(node))
                dns_nodes = list(G.successors(node))
                for dns in dns_nodes:
                    regulators = set(G.predecessors(dns))
                    new_regulators = regulators | ups_nodes
                    n = len(new_regulators)
                    if n > len(regulators) and n > max_in_degree:
                        stop = True
                        break
            if stop:
                continue
            if not any(node in p for p in reduced[node][1]):
                reduced, nc = delete_node(reduced, node)
                constants.update(nc)
                if len(reduced) > 0:
                    G = pyboolnet.interaction_graphs.primes2igraph(reduced)
                else:
                    G = nx.DiGraph()
                change = True
                # print("\n\n\nNEW PRIMES")
                # sm_format.pretty_print_prime_rules(reduced)
                break
        cur_order = sorted(reduced,key=lambda x: G.in_degree(x))

        pyboolnet.prime_implicants.percolate(reduced, add_constants=constants, remove_constants=False, copy=False)
        constants.update(pyboolnet.prime_implicants.find_constants(reduced))
        pyboolnet.prime_implicants.remove_all_constants(reduced, copy=False)
    return reduced, constants

def mediator_reduction(primes):
    """Network reduction method of Saadadtpour, Albert, Reluga (2013)
    Preserves fixed points. Number of complex attractors is often, but
    not always conserved (despite inital claims). Can be viewed as a more
    restrictive version of the deletion reduction method of Veliz-Cuba (2011).

    Parameters
    ----------
    primes : pyboolnet primes dictionary
        Update rules.

    Returns
    -------
    reduced : pyboolnet primes dictionary
        The reduced primes.
    constants : partial state dictionary
        Node states that became logically fixed during reduction.

    """
    return primes, {}
    reduced, constants = remove_outdag(primes)
    cur_order = sorted(reduced)
    G = pyboolnet.interaction_graphs.primes2igraph(reduced)
    candidates = [v for v in reduced if G.in_degree(v) == G.out_degree(v) == 1 and not G.has_edge(v,v)]

    for node in candidates:
        u = list(G.predecessors(node))[0]
        w = list(G.successors(node))[0]
        if not w in G.successors(u) and not w in G.predecessors(u):
            reduced, nc = delete_node(reduced, node)
            constants.update(nc)
            G = pyboolnet.interaction_graphs.primes2igraph(reduced)
            candidates = [v for v in reduced if G.in_degree(v) == G.out_degree(v) == 1 and not G.has_edge(v,v)]

    return reduced, constants


class MotifReduction:
    """Class to generate and store data about a network reduction that arises
    during the stable motif succession diagram construction algorithm.

    Parameters
    ----------
    motif_history : list of partial state dictionaries
        Stable motifs that can lock in to give the reduced network (in order).
    fixed : partial state dictionary
        Nodes values that have been fixed and reduced by stable motifs and their
        logical domain of influence.
    reduced_primes : pyboolnet primes dictionary
        Update rules for the reduced network.
    max_simulate_size : int
        Maximum number of variables for which to brute-force build a state
        transition graph (the default is 20).
    max_simulate_size_vc : int
        Maximum number of variables for which to brute-force build a state
        transition graph for the vc-reduced space (the default is the same as max_simulate_size).
    prioritize_source_motifs : bool
        Whether source nodes should be considered first (the default is True).
    max_stable_motifs : int
        Maximum number of output lines for pyboolnet to process from the
        AspSolver (the default is 10000).
    max_in_degree : int or float
        Will not try to delete nodes that will result an increase in the
        in-degree of the downstream node so that it has in-degree larger than this.
        Deleting nodes with large in-degree can be computationally expensive (the default
        is float('inf')).
    MPBN_update : bool
        Whether MPBN update is used instead of general asynchronous update
        (see Pauleve et al. 2020)(the default is False).

    Attributes
    ----------
    merged_history_permutations : list of lists of int
        Permutations of motif_history (by index) that also yield this reduction.
    logically_fixed_nodes : partial state dictionary
        Nodes values that have been fixed and reduced by stable motifs and their
        logical domain of influence.
    time_reverse_primes : pyboolnet primes dictionary
        Update rules of the time reversed reduced system.
    stable_motifs : list of partial state dictionaries
        Stable motifs of the reduced system.
    time_reverse_stable_motifs : list of partial state dictionaries
        Stable motifs of the time reversed system.
    merged_source_motifs : list of partial state dictionaries
        List of source-like stable motifs that have been merged into a single
        motif to avoid redundancy.
    source_independent_motifs : list of partial state dictionaries
        Stable motifs that exist independent of the values of the source nodes
    merge_source_motifs : list of partial state dictionaries
        Stable motifs generated by merging the stable motifs corresponding to
        source nodes.
    rspace : rspace list
        The rspace, or "restrict space" of the reduced network, describing a
        necessary condition for the system to avoid activating additional
        stable motifs (see restrict_space.py for further details).
    motif_history : list of partial state dictionaries
        Stable motifs that can lock in to give the reduced network (in order)
    reduced_primes : pyboolnet primes dictionary
        Update rules for the reduced network.
    fixed_rspace_nodes : partial state dictionary
        Nodes values that are fixed in the rspace.
    rspace_constraint : str
        BNET expression that is true in and only in the rspace.
    reduced_rspace_constraint : str
        S simplification of the rspace_constraint given the fixed_rspace_nodes
        states are satisfied
    rspace_update_primes : pyboolnet primes dictionary
        The update rules obtained from simplifying under the assumption that the
        fixed_rspace_nodes are fixed
    conserved_functions : list of pyboolnet expressions
        Boolean functions that are constant within every attractor, in pyboolnet
        update rule format
    rspace_attractor_candidates : list of str
        Attractors (lists of statestrings) in the rspace_update_primes that
        satisfy the reduced_rspace_constraint
    partial_STG : networkx.DiGraph
        Subgraph of the state transition graph of the reduced network that
        contains any and all attractors that do not lie in any of the reduced
        network's stable motifs.
    no_motif_attractors : list of str
        Complex attractors that do not "lock in" any additional stable motifs,
        stored as collections of state strings.
    attractor_dict_list : list of dictionaries
        Dictionaries corresponding to attractors that are in this reductions, but
        not in any of its subreductions (if it has any). Each describes the node
        states in the attractors according to the following
        1 variable is "ON"
        0 variable is "OFF"
        X variable is known to oscillate
        ? at least one such variable must oscillate
        ! the attractor may be false; if it is genuine, at least
          one such variable must oscillate
    terminal : str
        One of "yes", "no", or "possible", indicating whether the reduction contains
        attractors that are not in any of its subreductions.
    delprimes : pyboolnet prime dictionary
        Update rules for the system's deletion projection. Steady states and stable
        motif activation are preserved. These rules may yield additional, spurious
        complex attractors.
    deletion_STG : networkx.DiGraph
        Portion of the deletion projection's STG that contains all motif-avoidant
        attractors.
    deletion_no_motif_attractors : list of str
        Motif avoidant attractors of the deletion projection. The number of these is
        an upper bound on the number of motif avoidant attractors in the reduction.
    """

    def __init__(self,motif_history,fixed,reduced_primes,max_simulate_size=20,max_simulate_size_vc=None,prioritize_source_motifs=True,max_stable_motifs=10000,max_in_degree=float('inf'),MPBN_update=False):
        if max_simulate_size_vc == None:
            max_simulate_size_vc = max_simulate_size

        if motif_history is None:
            self.motif_history = []
        else:
            self.motif_history = motif_history.copy()
        self.merged_history_permutations = []
        self.logically_fixed_nodes = fixed
        self.reduced_primes = reduced_primes.copy()
        self.stable_motifs = pyboolnet.trap_spaces.compute_trap_spaces(self.reduced_primes, "max",max_output=max_stable_motifs)
        if MPBN_update==False:
            self.time_reverse_primes =sm_time.time_reverse_primes(self.reduced_primes)
            self.time_reverse_stable_motifs = pyboolnet.trap_spaces.compute_trap_spaces(self.time_reverse_primes, "max",max_output=max_stable_motifs)
        self.merged_source_motifs=None
        self.source_independent_motifs=None
        if self.motif_history == [] and prioritize_source_motifs:
            if MPBN_update==False:
                self.merge_source_motifs()
            else:
                self.simple_merge_source_motifs(reduced_primes,MPBN_update=MPBN_update)

        # These may or may not get calculated.
        # Sensible default values are in comments, but we will just use None for now.
        self.fixed_rspace_nodes=None # {}
        self.rspace_constraint=None # ""
        self.reduced_rspace_constraint=None # ""
        self.rspace_update_primes=None # {}
        self.conserved_functions=None # [[{}]]
        self.rspace_attractor_candidates=None # []
        self.partial_STG=None # nx.DiGraph()
        self.no_motif_attractors=None # []
        self.deletion_STG = None
        self.deletion_no_motif_attractors = None
        self.attractor_constants = None
        self.attractor_dict_list=None # []
        # self.attractor_dict={}

        # skips finding motif avoidant attractors when using MPBN_update
        if MPBN_update:
            if len(self.stable_motifs) == 0: # terminal reduction iff no more stable motifs
                self.terminal = "yes"
            else:
                self.terminal = "no"
            self.attractor_dict_list = self.generate_attr_dict_list(MPBN_update=MPBN_update)
            return

        self.rspace = sm_rspace.rspace(self.stable_motifs, self.time_reverse_stable_motifs,self.reduced_primes)

        study_possible_oscillation = False

        if not self.merged_source_motifs is None:
            self.terminal = "no"
        elif self.rspace == [[{'0':1}]] and len(self.stable_motifs) > 0: # a stable motif must lock in
            self.terminal = "no"
        elif self.rspace == [[{}]] and len(self.stable_motifs) > 0: # could not find 1-node drivers
            self.terminal = "possible"
            study_possible_oscillation = True
        elif len(self.stable_motifs) == 0: # necessarily terminal
            self.terminal = "yes"
            if len(self.reduced_primes) > 0: # Terminates in oscillation, else, fixed point
                study_possible_oscillation = True
        else: # found 1-node drivers, so we can investigate further
            self.terminal = "possible" # TODO: implement case-checking based on rspace

            #self.fixed_rspace_nodes =sm_rspace.fixed_rspace_nodes(self.rspace,self.reduced_primes)
            self.rspace = sm_rspace.reduce_rspace(self.rspace,self.reduced_primes)
            self.fixed_rspace_nodes = self.rspace[0][0]

            for motif in self.stable_motifs:
                if motif.items() <= self.fixed_rspace_nodes.items():
                    self.terminal = "no"
                    break
            if self.terminal == "possible":
                self.rspace_constraint = sm_format.pretty_print_rspace(self.rspace,simplify=False)
                #self.reduced_rspace_constraint = sm_rspace.reduce_rspace_string(self.rspace_constraint,self.fixed_rspace_nodes,simplify=False)
                self.reduced_rspace_constraint=sm_format.pretty_print_rspace(self.rspace[1:],simplify=False)
                self.rspace_update_primes = reduce_primes(self.fixed_rspace_nodes,self.reduced_primes)[0]
                #self.test_rspace(search_partial_STGs = search_partial_STGs)
            study_possible_oscillation = self.terminal == "possible" # value may be changed by test_rspace

        if study_possible_oscillation and max_simulate_size > 0:
            # Now, we check to see if we can afford to simulate the proper STG
            # to actually find the attractors. If we can't, we'll simulate the
            # reduction, which will give bounds on the number of motif-free
            # complex attractors.
            if self.rspace_update_primes is not None:
                simulate_size = len(self.rspace_update_primes)
            else:
                simulate_size = len(self.reduced_primes)

            if simulate_size < max_simulate_size:
                self.find_no_motif_attractors()
                if len(self.no_motif_attractors) == 0:
                    self.terminal = "no"
                else:
                    self.terminal = "yes"
                    self.conserved_functions = sm_rspace.attractor_space_candidates(self.stable_motifs,
                                                                         self.time_reverse_stable_motifs)
            else:
                #print("STG is too large to simulate ("+
                #str(simulate_size)+"/"+str(max_simulate_size)+
                #"). We will attempt reduction methods. Increase max_simulate_size to force simulation.")
                self.delprimes, self.attractor_constants = mediator_reduction(self.reduced_primes)
                self.delprimes, nc = deletion_reduction(self.delprimes, max_in_degree=max_in_degree)
                self.attractor_constants.update(nc)

                # Before building any STGs, let's see if we've already identified
                # that a stable motif must stabilize based on the reduction.
                for sm in self.stable_motifs:
                    sat = True
                    for k,v in sm.items():
                        if (k,v) in self.attractor_constants.items():
                            continue
                        elif k not in self.delprimes:
                            continue
                        sat = False
                        break
                    if sat:
                        self.terminal = "no"
                        #print("The reduction indicates that the branch is not terminal. No need to simulate.")
                        self.attractor_dict_list = self.generate_attr_dict_list()
                        return
                if len(self.delprimes) < max_simulate_size_vc:
                    #print("Simulating deletion reduction ("+str(len(self.delprimes))+" nodes)...")
                    self.find_deletion_no_motif_attractors(max_stable_motifs=max_stable_motifs)
                    if len(self.deletion_no_motif_attractors) == 0 and self.terminal != "yes":
                        self.terminal = "no"
                    # else:
                    #     self.terminal = "possible"
                else:
                    print("The STG is still too large ("+str(len(self.delprimes))+").")
                    print("Further analysis of this branch is needed.")
        self.attractor_dict_list = self.generate_attr_dict_list()

    def merge_source_motifs(self):
        """Merges stable motifs (and time-reversal stable motifs) that correspond to source nodes, e.g. A*=A, into combined motifs to
        avoid combinatorial explosion. For example, A*=A, B*=B, C*=C produces six motifs that can stabilize in 8 ways; without
        merging, these 8 combinations lead to 8*3!=48 successions because they can be considered in any order. This is silly because
        source nodes all stabilize simultaneously.

        We will assume that stable motifs and time reverse stable motifs have already been computed.

        Note that a source node in the forward time system is a source node in the time reverse system as well.
        This follows from A* = A => A- = ~(A*(A=~A)) = ~(~A) = A.

        If A* = A or X (i.e., A=1 is a stable motif), then A- = ~(~A | X) = A & ~X, so A=0 is a time-reverse stable motif. A similar
        argument applies for the A=0 stable motif. Thus, a motif is only a source motif if it is also a time-reverse motif.
        """
        source_motifs = [x for x in self.stable_motifs if len(x) == 1 and x in self.time_reverse_stable_motifs]
        if source_motifs == []:
            return
        self.source_independent_motifs = [x for x in self.stable_motifs if not x in source_motifs]

        source_vars = list(set([next(iter(x.keys())) for x in source_motifs])) # a list of source nodes

        self.merged_source_motifs = []
        for state in it.product([0,1],repeat=len(source_vars)):
            self.merged_source_motifs.append({v:x for v,x in zip(source_vars,state)})

    def simple_merge_source_motifs(self,primes,MPBN_update=False):
        """Merges stable motifs (and time-reversal stable motifs) that correspond to source nodes, e.g. A*=A, into combined motifs to
        avoid combinatorial explosion. For example, A*=A, B*=B, C*=C produces six motifs that can stabilize in 8 ways; without
        merging, these 8 combinations lead to 8*3!=48 successions because they can be considered in any order. This is silly because
        source nodes all stabilize simultaneously.

        Assumes that stable_motifs have already been computed,
        but time_reverse_primes and time_reverse_stable_motifs are not.

        To be used in the case of MPBN update.

        Parameters
        ----------
        primes : pyboolnet primes dictionary
            pyboolnet update rules whose source node stable motifs are to be merged.

        MPBN_update : bool
            Whether MPBN update is used instead of general asynchronous update
            (see Pauleve et al. 2020)(the default is False).

        Returns
        -------
        self.source_independent_motifs : list of dictionaries
            list of stable motifs that are not source motifs
            [{'node1':bool,'node2':bool, ...}, {'node3':bool,'node4':bool, ...}, ...]
        self.merged_source_motifs : list of dictionaries
            list of group of source motifs fixed at the same time
            [{'source_node1':bool,'source_node2':bool, ...}, ...]

        """
        assert MPBN_update == True, "This function is for MPBN update only"

        # source nodes will have update rule of the form 'A':[[{'A':0}],[{'A':1}]]
        source_vars = []
        for x in primes.keys():
            if primes[x] == [[{x:0}],[{x:1}]]:
                source_vars.append(x)

        source_motifs = []
        for x in source_vars:
            source_motifs.append({x:0})
            source_motifs.append({x:1})

        if source_motifs == []:
            return

        self.source_independent_motifs = [x for x in self.stable_motifs if not x in source_motifs]

        self.merged_source_motifs = []
        for state in it.product([0,1],repeat=len(source_vars)):
            self.merged_source_motifs.append({v:x for v,x in zip(source_vars,state)})

    # def test_rspace(self, search_partial_STGs=True):
    #     STG=pyboolnet.state_transition_graphs.primes2stg(self.rspace_update_primes,"asynchronous")
    #     steady_states,complex_attractors=pyboolnet.attractors.compute_attractors_tarjan(STG)
    #     names = sorted(self.rspace_update_primes)
    #     attractors = complex_attractors+[[s] for s in steady_states]
    #     self.build_rspace_attractor_candidates(attractors)

    # def build_rspace_attractor_candidates(self,attractors):
    #     self.rspace_attractor_candidates = []
    #     for attractor in attractors:
    #         possible_rspace_attractor = True
    #         for state in attractor:
    #             # state_dict = {** sm_format.statestring2dict(state,names),**self.fixed_rspace_nodes}
    #             # if pyboolnet.boolean_logic.are_mutually_exclusive(self.rspace_constraint,
    #             #                                                  sm_format.implicant2bnet(state_dict)):
    #             if sm_rspace.partial_state_contradicts_rspace(sm_format.statestring2dict(state,names),self.rspace):
    #                 possible_rspace_attractor = False
    #                 break
    #         if possible_rspace_attractor:
    #             self.rspace_attractor_candidates.append(attractor)
    #
    #     if len(self.rspace_attractor_candidates) == 0:
    #         self.terminal = "no"

    def build_K0(self):
        """Helper function for smart STG building. Builds initial set of nodes
        that are not part of any motif-avoidant attractor.

        Returns
        -------
        set of str
            Statestrings that do not belong to any motif-avoidant attractor.

        """
        K = set()
        for sm in self.stable_motifs:
            fill_vars = [k for k in self.reduced_primes if not k in sm]
            for fills in it.product(['0','1'],repeat = len(fill_vars)):
                s = ''
                fi = 0
                for k in self.reduced_primes:
                    if k in sm:
                        s += str(sm[k])
                    else:
                        s += fills[fi]
                        fi += 1
                K.add(s)
        return K

    def in_motif(self,ss,names):
        """Tests whether the (partial) state ss is in any stable motifs

        Parameters
        ----------
        ss : str
            Statestring (possibly on a subspace).
        names : list of str
            Variable names ordered to correspond to the positions of ss.

        Returns
        -------
        bool
            Whether ss is in any stable motif of the reduced system.

        """
        for sm in self.stable_motifs:
            smin = True
            for i,r in enumerate(names):
                if r in sm and not int(ss[i]) == sm[r]:
                    smin = False
            if smin: return True
        return False

    def build_inspace(self,ss,names, tr_stable_motifs = None):
        """Helper function for smart STG building. List all time reversal stable
        motifs to which (partial) state ss belongs.

        Parameters
        ----------
        ss : str
            Statestring (possibly on a subspace).
        names : list of str
            Variable names ordered to correspond to the positions of ss.
        tr_stable_motifs : list of partial state dictionaries
            Time reverse stable motifs. If None, use all time reverse stable
            motifs in the reduced system (the default is None).

        Returns
        -------
        list of partial state dictionaries
            Time reverse stable motifs that are active in the state ss.

        """
        inspaces = []
        if tr_stable_motifs is None:
            tr_stable_motifs = self.time_reverse_stable_motifs

        for ts in tr_stable_motifs:
            tsin = True
            for i,r in enumerate(names):
                if r in ts and not int(ss[i]) == ts[r]:
                    tsin = False
            if tsin: inspaces.append(ts)
        return inspaces

    def build_deletion_STG(self,max_stable_motifs=10000):
        """Build a piece of the STG that is guaranteed to contain all
        motif-avoidant attractors of the deletion projection. Complex attractors
        found here may be spurious.

        Parameters
        ----------
        max_stable_motifs : int
            Maximum number of output lines for pyboolnet to process from the
            AspSolver (the default is 10000).

        """
        names = sorted(self.delprimes)
        name_ind = {n:i for i,n in enumerate(names)}
        trprimes = sm_time.time_reverse_primes(self.delprimes)
        trsms = pyboolnet.trap_spaces.compute_trap_spaces(trprimes,"max",max_output=max_stable_motifs)

        if self.rspace_update_primes is not None:
            delrnames = [x for x in sorted(self.rspace_update_primes) if x in self.delprimes]
            rname_ind = {n:i for i,n in enumerate(names) if n in delrnames}
            fixed = {k:v for k,v in self.fixed_rspace_nodes.items()}
        else:
            rnames = names.copy()
            rname_ind = name_ind.copy()
            fixed = {}
        sim_names = [x for x in names if not x in fixed]

        #K = self.build_K0()
        K = set()
        self.deletion_STG = nx.DiGraph()

        inspace_dict = {}
        t = 0
        T = 1
        # note: product order gives s counting up in binary from 00..0 to 11..1
        for s in it.product(['0','1'],repeat=len(sim_names)):
            sl = ['']*len(names)
            j = 0
            for i in range(len(names)):
                if names[i] in fixed:
                    sl[i] = str(fixed[names[i]])
                else:
                    sl[i] = s[j]
                    j += 1

            ss = ''.join(sl)

            if ss in K: continue
            if self.in_motif(ss,names): continue

            simstate = True

            inspace = self.build_inspace(ss,names,tr_stable_motifs = trsms)
            inspace_dict[ss] = inspace

            self.deletion_STG.add_node(ss) # might end up removing the node later
            for i,r in enumerate(names):
                nri = int(not int(ss[i]))
                # if any p below is satisfied, we get a change of state
                # the value of the new r will be equal to nri
                for p in self.delprimes[r][nri]:
                    psat = True
                    for k,v in p.items():
                        if not int(ss[name_ind[k]]) == v:
                            psat = False
                            break
                    if psat: # state change verified
                        child_state_list = list(ss)
                        child_state_list[i] = str(nri)
                        child_state = ''.join(child_state_list)

                        # Check if changed something that should be fixed or landed in K
                        # If not, check if we left a TR stable motif
                        prune = r in fixed or child_state in K
                        # next we check to see if we've left the rspace
                        # note that we don't have to check rspace[0], as this
                        # is handled by checking r in fixed
                        if not prune:
                            stdict = sm_format.statestring2dict(child_state,names)
                            prune = not sm_rspace.partial_state_contradicts_rspace(stdict, self.rspace[1:])

                        # next, we check to see if we left a TR motif
                        if not prune:
                            if not child_state in inspace_dict:
                                inspace_dict[child_state] = self.build_inspace(child_state,names,tr_stable_motifs = trsms)
                            prune = not inspace_dict[child_state] == inspace


                        # By here, prune is TRUE if we left a TR motif or are in K
                        if prune:
                            # prune the STG and stop simulating ss
                            simstate = False
                            rnodes = list(nx.bfs_tree(self.deletion_STG,ss,reverse=True).nodes())
                            K.update(rnodes)
                            self.deletion_STG.remove_nodes_from(rnodes)
                        else:
                            self.deletion_STG.add_edge(ss,child_state)
                        break # we know the ss at r changed, no need to check more primes
                if not simstate: break # don't check other vars: already found ss -> K

    def find_deletion_no_motif_attractors(self,max_stable_motifs=10000):
        """Identify motif-avoidant attractors in the deletion projection.

        Parameters
        ----------
        max_stable_motifs : int
            Maximum number of output lines for pyboolnet to process from the
            AspSolver (the default is 10000).

        """
        if self.deletion_STG is None:
            self.build_deletion_STG(max_stable_motifs=max_stable_motifs)

        # Note: fixed points of the deletion system are fixed points of the
        # undeleted system, so we ignore these as they must contain stable motifs
        if len(list(self.deletion_STG.nodes())) > 0:
            candidates = [x for x in nx.attracting_components(self.deletion_STG) if len(x) > 1]
        else:
            candidates = []
        self.deletion_no_motif_attractors = []

        # next, we see if any of these activate stable motifs
        names = sorted(self.delprimes)
        for att in candidates:
            no_motif = True
            for s in att:
                # The following check stems from the result that a stable motif
                # is active in an attractor of the original system iff its projection
                # is active in the projected attractor in the deletion-reduced system
                st = sm_format.statestring2dict(s,names)
                st.update(self.attractor_constants)
                if any(not sm_doi.fixed_excludes_implicant(st,sm) for sm in self.stable_motifs):
                    no_motif = False
                    break

            if no_motif:
                self.deletion_no_motif_attractors.append(att)

    def build_partial_STG(self):
        """Build a piece of the STG that is guaranteed to contain all
        motif-avoidant attractors of the reduction.

        """
        names = sorted(self.reduced_primes)
        name_ind = {n:i for i,n in enumerate(names)}

        if self.rspace_update_primes is not None:
            rnames = sorted(self.rspace_update_primes)
            rname_ind = {n:i for i,n in enumerate(names) if n in rnames}
            fixed = self.fixed_rspace_nodes
        else:
            rnames = names.copy()
            rname_ind = name_ind.copy()
            fixed = {}

        sim_names = [x for x in names if not x in fixed]


        #K = self.build_K0()
        K = set()
        self.partial_STG = nx.DiGraph()

        inspace_dict = {}
        t = 0
        T = 1
        # note: product order gives s counting up in binary from 00..0 to 11..1
        for s in it.product(['0','1'],repeat=len(sim_names)):
            sl = ['']*len(names)
            j = 0
            for i in range(len(names)):
                if names[i] in fixed:
                    sl[i] = str(fixed[names[i]])
                else:
                    sl[i] = s[j]
                    j += 1

            ss = ''.join(sl)

            if ss in K: continue
            if self.in_motif(ss,names): continue

            simstate = True

            inspace = self.build_inspace(ss,names)
            inspace_dict[ss] = inspace

            self.partial_STG.add_node(ss) # might end up removing ss later
            for i,r in enumerate(names):
                nri = int(not int(ss[i]))
                # if any p below is satisfied, we get a change of state
                # the value of the new r will be equal to nri
                for p in self.reduced_primes[r][nri]:
                    psat = True
                    for k,v in p.items():
                        if not int(ss[name_ind[k]]) == v:
                            psat = False
                            break
                    if psat: # state change verified
                        child_state_list = list(ss)
                        child_state_list[i] = str(nri)
                        child_state = ''.join(child_state_list)

                        # Check if changed something that should be fixed or landed in K
                        # If not, check if we left a TR stable motif
                        prune = r in fixed or child_state in K

                        # next we check to see if we've left the rspace
                        # note that we don't have to check rspace[0], as this
                        # is handled by checking r in fixed
                        if not prune:
                            stdict = sm_format.statestring2dict(child_state,names)
                            prune = not sm_rspace.state_in_rspace(stdict, self.rspace[1:])

                        # next, we check to see if we left a TR motif
                        if not prune:
                            if not child_state in inspace_dict:
                                inspace_dict[child_state] = self.build_inspace(child_state,names)
                            prune = not inspace_dict[child_state] == inspace


                        # By here, prune is TRUE if we left a TR motif or are in K
                        if prune:
                            # prune the STG and stop simulating ss
                            simstate = False
                            rnodes = list(nx.bfs_tree(self.partial_STG,ss,reverse=True).nodes())
                            K.update(rnodes)
                            self.partial_STG.remove_nodes_from(rnodes)
                        else:
                            self.partial_STG.add_edge(ss,child_state)
                        break # we know the ss at r changed, no need to check more primes
                if not simstate: break # don't check other vars: already found ss -> K

    def find_no_motif_attractors(self):
        """Find attractors of the reduction that are not present in any of its
        subreductions.

        """
        if self.partial_STG is None:
            self.build_partial_STG()
        if len(list(self.partial_STG.nodes())) > 0:
            self.no_motif_attractors = list(nx.attracting_components(self.partial_STG))
        else:
            self.no_motif_attractors = []

    def generate_attr_dict_list(self, MPBN_update=False):
        """Generate a list of attractors that are present in the reduction, but
        not in any of its subreductions.

        Parameters
        ----------
        MPBN_update : bool
            Whether MPBN update is used instead of general asynchronous update
            (see Pauleve et al. 2020)(the default is False).

        Returns
        -------
        list of dictionaries
            Dictionaries corresponding to attractors that are in this reductions, but
            not in any of its subreductions (if it has any). Each describes the node
            states in the attractors according to the following
            1 variable is "ON"
            0 variable is "OFF"
            X variable is known to oscillate
            ? at least one such variable must oscillate
            ! the attractor may be false; if it is genuine, at least
              one such variable must oscillate
        """
        attractors_dict={}

        #the reduction is not terminal --> no attractor
        if self.terminal == 'no':
            return []#'not terminal reduction' #I should replace this with an empty dict

        nodes_sorted = sorted(list(set(self.logically_fixed_nodes.keys()) | set(self.reduced_primes.keys()))) #steady state

        node_state_dict = self.logically_fixed_nodes.copy()
        if self.fixed_rspace_nodes is not None:
            node_state_dict.update(self.fixed_rspace_nodes)

        # Found a steady state (will always be terminal)
        if len(node_state_dict) == len(nodes_sorted):
            assert self.terminal == 'yes', "Found non-terminal steady state. This should not be possible!"
            node_state_dict = {k:v for k,v in sorted(node_state_dict.items())}
            return [node_state_dict]
        #non_fixed_nodes = sorted(list(set(nodes_sorted)-set(node_state_dict.keys())))
        non_fixed_nodes = [x for x in nodes_sorted if x not in node_state_dict]

        #the reduction is only possibly terminal
        if self.terminal=='possible':
            for n in non_fixed_nodes: #non-stabilized nodes
                node_state_dict[n] = '!'
                node_state_dict = {k:v for k,v in sorted(node_state_dict.items())}
            return [node_state_dict]

        #the reduction is definitely terminal
        elif self.terminal=='yes':
            #the reduction is terminal, not all nodes are fixed
            #and it is MPBN update.
            if MPBN_update == True:
                for n in non_fixed_nodes:
                    node_state_dict[n]='X'
                node_state_dict = {k:v for k,v in sorted(node_state_dict.items())}
                return [node_state_dict]

            #the reduction is terminal, not all nodes are fixed
            #and there is NO complex attractor mapped out
            elif self.no_motif_attractors is None:
                for n in non_fixed_nodes:
                    node_state_dict[n]='?'
                node_state_dict = {k:v for k,v in sorted(node_state_dict.items())}
                return [node_state_dict]

            #the reduction is terminal, not all nodes are fixed
            #there are complex attractors mapped out:
            attr_list=[]
            for complex_attractor in self.no_motif_attractors:
                ca_dict = sm_format.statelist2dict(sorted(self.reduced_primes),complex_attractor)
                ns = node_state_dict.copy()
                ns.update(ca_dict)
                ns = {k:v for k,v in sorted(ns.items())}
                attr_list.append(ns)
                # # we check if there are stabilized nodes within the complex attractor
                # ca=self.find_constants_in_complex_attractor(complex_attractor)
                # node_state_dict=self.logically_fixed_nodes.copy()
                # for i in range(len(non_fixed_nodes)): #non-stabilized nodes
                #     node_state_dict[non_fixed_nodes[i]]=ca[i]
                #
                # attr_list.append(node_state_dict)
            return attr_list

    def find_constants_in_complex_attractor(self,c):
        """Given a set of strings representing the states of a complex attractor the function finds the nodes
        that are constant in the full complex attractor.

        Parameters
        ----------
        c : a set of binary strings
            Set of statestrings, e.g. set(['000', '010', '100']).

        Returns
        -------
        list of str
            An array consisting of 0s, 1s, and Xs. X represents an oscillating
            node, and the 0s and 1s represent nodes stabilized to those states.

        """
        import numpy as np
        ca=np.array([np.fromiter(i, int, count=len(i)) for i in c])
        attr=np.array(['X' for i in range(len(ca[0]))])
        sum_a0=ca.sum(axis=0)
        attr[np.where(sum_a0==0)[0]]=0
        attr[np.where(sum_a0==len(ca))[0]]=1
        return attr

    def summary(self,show_original_rules=True,hide_rules=False,show_explicit_permutations=False):
        """Print a summary of the reduction.

        Parameters
        ----------
        show_original_rules : bool
            Show rules of the unreduced system (the default is True)?
        hide_rules : bool
            Hide rules of the reduced system (the default is False)?
        show_explicit_permutations : bool
            Show motif permutations explicitly, instead of by index (the default is False)?

        """
        print("Motif History:",self.motif_history)
        print()
        print("Logically Fixed Nodes:",self.logically_fixed_nodes)
        print()
        if hide_rules:
            pass
        elif not self.motif_history == []:
            print("Reduced Update Rules:")
            sm_format.pretty_print_prime_rules(self.reduced_primes)
        else:
            if show_original_rules:
                print("Original Update Rules:")
                sm_format.pretty_print_prime_rules(self.reduced_primes)
            else:
                print("The update rules are not reduced.")
        print()
        if self.terminal == "no":
            if self.merged_source_motifs is None:
                print("At least one additional stable motif must stabilize.")
                print()
                print("Stable motifs:", self.stable_motifs)
            else:
                print("Source node values are not yet specified for the following nodes:",
                       ', '.join(sorted([k for k in self.merged_source_motifs[0]])))
                print()
                if self.source_independent_motifs == []:
                    print("There are no source-independent stable motifs.")
                else:
                    print("The following stable motifs exist independently of the source configuration:")
                    print(self.source_independent_motifs)
        elif self.terminal == "yes":
            if len(self.reduced_primes) > 0:
                print("There is a complex attractor in this reduced system in which no additional stable motifs activate.")
                print("At least some of the following must oscillate in such an attractor:")
                print(list(self.reduced_primes.keys()))
            else:
                print("This branch terminates in a steady state.")
        elif self.terminal == "possible":
            print("Some or none of these stable motifs may stabilize:",
                  self.stable_motifs)
            print()
            if not self.fixed_rspace_nodes is None:
                print("If no more stable motifs stabilize, these node states must be fixed:",
                      self.fixed_rspace_nodes)
                print()
                print("In addition, the following must stabilize to TRUE:")
                print(self.reduced_rspace_constraint)
                print()
                print("In this case, the unfixed nodes update according to the following rules:")
                sm_format.pretty_print_prime_rules(self.rspace_update_primes)

        if not self.conserved_functions is None:
            print()
            if len(self.conserved_functions) > 0:
                print("Found the following functions that are constant on attractors in this branch:")
                for x in self.conserved_functions:
                    if len(x) > 0:
                        sm_format.pretty_print_rspace([x],silent=False)
                        print()
            else:
                print("Unable to find non-trivial conserved functions for attractors in this branch.")
                print()
            if not self.no_motif_attractors is None:
                if len(self.no_motif_attractors) > 0:
                    print("Found the following complex attractors that do not lock in additional stable motifs:")
                    for x in self.no_motif_attractors:
                        print(x)

        if len(self.merged_history_permutations) > 0:
            print()
            print("This branch contains the following motif_history permutation(s):")
            if show_explicit_permutations:
                for x in self.merged_history_permutations: print([self.motif_history[i] for i in x])
            else:
                for x in self.merged_history_permutations: print(x)
