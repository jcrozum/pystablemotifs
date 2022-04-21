import networkx as nx
import itertools as it
import matplotlib
import matplotlib.pyplot as plt

from copy import deepcopy

import pystablemotifs.reduction as sm_reduction
import pystablemotifs.format as sm_format
import pystablemotifs.drivers as sm_doi

class SuccessionDiagram:
    """Class describing the succession diagram of a Boolean system. See, e.g.,
    Zanudo and Albert (2015) or Rozum et al. (2021).

    Attributes
    ----------
    motif_reduction_dict : dictionary
        MotifReduction-valued dictionary with integer (index) keys (see reduction.py).
    digraph : networkx.DiGraph
        Topological structure of hte succession diagram. Nodes are integers that
        align with the enteries of motif_reduction_dict.

    """

    def __init__(self):
        self.motif_reduction_dict = {}
        self.digraph = nx.DiGraph()

    def get_motifs(self):
        """Extract the stable motifs of a system and its reduced networks from its
        attractor repertoire. Notably, these include both the system's primary stable
        motifs and conditionally stable motifs (see, e.g., Deritei et al. 2019).

        Returns
        -------
        list of dictionaries
            Stable motifs that appear in the system or during reduction (in no
            particular order).

        """
        SM_set = set([])
        for reduction in self.motif_reduction_dict.values():
            if len(reduction.motif_history)>0:
                SM_set.add(frozenset(reduction.motif_history[-1].items()))
        return [dict(s) for s in SM_set]

    def find_motif_permutation(self,motif_history):
        """Check whether some permutation of the input motif_history is already
        represented in the succession diagram. If so, return the preexisting
        reduction's index and the permutation that maps between the two histories.

        Parameters
        ----------
        motif_history : list of partial state dictionaries
            Stable motifs that can lock in to give a given reduced network (in
            order).

        Returns
        -------
        reduction_index : int
            Index of the preexisting reduced network. This value is None if no such
            reduced network exists.
        permutation : list of int
            Permutation that maps the preexisting history to the input history.
            This value is None if no such history exists.

        """
        for reduction_index,mr in self.motif_reduction_dict.items():
            if len(mr.motif_history) == len(motif_history):
                if all([x in mr.motif_history for x in motif_history]):
                    permutation = []
                    for x in motif_history:
                        permutation.append(mr.motif_history.index(x))
                    return reduction_index,permutation # We already have an equivalent motif in the list
        return None,None

    def add_motif_permutation(self,reduction_index,permutation):
        """Adds a permutation of a preexisting stable motif history to a precomputed
        MotifReduction object.

        Parameters
        ----------
        reduction_index : int
            Index of the preexisting reduced network.
        permutation : list of int
            Permutation that maps the preexisting history to the input history.

        """
        self.motif_reduction_dict[reduction_index].merged_history_permutations.append(permutation)
        for child in nx.topological_sort(self.digraph):
            for parent in self.digraph.predecessors(child):
                for parent_perm,child_perm in it.product(
                    self.motif_reduction_dict[parent].merged_history_permutations,
                    self.motif_reduction_dict[child].merged_history_permutations):
                    new_perm = child_perm.copy()
                    for i,p in enumerate(parent_perm):
                        new_perm[i] = child_perm[p]
                    if not new_perm in self.motif_reduction_dict[child].merged_history_permutations:
                        self.motif_reduction_dict[child].merged_history_permutations.append(new_perm)

    def find_equivalent_reduction(self,fixed):
        """Extracts the MotifReduction object that has the frozen node values
        specified by fixed, if such an object exists (returns None otherwise).

        Parameters
        ----------
        fixed : partial state dictionary
            Nodes values that have been fixed and reduced by stable motifs and their
            logical domain of influence.

        Returns
        -------
        MotifReduction
            Reduced network that has the frozen node values specified by fixed,
            if such an object exists (returns None otherwise).

        """
        for reduction in self.motif_reduction_dict.values():
            if reduction.logically_fixed_nodes == fixed:
                return reduction
        return None

    def add_motif_reduction(self,motif_reduction):
        """Inserts a given MotifReduction into the succession diagram. Does not
        check for consistency, but will insert a properly constructed MotifReduction
        into the correct place in the diagram.

        Parameters
        ----------
        motif_reduction : MotifReduction
            Reduced network to be appended to the succession diagram.

        """
        if self.motif_reduction_dict == {}:
            self.unreduced_primes = motif_reduction.reduced_primes

        # note: N is computed BEFORE the new reduction is added,
        # so it will be the reduction index AFTER the reduction is added.
        N = len(self.motif_reduction_dict)
        new_set = set([frozenset(tuple(x.items())) for x in motif_reduction.motif_history])
        if N == 0:
            self.digraph.add_node(0)
        else:
            for i,reduction in self.motif_reduction_dict.items():
                old_set = set([frozenset(tuple(x.items())) for x in reduction.motif_history])

                # see if we're adding a parent of an existing reduction
                if len(old_set) == len(new_set) + 1:
                    diff_set = old_set - new_set
                    if len(diff_set) == 1:
                        missing_motif = dict(diff_set.pop())
                        if missing_motif in motif_reduction.stable_motifs:
                            self.digraph.add_edge(N,i)
                # see if we're adding a child of an existing reduction
                elif len(old_set) == len(new_set) - 1:
                    diff_set = new_set - old_set
                    if len(diff_set) == 1:
                        missing_motif = dict(diff_set.pop())
                        if missing_motif in reduction.stable_motifs or len(new_set) == 1: # valid motif OR source combination
                            self.digraph.add_edge(i,N)

        self.motif_reduction_dict[N] = motif_reduction
        self.add_motif_permutation(N,list(range(len(motif_reduction.motif_history))))

    def reductions_indices_with_states(self,logically_fixed,optimize=True):
        """Find all reductions (by index) that have the nodes states specified
        logically fixed.

        Parameters
        ----------
        logically_fixed : partial state dictionary
            Nodes states that should be fixed in all returned network reductions.
        optimize : bool
            Whether to remove reduced networks that are subnetworks of valid
            reductions. This is generally recommended so as to obtain the most
            parsimonious control strategies (the default is True).

        Returns
        -------
        list of int
            Indices of reduced networks that have the appropriate fixed states.

        """
        if not optimize:
            # NOTE: This finds all reductions, not just those closest to the root
            target_indices = []
            for i,reduction in self.motif_reduction_dict.items():
                if logically_fixed.items() <= reduction.logically_fixed_nodes.items():
                    target_indices.append(i)
            return target_indices
        else:
            # The unoptimized result
            target_indices_unoptimized = self.reductions_indices_with_states(logically_fixed,optimize=False)
            target_indices_all = target_indices_unoptimized.copy()

            # Add nodes that inevitably reach nodes in the unoptimized result
            nodes_to_consider = [x for x in self.digraph if (not x in target_indices_all
                and self.motif_reduction_dict[x].terminal == "no")]
            for i in nodes_to_consider:
                bad_sinks = [x for x in nx.descendants(self.digraph,i) | set([i]) if (
                    not self.motif_reduction_dict[x].terminal == "no" # is a sink and . . .
                    and not x in target_indices_unoptimized)] # does not have the right nodes fixed

                if len(bad_sinks) == 0:
                    target_indices_all.append(i)

            # remove nodes who have ancestors that would work instead, as those
            # will be closer to the root and easier to control
            target_indices = []
            for target_index in target_indices_all:
                if set(nx.ancestors(self.digraph,target_index)) & set(target_indices_all) == set():
                    target_indices.append(target_index)

            return target_indices

    def reduction_drivers(self,target_index,method='internal',max_drivers=None,GRASP_iterations=None):
        """Find control strategies that lead to the reduced network specified by
        the target index. Several control strategies are implemented. See
        succession.SuccessionDiagram.reprogram_to_trap_spaces for a detailed
        description of control methods available. Generally, this method should
        not be used directly. Instead, use reprogram_to_trap_spaces.

        Parameters
        ----------
        target_index : int
            Index of the target reduced network.
        method : str
            One of 'internal', 'minimal', or 'GRASP'. See
            succession.SuccessionDiagram.reprogram_to_trap_spaces for details.
        max_drivers : int
            Maximum number of driver nodes to consider (not used in GRASP methods).
            If none, the upper limit is given by the number of free variables
            (the default is None).
        GRASP_iterations : int
            Number of times to construct GRASP driver sets; only used in GRASP
            methods. If none, the number of iterations is chosen based on the
            network size (the default is None).

        Returns
        -------
        list
            Control strategies found; interpretation depends on method selected
            See succession.SuccessionDiagram.reprogram_to_trap_spaces for details.

        """
        methods = ['internal','minimal','GRASP']
        assert method in methods, ' '.join(["method argument of reduction_drivers must be among",str(methods)])
        drivers = []
        for path in nx.all_simple_paths(self.digraph,0,target_index):
            path_motif_history=[]
            path_drivers = []
            for ind in path:
                if ind == 0:
                    ind_prev = ind
                    continue
                path_motif_history += [x for x in self.motif_reduction_dict[ind].motif_history if not x in path_motif_history]

                if method == 'internal':
                    history_drivers = sm_doi.internal_drivers(path_motif_history[-1],
                        self.motif_reduction_dict[ind_prev].reduced_primes,
                        max_drivers=max_drivers)
                elif method == 'GRASP':
                    history_drivers = sm_doi.GRASP(path_motif_history[-1],
                        self.motif_reduction_dict[ind_prev].reduced_primes,
                        GRASP_iterations = GRASP_iterations)
                    if len(history_drivers) == 0:
                        history_drivers = [path_motif_history[-1]]
                elif method == 'minimal':
                    history_drivers = sm_doi.minimal_drivers(path_motif_history[-1],
                        self.motif_reduction_dict[ind_prev].reduced_primes,
                        max_drivers=max_drivers)

                path_drivers.append(history_drivers)

                ind_prev = ind
            if method == 'internal':
                # merge control sets along the path
                for control_sequence in it.product(*path_drivers):
                    control_set = {k:v for x in control_sequence for k,v in x.items()}
                    if not control_set in drivers:
                        drivers.append(control_set)
            elif method == 'minimal' or method == 'GRASP':
                drivers.append(path_drivers)
        return drivers

    def reprogram_to_trap_spaces(self,logically_fixed,target_method='history',driver_method='internal',max_drivers=None,GRASP_iterations=None,GRASP_score_override=None):
        """Find driver sets that lead to fixing the node states specified.

        Parameters
        ----------
        logically_fixed : partial state dictionary
            Targeted fixed nodes.
        target_method : str
            Either 'history' or 'merge'; see Notes below for details.
        driver_method : str
            Either 'internal', 'minimal', or 'GRASP' see Notes below for details.
        max_drivers : int
            Maximum number of driver nodes to consider (not used in GRASP methods).
            If none, the upper limit is given by the number of free variables
            (the default is None).
        GRASP_iterations : int
            Number of times to construct GRASP driver sets; only used in GRASP
            methods. If none, the number of iterations is chosen based on the
            network size (the default is None).
        GRASP_score_override : function
            Optional heuristic score function override (see drivers.GRASP
            for details). Only used in GRASP methods (the default is None).

        Returns
        -------
        list
            Control strategies found; interpretation depends on method selected
            See Notes below for details.

        Notes
        -----
        The various combinations of target_method and driver_method options result
        in different control strategies, which are outlined below.

        target_method = history, driver_method = internal:
        Finds all shortest stable motif histories that result in the target node states
        being logically fixed. Each stable motif is searched for internal driver nodes.
        The resulting internal drivers are combined into a single  control set. The
        return value consists of all such control sets for all  stable motif histories.
        Each control set eventually becomes self-sustaining.

        target_method = history, driver_method = minimal:
        Similar to the history method, except the search for stable motif drivers
        includes external driver nodes for the motif and does not extend to driver sets
        of larger size once one driver set has been found for a motif. Because the
        search includes external driver nodes, special care must be taken in interpreting
        the effect of the drivers, as their influence may impact the effect of motifs
        stabilizing. Thus, the control is only guaranteed to work if the interventions
        are temporary and implemented in the order specified by the motif history.

        For this reason, the output consists of lists of ordered interventions.
        Each element of the return value is a list of lists of dictionaries. Each
        element of the return value represents a control strategy. To implement such
        a strategy, select a dictionary from the first element of the strategy and
        fix the node states it specifies until their influence has propagated through
        the system. Then repeat this process iteratively for each element of the strategy
        list, in order. For example, if
        nonredundant_drivers = [ [[{'xD':1,'xE=1'}]], [[{'xA':1},{'xB':1}],[{'xC':1}]] ]
        then there are two control strategies available:
        1) fix xD=xE=1 temporarily and
        2) first fix either xA=1 or xB=1 temporarily, then fix xC=1 temporarily.

        target_method = history, driver_method = GRASP:
        The same as history, minimal, except external driver nodes are searched for
        using the GRASP algorithm using GRASP_iterations iterations.

        target_method = merge, driver_method = internal:
        Finds all shortest stable motif histories that result in the target node states
        being logically fixed. All node states in the motifs in the history are merged
        into a stable module dictionary. This is then searched for internal driver
        nodes. Each element of the return value is a dictionary corresponding to a
        control set. Each control set eventually becomes self-sustaining.

        target_method = merge, driver_method = minimal:
        Similar to the merge method, except the search for drivers is conducted over
        all nodes, not just those internal to the merged stable module. Furthermore,
        the search is truncated when a control set is found such that the search does
        not proceed to driver sets larger than the smallest found. Each element of
        the return value is a dictionary corresponding to a control set. The control
        sets are only guaranteed to result in activation of the target if they are
        temporary interventions.

        target_method = merge, driver_method = GRASP:
        The same as merge, minimal, except external driver nodes are searched for
        using the GRASP algorithm using GRASP_iterations iterations.

        """

        #methods = ['history','merge','minimal_history','minimal_merge']
        target_methods = ['history','merge']
        driver_methods= ['internal','minimal','GRASP']
        assert target_method in target_methods, ' '.join(["target_method argument of reprogram_to_trap_spaces must be among",str(target_methods)])
        assert driver_method in driver_methods, ' '.join(["driver_method argument of reprogram_to_trap_spaces must be among",str(driver_methods)])
        drivers = []

        if driver_method == 'GRASP' and GRASP_iterations is None:
            if target_method == 'merge':
                GRASP_iterations = len(self.unreduced_primes)**2
            if target_method == 'history':
                GRASP_iterations = 2*len(self.unreduced_primes)

        target_indices = self.reductions_indices_with_states(logically_fixed)

        if target_method == 'history':
            for target_index in target_indices:
                drivers += self.reduction_drivers(target_index,max_drivers=max_drivers,GRASP_iterations=GRASP_iterations,method=driver_method)
        elif target_method == 'merge':
            for target_index in target_indices:
                target_history = self.motif_reduction_dict[target_index].motif_history
                motif_merger = {k:v for d in target_history for k,v in d.items()}

                # Because we are potentially dealing with external drivers, we
                # want to make sure the external drivers do not interfere with
                # the motif's ability to drive downstream targets
                # if driver_method != 'internal':
                #     motif_merger.update(logically_fixed)

                if driver_method == 'GRASP':
                    merger_drivers = sm_doi.GRASP(motif_merger,self.unreduced_primes,GRASP_iterations)
                    if len(merger_drivers) == 0:
                        merger_drivers = [motif_merger.copy()]
                elif driver_method == 'minimal':
                    merger_drivers = sm_doi.minimal_drivers(motif_merger,
                        self.unreduced_primes,max_drivers=max_drivers)
                elif driver_method == 'internal':
                    merger_drivers = sm_doi.internal_drivers(motif_merger,
                        self.unreduced_primes,max_drivers=max_drivers)

                drivers += [x for x in merger_drivers if not x in drivers]

        drivers = sorted(drivers,key=lambda x: len(x))

        # Next, we remove redundant control sets.
        # In the minimal_history scheme, a set x is redundant if there is a y
        # that contains all control sets of x in the same order.
        # In all other schemes, x is redundant if it is a subset of some y.
        nonredundant_drivers = []
        if target_method == 'history' and (driver_method == 'minimal' or driver_method == 'GRASP'):
            nonredundant_drivers = []
            for x in drivers: # drivers is sorted from fewest timesteps to most
                add_x = True
                for y in nonredundant_drivers: # note: len(y) <= len(x)
                    for offset in range(len(x)-len(y)+1):
                        for i in range(len(y)):
                            # want to check that x[i+offset] is a special case of y[i]
                            for j in range(len(y)):
                                xset = set([frozenset(tuple(xdict.items())) for xdict in x[i+offset]])
                                yset = set([frozenset(tuple(xdict.items())) for xdict in y[i]])
                                if yset <= xset:
                                    add_x = False
                                    break
                            if not add_x: break
                        if not add_x: break
                    if not add_x: break
                if add_x:
                    nonredundant_drivers.append(x)

        else:
            for i in range(len(drivers)):
                use_i = True
                for j in range(i): # i > j
                    if drivers[i].items() >= drivers[j].items():
                        use_i = False
                        break
                if use_i: nonredundant_drivers.append(drivers[i])

        return nonredundant_drivers

def build_succession_diagram(primes, fixed=None, motif_history=None, diagram=None, merge_equivalent_motifs=True,max_simulate_size=20,max_simulate_size_vc=None,prioritize_source_motifs=True,max_stable_motifs=10000,max_in_degree=float('inf'),MPBN_update=False):
    """Recursively construct a succession diagram from the input update rules.
    Generally, it is preferable to construct this from within the AttractorRepertoire
    class (using, e.g., AttractorRepertoire.from_primes).

    Parameters
    ----------
    primes : pyboolnet primes dictionary
        Update rules.
    fixed : partial state dictionary
        Used only for recursion. Specifies nodes to be fixed in the next reduced
        network to be added to the diagram.
    motif_history : list of partial state dictionaries
        Used only for recursion. Specifies stable motif history for the next reduced
        network to be added to the diagram.
    diagram : SuccessionDiagram
        Used only for recursion. The SuccessionDiagram object that is under
        construction.
    merge_equivalent_motifs : bool
        If False, equivalent reduced networks have their data recomputed and
        copied. Making this False is only recommended if the succession diagram
        must be represented in a form that has no feedforward loops; making this
        True provides large computational advantages, both in terms of speed and
        memory usage (the default is True).
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
        Whether MBPN update is used instead of general asynchronous update
        (the default is False).

    Returns
    -------
    SuccessionDiagram
        The succession diagram for the input update rules.

    """

    if diagram is None:
        diagram = SuccessionDiagram()

    if fixed is None:
        fixed = {}

    myMotifReductionToCopy = diagram.find_equivalent_reduction(fixed)

    if myMotifReductionToCopy is None:
        myMotifReduction=sm_reduction.MotifReduction(motif_history,fixed.copy(),primes,max_simulate_size=max_simulate_size,max_simulate_size_vc=max_simulate_size_vc,prioritize_source_motifs=prioritize_source_motifs,max_stable_motifs=max_stable_motifs,MPBN_update=MPBN_update)
    else:
        myMotifReduction = deepcopy(myMotifReductionToCopy)
        myMotifReduction.motif_history = motif_history.copy()
        myMotifReduction.merged_history_permutations = []

    diagram.add_motif_reduction(myMotifReduction)

    # Prioritize source nodes
    if myMotifReduction.merged_source_motifs is None:
        stable_motif_list = myMotifReduction.stable_motifs
    else:
        stable_motif_list = myMotifReduction.merged_source_motifs

    for sm in stable_motif_list:
        if merge_equivalent_motifs:
            perm_index,perm = diagram.find_motif_permutation(myMotifReduction.motif_history+[sm])
            if not perm_index is None:
                diagram.add_motif_permutation(perm_index,perm)
        if not merge_equivalent_motifs or perm_index is None:
            np,fixed2 = sm_reduction.reduce_primes(sm,primes)
            fixed3 = fixed.copy()
            fixed3.update(fixed2)
            diagram = build_succession_diagram(np,fixed3,myMotifReduction.motif_history+[sm],
                diagram, merge_equivalent_motifs=merge_equivalent_motifs,
                max_simulate_size=max_simulate_size,
                max_simulate_size_vc=max_simulate_size_vc,
                max_stable_motifs=max_stable_motifs,
                prioritize_source_motifs=prioritize_source_motifs,
                max_in_degree=max_in_degree,
                MPBN_update=MPBN_update)
    return diagram
