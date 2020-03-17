import networkx as nx
import itertools as it

from StableMotifs.Reduction import MotifReduction, reduce_primes
from StableMotifs.Format import pretty_print_prime_rules
from StableMotifs.DomainOfInfluence import internal_drivers, minimal_drivers

class SuccessionDiagram:
    """
    Representation of a succession diagram of a Boolean system

    Variables:
    motif_reduction_list - a list of MotifReductions (see Reduction.py)

    Functions:
    __init__(self)
    find_motif_permutation(self,motif_history)
    add_motif_permutation(self,reduction_index,permutation)
    add_motif_reduction(self,motif_reduction)
    summary(self,terminal_keys=None,show_original_rules=False) - prints a summary of the succession diagram to screen
    attractor_candidate_summary(self) - prints a summary of found or potential attractors
    """

    def __init__(self):
        self.motif_reduction_list = []
        self.digraph = nx.DiGraph()
        self.attractor_fixed_nodes_list = []
        self.attractor_reduced_primes_list = []
        self.attractor_guaranteed_list = []
        self.reduced_complex_attractor_list = []
        self.unreduced_primes = None
    def find_motif_permutation(self,motif_history):
        for i,mr in enumerate(self.motif_reduction_list):
            if len(mr.motif_history) == len(motif_history):
                if all([x in mr.motif_history for x in motif_history]):
                    permutation = []
                    for x in motif_history:
                        permutation.append(mr.motif_history.index(x))
                    return i,permutation # We already have an equivalent motif in the list
        return None,None

    def add_motif_permutation(self,reduction_index,permutation):
        self.motif_reduction_list[reduction_index].merged_history_perumutations.append(permutation)
        for child in nx.topological_sort(self.digraph):
            for parent in self.digraph.predecessors(child):
                for parent_perm,child_perm in it.product(self.motif_reduction_list[parent].merged_history_perumutations,self.motif_reduction_list[child].merged_history_perumutations):
                    new_perm = child_perm.copy()
                    for i,p in enumerate(parent_perm):
                        new_perm[i] = child_perm[p]
                    if not new_perm in self.motif_reduction_list[child].merged_history_perumutations:
                        self.motif_reduction_list[child].merged_history_perumutations.append(new_perm)

    def add_motif_reduction(self,motif_reduction):
        if self.motif_reduction_list == []:
            self.unreduced_primes = motif_reduction.reduced_primes

        # note: N is computed BEFORE the new reduction is added,
        # so it will be the reduction index AFTER the reduction is added.
        N = len(self.motif_reduction_list)
        new_set = set([frozenset(tuple(x.items())) for x in motif_reduction.motif_history])
        for i,reduction in enumerate(self.motif_reduction_list):
            old_set = set([frozenset(tuple(x.items())) for x in reduction.motif_history])

            # see if we're adding a parent of an existing reduction
            if len(old_set) == len(new_set) + 1:
                diff_set = old_set - new_set
                if len(diff_set) == 1:
                    missing_motif = dict(diff_set.pop())#[dict(s) for s in old_set - new_set][0]
                    if missing_motif in motif_reduction.stable_motifs:
                        self.digraph.add_edge(N,i)
            # see if we're adding a child of an existing reduction
            elif len(old_set) == len(new_set) - 1:
                diff_set = new_set - old_set
                if len(diff_set) == 1:
                    missing_motif = dict(diff_set.pop())#[dict(s) for s in new_set - old_set][0]
                    if missing_motif in reduction.stable_motifs:
                        self.digraph.add_edge(i,N)

        self.motif_reduction_list.append(motif_reduction)
        self.add_motif_permutation(N,list(range(len(motif_reduction.motif_history))))
        if not motif_reduction.terminal == "no":
            if not motif_reduction.logically_fixed_nodes in self.attractor_fixed_nodes_list:
                self.attractor_fixed_nodes_list.append(motif_reduction.logically_fixed_nodes)
                self.attractor_reduced_primes_list.append(motif_reduction.reduced_primes)
                self.attractor_guaranteed_list.append(motif_reduction.terminal)
                self.reduced_complex_attractor_list.append(motif_reduction.no_motif_attractors)


    def attractor_candidate_summary(self):
        print("Found", len([x for x in self.attractor_guaranteed_list if x == "yes"]), "guaranteed attractor space(s) and",
            len([x for x in self.attractor_guaranteed_list if x == "possible"]), "possible attractor space(s).")
        print("Found", len([x for x in self.attractor_reduced_primes_list if len(x)==0]), "steady state(s) and",
            sum([len(x) for x in self.reduced_complex_attractor_list if not x is None]), "complex attractor(s) in the guaranteed attractor space(s).")
        for fn,rp,tr,at in zip(self.attractor_fixed_nodes_list,
                               self.attractor_reduced_primes_list,
                               self.attractor_guaranteed_list,
                               self.reduced_complex_attractor_list):
            print("__________________")
            if tr == "possible":
                print("Space May Contain Attractor")
                print()
            else:
                print("Space Guaranteed to Contain Attractor(s)")
            print("Logically Fixed Nodes:",{k:v for k,v in sorted(fn.items())})
            print()
            if len(rp) > 0:
                print("Reduced Rules:")
                pretty_print_prime_rules(rp)
                if not at is None:
                    print()
                    print("Complex Attractors in Reduced Network (Alphabetical Node Ordering):")
                    for x in at:
                        print(x)
            else:
                print("No Free Nodes Remain.")

    def summary(self,terminal_keys=None,show_original_rules=True,show_explicit_permutations=False):
        for motif_reduction in self.motif_reduction_list:
            if terminal_keys is None or motif_reduction.terminal in terminal_keys:
                print("__________________")
                motif_reduction.summary(show_original_rules=show_original_rules,show_explicit_permutations=show_explicit_permutations)

    def merge_reduction_motifs(self,target_reductions):
        """
        Input:
        target_reductions - a list of MotifReductions that we want to reprogram to (we want to reach any, not necessarily all)

        Output:
        stable motif sequences that achieve the desired reprogramming
        """
        target_motif_mergers = [] # we will think of members as unordered for now; might consider orders later
        for reduction in target_reductions:
            target_motif_mergers.append({k:v for d in reduction.motif_history for k,v in d.items()})
        return target_motif_mergers

    def reductions_indices_with_states(self,logically_fixed,optimize=True):
        if not optimize:
            # NOTE: This finds all reductions, not just those closest to the root
            target_indicies = []
            for i,reduction in enumerate(self.motif_reduction_list):
                if logically_fixed.items() <= reduction.logically_fixed_nodes.items():
                    target_indicies.append(i)
            return target_indicies
        else:
            # The unoptimized result
            target_indices_unoptimized = self.reductions_indices_with_states(logically_fixed,optimize=False)
            target_indices_all = target_indices_unoptimized.copy()

            # Add nodes that inevitably reach nodes in the unoptimized result
            nodes_to_consider = [x for x in self.digraph if (not x in target_indices_all
                and self.motif_reduction_list[x].terminal == "no")]
            for i in nodes_to_consider:
                bad_sinks = [x for x in nx.descendants(self.digraph,i) | set([i]) if (
                    not self.motif_reduction_list[x].terminal == "no" # is a sink and . . .
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


    def reduction_drivers(self,target_index,method='internal',max_drivers=None):
        methods = ['internal','minimal']
        assert method in methods, ' '.join(["method argument of reduction_drivers must be among",str(methods)])
        drivers = []
        for path in nx.all_simple_paths(self.digraph,0,target_index):
            path_motif_history=[]
            path_drivers = []
            for ind in path:
                if ind == 0:
                    ind_prev = ind
                    continue
                path_motif_history += [x for x in self.motif_reduction_list[ind].motif_history if not x in path_motif_history]

                if method == 'internal':
                    history_drivers = internal_drivers(path_motif_history[-1],
                        self.motif_reduction_list[ind_prev].reduced_primes,
                        max_drivers=max_drivers)
                elif method == 'minimal':
                    history_drivers = minimal_drivers(path_motif_history[-1],
                        self.motif_reduction_list[ind_prev].reduced_primes,
                        max_drivers=max_drivers)

                path_drivers.append(history_drivers)

                ind_prev = ind
            if method == 'internal':
                # merge control sets along the path
                for control_sequence in it.product(*path_drivers):
                    control_set = {k:v for x in control_sequence for k,v in x.items()}
                    if not control_set in drivers:
                        drivers.append(control_set)
            elif method == 'minimal':
                drivers.append(path_drivers)
        return drivers

    def reprogram_to_trap_spaces(self,logically_fixed,max_drivers=None,method='history'):
        """
        Find driver sets that lead to the node states specified by logically_fixed

        Inputs:
        logically_fixed - state dictionary specifying the control target
        max_drivers - the maximum number of driver nodes to attempt when looking
                      for stable motif driver nodes before specifying the entire
                      stable motif as part fo the driver set
        method - One of the following: 'history', 'merge', 'minimal',
                 'minimal_history'; specifies the reprogramming method to use
                 (see below for further details)

        Output:
        nonredundant_drivers - control strategies found; interpretation depends
                               on method selected (see below)


        - Methods -
        history:
        Finds all shortest stable motif histories that result in the target node states
        being logically fixed. Each stable motif is searched for internal driver nodes.
        The resulting internal drivers are combined into a single  control set. The
        return value consists of all such control sets for all  stable motif histories.
        Each control set eventually becomes self-sustaining.

        minimal_history:
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

        merge:
        Finds all shortest stable motif histories that result in the target node states
        being logically fixed. All node states in the motifs in the history are merged
        into a stable module dictionary. This is then searched for internal driver
        nodes. Each element of the return value is a dictionary corresponding to a
        control set. Each control set eventually becomes self-sustaining.

        minimal:
        Similar to the merge method, except the search for drivers is conducted over
        all nodes, not just those internal to the merged stable module. Furthermore,
        the search is truncated when a control set is found such that the search does
        not proceed to driver sets larger than the smallest found. Each element of
        the return value is a dictionary corresponding to a control set. The control
        sets are only guaranteed to result in activation of the target if they are
        temporary interventions.

        """

        #methods = ['history','merge','minimal_history','minimal_merge']
        methods = ['history','merge','minimal','minimal_history']
        assert method in methods, ' '.join(["method argument of reprogram_to_trap_spaces must be among",str(methods)])

        drivers = []
        #target_indices_all = self.reductions_indices_with_states(logically_fixed,optimize=False)
        target_indices = self.reductions_indices_with_states(logically_fixed)
        # # We only look for drivers if the reduction doesn't have any ancestors that would also work
        # target_indices = []
        # for target_index in target_indices_all:
        #     if set(nx.ancestors(self.digraph,target_index)) & set(target_indices_all) == set():
        #         target_indices.append(target_index)

        if method == 'history':
            for target_index in target_indices:
                drivers += self.reduction_drivers(target_index,max_drivers=max_drivers)
        elif method == 'minimal_history':
            for target_index in target_indices:
                drivers += self.reduction_drivers(target_index,max_drivers=max_drivers,method='minimal')
        elif method == 'merge':
            for target_index in target_indices:
                target_history = self.motif_reduction_list[target_index].motif_history
                motif_merger = {k:v for d in target_history for k,v in d.items()}
                merger_drivers = internal_drivers(motif_merger,
                    self.unreduced_primes,max_drivers=max_drivers)
                drivers += [x for x in merger_drivers if not x in drivers]
        elif method == 'minimal':
            for target_index in target_indices:
                target_history = self.motif_reduction_list[target_index].motif_history
                motif_merger = {k:v for d in target_history for k,v in d.items()}

                # Because we are potentially dealing with external drivers, we
                # want to make sure the external drivers do not interfere with
                # the motif's ability to drive downstream targets, or the
                motif_merger.update(logically_fixed)
                merger_drivers = minimal_drivers(motif_merger,
                    self.unreduced_primes,max_drivers=max_drivers)
                drivers += [x for x in merger_drivers if not x in drivers]

        drivers = sorted(drivers,key=lambda x: len(x))

        # Next, we remove redundant control sets.
        # In the minimal_history scheme, a set x is redundant if there is a y
        # that contains all control sets of x in the same order.
        # In all other schemes, x is redundant if it is a subset of some y.
        nonredundant_drivers = []
        if method == 'minimal_history':
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

def build_succession_diagram(primes, fixed=None, motif_history=None, diagram=None, merge_equivalent_motifs=True,search_partial_STGs=True,prioritize_source_motifs=True):
    """
    Constructs a succession diagram recursively from the rules specified by primes

    Inputs:
    primes - PyBoolNet primes dictionary specifying the Boolean update rules

    Inputs used in recursion only:
    fixed - dictionary with node names as keys and fixed node states as values
    motif_history: list of stable motifs that have been "locked in" so far in the recursion
    diagram - the succession diagram being constructed by the recursion

    Outputs:
    diagram - SuccessionDiagram object describing the succession diagram for the system
    """
    if fixed is None:
        fixed = {}
    myMotifReduction=MotifReduction(motif_history,fixed.copy(),primes,search_partial_STGs=search_partial_STGs,prioritize_source_motifs=prioritize_source_motifs)
    if diagram is None:
        diagram = SuccessionDiagram()
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
            np,fixed2 = reduce_primes(sm,primes)
            fixed3 = fixed.copy()
            fixed3.update(fixed2)
            diagram = build_succession_diagram(np,fixed3,myMotifReduction.motif_history+[sm],
                diagram, merge_equivalent_motifs=merge_equivalent_motifs,
                search_partial_STGs=search_partial_STGs,
                prioritize_source_motifs=prioritize_source_motifs)
    return diagram
