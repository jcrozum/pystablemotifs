import networkx as nx
import itertools as it

from StableMotifs.Reduction import MotifReduction, reduce_primes
from StableMotifs.Format import pretty_print_prime_rules
from StableMotifs.DomainOfInfluence import find_internal_motif_drivers

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
                missing_motif = [dict(s) for s in old_set - new_set][0]
                if missing_motif in motif_reduction.stable_motifs:
                    self.digraph.add_edge(N,i)
            # see if we're adding a child of an existing reduction
            elif len(old_set) == len(new_set) - 1:
                missing_motif = [dict(s) for s in new_set - old_set][0]
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

    def summary(self,terminal_keys=None,show_original_rules=True):
        for motif_reduction in self.motif_reduction_list:
            if terminal_keys is None or motif_reduction.terminal in terminal_keys:
                print("__________________")
                motif_reduction.summary(show_original_rules=show_original_rules)

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

    def find_reductions_with_states(self,logically_fixed):
        # NOTE: This finds all reductions, not just those closest to the root
        target_reductions = []
        target_indicies = []
        for i,reduction in enumerate(self.motif_reduction_list):
            if logically_fixed.items() <= reduction.logically_fixed_nodes.items():
                target_indicies.append(i)
                target_reductions.append(reduction)
        return target_reductions,target_indicies

    def find_reduction_drivers(self,target_index,max_drivers=None):
        drivers = []
        for path in nx.all_simple_paths(self.digraph,0,target_index):
            path_motif_history=[]
            path_drivers = []
            for ind in path:
                if ind == 0:
                    ind_prev = ind
                    continue
                path_motif_history += [x for x in self.motif_reduction_list[ind].motif_history if not x in path_motif_history]
                path_drivers.append(find_internal_motif_drivers(path_motif_history[-1],
                    self.motif_reduction_list[ind_prev].reduced_primes,
                    max_drivers=max_drivers))
                ind_prev = ind
            for control_sequence in it.product(*path_drivers):
                drivers.append({k:v for x in control_sequence for k,v in x.items()})
        return drivers

    def reprogram_to_trap_spaces(self,logically_fixed,max_drivers=None,method='history'):
        # TODO: consider motifs separately for better scaling. Maybe find middle
        # ground with motif merger approach.
        methods = ['history','merge']
        assert method in methods, ' '.join(["method argument of reprogram_to_trap_spaces must be among",str(methods)])

        drivers = []
        target_reductions,target_indices = self.find_reductions_with_states(logically_fixed)

        if method == 'history':
            for target_index in target_indices:
                if set(nx.ancestors(self.digraph,target_index)) & set(target_indices) == set():
                    drivers += self.find_reduction_drivers(target_index,max_drivers=max_drivers)
                else: print("ignoring",target_index)
        elif method == 'merge':
            target_motif_mergers = self.merge_reduction_motifs(target_reductions)
            for motif_merger in target_motif_mergers:
                merger_drivers = find_internal_motif_drivers(motif_merger,
                    self.unreduced_primes,max_drivers=max_drivers)
                drivers += [x for x in merger_drivers if not x in drivers]

        drivers = sorted(drivers,key=lambda x: len(x))

        nonredundant_drivers = []
        # remove redundant control sets
        for i in range(len(drivers)):
            use_i = True
            for j in range(i): # i > j
                if drivers[i].items() >= drivers[j].items():
                    use_i = False
                    break
            if use_i: nonredundant_drivers.append(drivers[i])

        return nonredundant_drivers

def build_succession_diagram(primes, fixed=None, motif_history=None, diagram=None, merge_equivalent_motifs=True):
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
    myMotifReduction=MotifReduction(motif_history,fixed.copy(),primes)
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
                diagram, merge_equivalent_motifs=merge_equivalent_motifs)
    return diagram
