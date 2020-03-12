from StableMotifs.Reduction import MotifReduction, reduce_primes
from StableMotifs.Format import pretty_print_prime_rules
from StableMotifs.DomainOfInfluence import find_internal_motif_drivers
class SuccessionDiagram:
    """
    Representation of a succession diagram of a Boolean system

    Variables:
    MotifReductionList - a list of MotifReductions (see Reduction.py)

    Functions:
    __init__(self)
    find_motif_permutation(self,motif_history)
    add_motif_permutation(self,reduction_index,permutation)
    add_motif_reduction(self,motif_reduction)
    summary(self,terminal_keys=None,show_original_rules=False) - prints a summary of the succession diagram to screen
    attractor_candidate_summary(self) - prints a summary of found or potential attractors
    """

    def __init__(self):
        self.MotifReductionList = []
        self.reduction_permutations = []
        self.attractor_fixed_nodes_list = []
        self.attractor_reduced_primes_list = []
        self.attractor_guaranteed_list = []
        self.reduced_complex_attractor_list = []
        self.unreduced_primes = None
    def find_motif_permutation(self,motif_history):
        for i,mr in enumerate(self.MotifReductionList):
            if len(mr.motif_history) == len(motif_history):
                if all([x in mr.motif_history for x in motif_history]):
                    permutation = []
                    for x in motif_history:
                        permutation.append(mr.motif_history.index(x))
                    return i,permutation # We already have an equivalent motif in the list
        return None,None

    def add_motif_permutation(self,reduction_index,permutation):
        self.reduction_permutations[reduction_index].append(permutation)

    def add_motif_reduction(self,motif_reduction,merge_equivalent_motifs=True):
        if self.MotifReductionList == []:
            self.unreduced_primes = motif_reduction.reduced_primes
        self.MotifReductionList.append(motif_reduction)
        self.reduction_permutations.append([])
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
        for motif_reduction, motif_permutations in zip(self.MotifReductionList,self.reduction_permutations):
            if terminal_keys is None or motif_reduction.terminal in terminal_keys:
                print("__________________")
                motif_reduction.summary(show_original_rules=show_original_rules)
                if len(motif_permutations) > 0:
                    print()
                    print("The following motif_history permutation(s) have been merged into this branch:")
                    for x in motif_permutations: print(x)

    def motif_sequence_to_reductions(self,target_motif_reductions):
        """
        Input:
        target_motif_reductions - a list of MotifReductions that we want to reprogram to (we want to reach any, not necessarily all)

        Output:
        stable motif sequences that achieve the desired reprogramming
        """
        target_motif_mergers = [] # we will think of members as unordered for now; might consider orders later
        for reduction in target_motif_reductions:
            target_motif_mergers.append({k:v for d in reduction.motif_history for k,v in d.items()})
        return target_motif_mergers

    def find_reductions_with_states(self,logically_fixed):
        # NOTE: This finds all reductions, not just those closest to the root
        target_motif_reductions = []
        for reduction in self.MotifReductionList:
            if logically_fixed.items() <= reduction.logically_fixed_nodes.items():
                target_motif_reductions.append(reduction)
        return target_motif_reductions

    def reprogram_to_trap_spaces(self,logically_fixed,max_internal_drivers=None,max_external_drivers=None):
        # TODO: consider motifs separately for better scaling.
        # Let all but one motif lock-in, then look for drivers for the last.
        # Do this for each motif, then we don't have to worry about order, but
        # we're not oversampling the driver space by such an obscene amount. 
        target_motif_reductions = self.find_reductions_with_states(logically_fixed)
        target_motif_mergers = self.motif_sequence_to_reductions(target_motif_reductions)

        drivers = []
        for motif_merger in target_motif_mergers:
            drivers += find_internal_motif_drivers(motif_merger,self.unreduced_primes)
        return sorted({k:v for d in drivers for k,v in d.items()},key=lambda x:len(x))

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
    diagram.add_motif_reduction(myMotifReduction,merge_equivalent_motifs=merge_equivalent_motifs)

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
            diagram = build_succession_diagram(np,fixed3,myMotifReduction.motif_history+[sm],diagram, merge_equivalent_motifs=merge_equivalent_motifs)
    return diagram
