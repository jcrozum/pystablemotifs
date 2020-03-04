from StableMotifs.Reduction import MotifReduction, reduce_primes
from StableMotifs.Format import pretty_print_prime_rules
class SuccessionDiagram:
    """
    Representation of a succession diagram of a Boolean system

    Variables:
    MotifReductionList - a list of MotifReductions (see Reduction.py)

    Functions:
    __init__(self)
    add_motif_reduction(self,motif_reduction)
    summary(self,terminal_keys=None) - prints a summary of the succession diagram to screen
    attractor_candidates(self) - prints a summary of found or potential attractors
    """

    def __init__(self):
        self.MotifReductionList = []
    def add_motif_reduction(self,motif_reduction):
        self.MotifReductionList.append(motif_reduction)
    def attractor_candidates(self):
        fixed_nodes_list = []
        reduced_primes_list = []
        terminal_list = []
        for motif_reduction in self.MotifReductionList:
            if not motif_reduction.terminal == "no":
                if not motif_reduction.logically_fixed_nodes in fixed_nodes_list:
                    fixed_nodes_list.append(motif_reduction.logically_fixed_nodes)
                    reduced_primes_list.append(motif_reduction.reduced_primes)
                    terminal_list.append(motif_reduction.terminal)
        for fn,rp,tr in zip(fixed_nodes_list,reduced_primes_list,terminal_list):
            print("__________________")
            if tr == "possible":
                print("Space May Contain Attractor")
                print()
            else:
                print("Space Guaranteed to Contain Attractor")
            print("Logically Fixed Nodes:",{k:v for k,v in sorted(fn.items())})
            print()
            print("Reduced Rules:")
            pretty_print_prime_rules(rp)

    def summary(self,terminal_keys=None):
        for motif_reduction in self.MotifReductionList:
            if terminal_keys is None or motif_reduction.terminal in terminal_keys:
                print("__________________")
                motif_reduction.summary()

def build_succession_diagram(primes, fixed=None, motif_history=None, diagram=None):
    """
    Constructs a succession diagram recursively from the rules specified by primes

    Inputs:
    primes - PyBoolNet primes dictionary specifying the Boolean update rules

    Inputs used in recursion only:
    fixed - dictionary with node names as keys and fixed node states as values
    motif_hisory: list of stable motifs that have been "locked in" so far in the recursion
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

    for sm in myMotifReduction.stable_motifs:
        np,fixed2 = reduce_primes(sm,primes)
        fixed3 = fixed.copy()
        fixed3.update(fixed2)
        diagram = build_succession_diagram(np,fixed3,myMotifReduction.motif_history+[sm],diagram)
    return diagram
