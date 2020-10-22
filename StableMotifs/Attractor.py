import StableMotifs.Reduction as sm_reduction


class Attractor:
    """Stores attractor data for a reduced network. Automatically initialized by
    the AttractorRepertoire class.

    Parameters
    ----------
    reduction : Reduction.MotifReduction
        Motif reduction to use as the representative (see attributes).
    reduction_attractor_id : int
        Reduction id to use for the representative (see attributes)

    Attributes
    ----------
    logically_fixed_nodes : partial state dictionary
        The nodes that are fixed by percolation on the expanded network (i.e.,
        not by up-stream oscillations)
    representative : Reduction.MotifReduction, int tuple
        Entry 0 is a maximally reduced Reduction.MotifReduction object that
        contains the attractor. In general, other such objects conain the
        attractor, but they will correspond to equivalent reduced networks.
        Entry 1 is a unique identifier number (integer) for the attractor within
        the reduced network; this is necessary in cases when a fully reduced n
        etwork contains multiple (complex) attractors.
    reductions : list of Reduction.MotifReduction
        Maximally reduced MotifReductions that contain the attractor.
    attractor_dict : dictionary
        a dictionary describing the node states in the attractor according to
        the following key -
             1 variable is "ON"
             0 variable is "OFF"
             X variable is known to oscillate
             ? at least one such variable must oscillate
             ! the attractor may be false; if it is genuine, at least one such
               variable must oscillate
    stg : networkx.DiGraph
        The state transition graph corresponding to the attractor (if computed)
    fixed_nodes : type
        Description of attribute `fixed_nodes`.
    oscillation_fixed_nodes : type
        Description of attribute `oscillation_fixed_nodes`.
    reduced_primes : PyBoolNet primes dictionary
        Update rules for the maximally reduced network that contains the attractor.
    n_unfixed : int
        Number of nodes that are not logically fixed.
    size : type
        Description of attribute `size`.
    size_lower_bound : type
        Description of attribute `size_lower_bound`.
    size_upper_bound : type
        Description of attribute `size_upper_bound`.
    explored : bool
        True if all attractor states and transitions are explicitly computed.
    guaranteed : bool
        True if and only if the attractor is known to be genuine. If False, the
        attractor may not actually be stable.


    """

    def __init__(self, reduction, reduction_attractor_id):
        self.logically_fixed_nodes = reduction.logically_fixed_nodes
        self.representative = (reduction,reduction_attractor_id)
        self.reductions = [reduction]
        self.attractor_dict = reduction.attractor_dict_list[reduction_attractor_id]

        if reduction.no_motif_attractors is not None:
            self.stg = reduction.partial_STG.subgraph(reduction.no_motif_attractors[reduction_attractor_id])
        else:
            self.stg = None

        self.fixed_nodes = {k:v for k,v in self.attractor_dict.items() if v in [0,1]}

        # if reduction.fixed_rspace_nodes is not None:
        #     self.fixed_nodes.update(reduction.fixed_rspace_nodes)

        self.oscillation_fixed_nodes = {k:v for k,v in self.fixed_nodes.items() if k not in self.logically_fixed_nodes}

        if reduction.rspace_update_primes is not None:
            self.reduced_primes = reduction.rspace_update_primes
        else:
            self.reduced_primes = reduction.reduced_primes

        self.n_unfixed = len(self.attractor_dict) - len(self.fixed_nodes)

        if self.stg is not None:
            self.size = len(self.stg)
            self.size_lower_bound = self.size
            self.size_upper_bound = self.size
        else:
            if self.n_unfixed <= 1:
                self.size = self.n_unfixed + 1
                self.size_lower_bound = self.size
                self.size_upper_bound = self.size
            else:
                self.size = None
                self.size_lower_bound = 2
                self.size_upper_bound = 2**(self.n_unfixed)

        if '?' in self.attractor_dict.values():
            self.explored = False
            self.guaranteed = True
        elif '!' in self.attractor_dict.values():
            self.explored = False
            self.guaranteed = False
        else: # X's, 1's, 0's
            self.explored = True
            self.guaranteed = True

        self.reductions = [reduction]

    def add_reduction(self, reduction):
        """Add a reduction to the attractor. Does not check for compatibility.

        Parameters
        ----------
        reduction : Reduction.MotifReduction
            Motif reduction that also contains the attractor.

        """
        self.reductions.append(reduction)
