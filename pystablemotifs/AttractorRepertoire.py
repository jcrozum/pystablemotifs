import pystablemotifs.succession as sm_succession
import pystablemotifs.Attractor as Attractor
import itertools as it
import networkx as nx
import pyboolnet.prime_implicants

class AttractorRepertoire:
    """The class that stores information about attractors. Initialize using
    either from_primes or from_succession_diagram.

    Attributes
    ----------
    succession_diagram : succession.SuccessionDiagram
        Succession diagram summarizing the stable motif structure of the model.
    attractors : list of Attractor.Attractor
        List of (possible) attractors in the model.
    reduction_attractors : dictionary
        A dictionary with integer keys that correspond to the
        succession_diagram.digraph nodes. The dictionary values are lists of
        Attractor.Attractor objects that correspond to attractors that exist in
        the region of statespace corresponding to the reduced network
        represented by the key in the succession diagram.
    fewest_attractors : int
        A lower bound on the number of attractors in the model.
    most_attractors : int
        An upper bound on the number of attractors in the model.
    primes : pyboolnet primes dictionary
        The model rules.
    succession_digraph : networkx digraph
        Networkx digraph representation of the succession_diagram object. If
        AttractorRepertoire.simplify_diagram, it is equivalent to
        AttractorRepertoire.succession_diagram.digraph. Otherwise, several of its
        nodes may be contracted (depending on input parameters).
    attractor_equivalence_classes : list
        List of attractor equivalence classes. Each item is a dictionary with keys
        'states', 'attractors', and 'reductions'. The 'states' value is a dictionary
        of variable values that all attractors in the class share. The 'attractors'
        value is a list of Attractor objects (i.e., a sublist of self.attractors);
        all attractors in this list have all relevant nodes equivalently characterized.
        The 'reductions' value is a list of reduction_attractor keys that collectively
        contain all the attractors in the class (and therefore cannot differ in any
        relevant node).
    relevant_nodes : list
        List of nodes that are "relevant", i.e., if trap spaces differ in the values
        of these variables, then the corresponding succession diagram nodes and
        attractors will not be merged.

    """

    def __init__(self):
        self.succession_diagram = None
        self.attractors = []
        self.reduction_attractors = {}
        self.fewest_attractors = None
        self.most_attractors = None
        self.primes = None

        self.succession_digraph = None
        self.attractor_equivalence_classes = None
        self.relevant_nodes = None

    @classmethod
    def from_primes(cls,primes,max_simulate_size=20,max_simulate_size_vc=None,max_stable_motifs=10000,max_in_degree=float('inf'),MPBN_update=False):
        """Build the succession diagram and attractor repertoire from pyboolnet
        formatted update rules rules.

        Parameters
        ----------
        primes : pyboolnet primes dictionary
            The model rules.
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
            (see L Pauleve, J Kolcak, T Chatain, S Haar, "Reconciling qualitative, abstract, and scalable modeling of biological networks." Nat. Com. vol. 11, no. 4256 (2020))
            (the default is False).

        Returns
        -------
        AttractorRepertoire
            AttractorRepertoire object for the input primes.

        """
        x = cls()
        x.primes = primes
        x.analyze_system(primes,max_simulate_size=max_simulate_size,max_simulate_size_vc=max_simulate_size_vc,max_stable_motifs=max_stable_motifs,max_in_degree=max_in_degree,MPBN_update=MPBN_update)
        x.simplify_diagram([], merge_equivalent_reductions = False)
        return x

    @classmethod
    def from_succession_diagram(cls,succession_diagram):
        """Build the succession diagram and attractor repertoire from a
        precomputed succession diagram.

        Parameters
        ----------
        succession_diagram : succession.SuccessionDiagram
            Succession diagram summarizing the stable motif structure of the model.

        Returns
        -------
        AttractorRepertoire
            AttractorRepertoire object for the input succession diagram.

        """
        x = cls()
        x.succession_diagram = succession_diagram
        x.primes = succession_diagram.unreduced_primes
        x._get_attractors_from_succession_diagram()
        x._count_attractors()
        x.simplify_diagram([], merge_equivalent_reductions = False)
        return x


    def _get_attractors_from_succession_diagram(self):
        """Extract attractors from the succession diagram of the model.

        """
        for ri, reduction in self.succession_diagram.motif_reduction_dict.items():
            if reduction.terminal == "no": continue

            self.reduction_attractors[ri] = []

            duplicate = False # have we added the attractors for this reduciton yet?
            for attractor in self.attractors:
                if attractor.logically_fixed_nodes == reduction.logically_fixed_nodes:
                    duplicate = True
                    attractor.add_reduction(reduction)
                    self.reduction_attractors[ri].append(attractor)

            if duplicate: continue

            for id,att in enumerate(reduction.attractor_dict_list):
                new_attractor = Attractor(reduction,id)
                self.attractors.append(new_attractor)
                self.reduction_attractors[ri].append(new_attractor)

    def _count_attractors(self):
        """Place upper and lower bounds on the number of attractors.

        """
        self.fewest_attractors = 0
        self.most_attractors = 0
        for attractor in self.attractors:
            if attractor.guaranteed:
                self.fewest_attractors += 1
            # else: self.fewest_attractors += 0

            if attractor.explored:
                self.most_attractors += 1
            else:
                if attractor.representative[0].deletion_no_motif_attractors is not None:
                    self.most_attractors += len(attractor.representative[0].deletion_no_motif_attractors)
                else:
                    # ludicrously conservative upper bound; assumes STG is all 2-cycles
                    self.most_attractors += 2**(attractor.n_unfixed - 1)

    def analyze_system(self,primes,max_simulate_size=20,max_simulate_size_vc=None,max_stable_motifs=10000,max_in_degree=float('inf'),MPBN_update=False):
        """Build and process the succession diagram for the model.

        Parameters
        ----------
        primes : pyboolnet primes dictionary
            The model rules.
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

        """
        self.succession_diagram = sm_succession.build_succession_diagram(primes,max_simulate_size=max_simulate_size,max_simulate_size_vc=max_simulate_size_vc,max_stable_motifs=max_stable_motifs,max_in_degree=max_in_degree,MPBN_update=MPBN_update)
        self._get_attractors_from_succession_diagram()
        self._count_attractors()

    def simplify_diagram(self, projection_nodes, merge_equivalent_reductions = True, keep_only_projection_nodes = False, condense_simple_paths = False):
        """Simplify the succession diagram for the model. This is done in two ways.
        First, variables can be designated ignorable using the projection_nodes
        parameter. If keep_only_projection_nodes is False, these variables are
        ignorable, otherwise, all other nodes are ignorable. When
        merge_equivalent_reductions is True, all nodes of the succession diagram
        that correspond to trap spaces whose fixed variables differ only in ignorable
        variables are contracted (in the graph theory sense). After this process,
        if condense_simple_paths is True, then all succession diagram nodes with
        in-degree equal to one are contracted with their parent node. This function
        constructs the succession_digraph and attractor_equivalence_classes attributes,
        which are described in the class documentation.

        Parameters
        ----------
        projection_nodes : list of variable names
            These nodes will be ignored if keep_only_projection_nodes is False
            (default); otherwise, all nodes except these will be ignored.
        merge_equivalent_reductions : bool
            Whether to contract succession diagram nodes whose reductions differ
            only in ignorable nodes.
        keep_only_projection_nodes : bool
            Whether projection_nodes specifies non-ignorable nodes.
        condense_simple_paths : bool
            Whether to contract nodes with in-degree one.
        """

        if not keep_only_projection_nodes:
            keep = set(self.primes.keys()) - set(projection_nodes)
            ignore = set(projection_nodes)
        else:
            keep = set(projection_nodes)
            ignore = set(self.primes.keys()) - set(projection_nodes)

        G = self.succession_diagram.digraph.copy()

        # Merge equivalent nodes
        if merge_equivalent_reductions:
            for u,v in it.combinations(self.succession_diagram.digraph.nodes(),2):
                if u not in G.nodes(): continue # we've already merged u

                ru = self.succession_diagram.motif_reduction_dict[u]
                rv = self.succession_diagram.motif_reduction_dict[v]
                if pyboolnet.prime_implicants.primes_are_equal(ru.reduced_primes,rv.reduced_primes):
                    rud = ru.logically_fixed_nodes
                    rvd = rv.logically_fixed_nodes
                    if all([rud[k]==rvd[k] for k in keep if k in rud and k in rvd]):
                        G = nx.contracted_nodes(G,u,v,self_loops=False)

        # If a node v has out degree = 1, then merge it into its only child u
        if condense_simple_paths:
            old_N = float('inf')
            while G.number_of_nodes() < old_N:
                old_N = G.number_of_nodes()
                for v in G.nodes():
                    if G.out_degree(v) == 1:
                        u = next(G.successors(v))
                        G = nx.contracted_nodes(G,u,v,self_loops=False)


        self.succession_digraph = G

        self.attractor_equivalence_classes = []
        for a in self.attractors:
            keys = set()
            for n in self.succession_digraph.nodes():
                if n not in self.reduction_attractors: continue
                for b in self.reduction_attractors[n]:
                    if all([b.attractor_dict[k] == a.attractor_dict[k] for k in keep]):
                        keys.add(n)
                        break

            merged = False
            for c in self.attractor_equivalence_classes:
                if c['states'].items() <= a.attractor_dict.items():
                    c['attractors'].append(a)
                    c['reductions'] |= keys
                    merged = True
            if not merged:
                self.attractor_equivalence_classes.append({'states':{k:v for k,v in a.attractor_dict.items() if k in keep},'attractors':[a], 'reductions':keys})

        self.relevant_nodes = keep

    def summary(self):
        """Prints a summary of the attractors to standard output.

        """
        if self.fewest_attractors == 0:
            print("Unable to properly count attractors.")
        elif self.fewest_attractors == self.most_attractors:
            if self.fewest_attractors == 1:
                print("There is 1 attractor.")
            else:
                print("There are",self.fewest_attractors,"attractors.")
        else:
            print("There are between",self.fewest_attractors,"and",self.most_attractors,"attractors.")
        for att in self.attractors:
            print(att.attractor_dict)
            print()

    def reprogram_to_trap_spaces(self,logically_fixed,
        target_method='history',driver_method='internal',
        max_drivers=None,GRASP_iterations=None,GRASP_score_override=None):
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
        return self.succession_diagram.reprogram_to_trap_spaces(logically_fixed,
            target_method=target_method,driver_method=driver_method,
            max_drivers=max_drivers,
            GRASP_iterations=GRASP_iterations,GRASP_score_override=GRASP_score_override)
