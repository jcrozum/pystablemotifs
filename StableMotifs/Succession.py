import networkx as nx
import itertools as it
import matplotlib
import matplotlib.pyplot as plt

import StableMotifs.Reduction as sm_reduction
import StableMotifs.Format as sm_format
import StableMotifs.DomainOfInfluence as sm_doi

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
        self.G_reduced_network_based = nx.DiGraph()
        self.G_reduced_network_based_labeled = nx.DiGraph()
        self.G_motif_based = nx.DiGraph()
        self.G_motif_based_labeled = nx.DiGraph()
        self.pos_reduced_network_based = dict()
        self.pos_motif_based = dict()
        self.attractor_fixed_nodes_list = []
        self.attractor_reduced_primes_list = []
        self.attractor_guaranteed_list = []
        self.reduced_complex_attractor_list = []
        self.deletion_attractor_list = []
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
        self.motif_reduction_list[reduction_index].merged_history_permutations.append(permutation)
        for child in nx.topological_sort(self.digraph):
            for parent in self.digraph.predecessors(child):
                for parent_perm,child_perm in it.product(self.motif_reduction_list[parent].merged_history_permutations,self.motif_reduction_list[child].merged_history_permutations):
                    new_perm = child_perm.copy()
                    for i,p in enumerate(parent_perm):
                        new_perm[i] = child_perm[p]
                    if not new_perm in self.motif_reduction_list[child].merged_history_permutations:
                        self.motif_reduction_list[child].merged_history_permutations.append(new_perm)

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
                self.deletion_attractor_list.append(motif_reduction.deletion_no_motif_attractors)

    def attractor_candidate_summary(self, show_reduced_rules = True):
        guaranteed_spaces = len([x for x in self.attractor_guaranteed_list if x == "yes"])
        possible_spaces = len([x for x in self.attractor_guaranteed_list if x == "possible"])
        steady_states =  len([x for x in self.attractor_reduced_primes_list if len(x)==0])
        found_complex_attractors = sum([len(x) for x in self.reduced_complex_attractor_list if x is not None])
        deletion_split_oscillations = sum([len(x) for x in self.deletion_attractor_list if x is not None and len(x)>1])
        deletion_lone_oscillations = sum([len(x) for x in self.deletion_attractor_list if x is not None and len(x)==1])

        lbound_oscillations = found_complex_attractors + deletion_lone_oscillations
        ubound_oscillations = lbound_oscillations + deletion_split_oscillations

        print("Found", guaranteed_spaces, "guaranteed attractor space(s) and",
            possible_spaces, "possible attractor space(s).")
        print("Found", steady_states, "steady state(s) and explored", found_complex_attractors,
            "complex attractor(s) in the guaranteed attractor space(s).")
        if ubound_oscillations == found_complex_attractors:
            print("There are no additional attractors.")
        elif deletion_split_oscillations == 0:
            print("There are exactly",deletion_lone_oscillations,"additional complex attractor(s) that were not fully explored.")
        else:
            print("There are between",lbound_oscillations,"and",ubound_oscillations,"complex attractors in total.")

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
                if show_reduced_rules:
                    print("Reduced Rules:")
                    sm_format.pretty_print_prime_rules(rp)
                else: print("Logically Unfixed Nodes:", sorted(rp.keys()))
                if not at is None:
                    print()
                    print("Complex Attractors in Reduced Network (Alphabetical Node Ordering):")
                    for x in at:
                        print(x)
            else:
                print("No Free Nodes Remain.")

    def summary(self,terminal_keys=None,show_original_rules=True,hide_rules=False,show_explicit_permutations=False):
        for motif_reduction in self.motif_reduction_list:
            if terminal_keys is None or motif_reduction.terminal in terminal_keys:
                print("__________________")
                motif_reduction.summary(show_original_rules=show_original_rules,hide_rules=hide_rules,show_explicit_permutations=show_explicit_permutations)

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
            target_indices = []
            for i,reduction in enumerate(self.motif_reduction_list):
                if logically_fixed.items() <= reduction.logically_fixed_nodes.items():
                    target_indices.append(i)
            return target_indices
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


    def reduction_drivers(self,target_index,method='internal',max_drivers=None,GRASP_iterations=None):
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
                path_motif_history += [x for x in self.motif_reduction_list[ind].motif_history if not x in path_motif_history]

                if method == 'internal':
                    history_drivers = sm_doi.internal_drivers(path_motif_history[-1],
                        self.motif_reduction_list[ind_prev].reduced_primes,
                        max_drivers=max_drivers)
                elif method == 'GRASP':
                    history_drivers = sm_doi.GRASP(path_motif_history[-1],
                        self.motif_reduction_list[ind_prev].reduced_primes,
                        GRASP_iterations = GRASP_iterations)
                    if len(history_drivers) == 0:
                        history_drivers = [path_motif_history[-1]]
                elif method == 'minimal':
                    history_drivers = sm_doi.minimal_drivers(path_motif_history[-1],
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
            elif method == 'minimal' or method == 'GRASP':
                drivers.append(path_drivers)
        return drivers

    def reprogram_to_trap_spaces(self,logically_fixed,target_method='history',driver_method='internal',max_drivers=None,GRASP_iterations=None,GRASP_score_override=None):
        """
        Find driver sets that lead to the node states specified by logically_fixed

        Inputs:
        logically_fixed - state dictionary specifying the control target
        max_drivers - the maximum number of driver nodes to attempt when looking
                      for stable motif driver nodes before specifying the entire
                      stable motif as part fo the driver set; not used in GRASP
                      methods.
        GRASP_iterations - number of times to construct GRASP driver sets; only
                           used in GRASP methods.

        GRASP_score_override - optional heuristic score function override (see
                                GRASP function for details). Only used in GRASP
                                methods.
        target_method - either 'history' or 'merge'; see Methods below for details
        driver_method - either 'internal', 'minimal', or 'GRASP' see Methods below
                        for details

        Output:
        nonredundant_drivers - control strategies found; interpretation depends
                               on method selected (see below)


        - Methods -
        history, internal:
        Finds all shortest stable motif histories that result in the target node states
        being logically fixed. Each stable motif is searched for internal driver nodes.
        The resulting internal drivers are combined into a single  control set. The
        return value consists of all such control sets for all  stable motif histories.
        Each control set eventually becomes self-sustaining.

        history, minimal:
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

        history, GRASP:
        The same as history, minimal, except external driver nodes are searched for
        using the GRASP algorithm using GRASP_iterations iterations.

        merge, internal:
        Finds all shortest stable motif histories that result in the target node states
        being logically fixed. All node states in the motifs in the history are merged
        into a stable module dictionary. This is then searched for internal driver
        nodes. Each element of the return value is a dictionary corresponding to a
        control set. Each control set eventually becomes self-sustaining.

        merge, minimal:
        Similar to the merge method, except the search for drivers is conducted over
        all nodes, not just those internal to the merged stable module. Furthermore,
        the search is truncated when a control set is found such that the search does
        not proceed to driver sets larger than the smallest found. Each element of
        the return value is a dictionary corresponding to a control set. The control
        sets are only guaranteed to result in activation of the target if they are
        temporary interventions.

        merge, GRASP:
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
                target_history = self.motif_reduction_list[target_index].motif_history
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

    def process_networkx_succession_diagram(self,include_attractors_in_diagram=True):
        """
        Generate network x graphs for the succession diagrams labeled with the stable motifs and the motif history.
        The succession diagram graphs are stored in the Succession object. We also store some position dictionaries
        for used for plotting these graphs with matplotlib. There are 4 succession diagram graphs:

        G_reduced_network_based, G_reduced_network_based_labeled: Each node denotes either the original network model or,
            a reduced network obtained after fixing the state of a stable motif. The edges denote a stable motif that
            when fixed, transforms one reduced network to another. G_reduced_network_based only has numbers as the
            node and edge names/labels, and is mostly used to plot the network with matplotlib. G_reduced_network_based_labeled
            has the reduced network motif history as the node labels, the stable motifs as the edge labels, and the matplotlib
            x,y plotting positions as node and edge attributes.
        G_motif_based, G_motif_based_labeled: Analogous to the above, but nodes denote stable motifs and edges denote reduced networks.
            These are obtained through a line_graph transformation of the G_reduced_network_based networks.

        There are 2 types of position dictionaries, pos_reduced_network_based and pos_motif_based, which store the positions for
        matplotlib plotting.

        Inputs:
        self - Succession object in which Succession.build_succession_diagram has been run
        include_attractors_in_diagram - Boolean variable indicating whether the attractors
                      will be included as sink nodes in the succession diagram graph

        Output:
        None, the Network x graphs and the matplotlib plot positions are formatted and stored in the Succession object 
        """

        self.G_reduced_network_based,self.G_reduced_network_based_labeled,self.pos_reduced_network_based,h_dict,h_dict_edges=self.networkx_succession_diagram_reduced_network_based()
        self.G_motif_based,self.G_motif_based_labeled,self.pos_motif_based=self.networkx_succession_diagram_motif_based(h_dict,h_dict_edges)
        if(include_attractors_in_diagram):
            self.add_attractors_networkx_diagram(h_dict,h_dict_edges)

    def networkx_succession_diagram_reduced_network_based(self):
        G_reduced_network_based=self.digraph.copy()
        self.motif_reduction_list_dictionary={i:set([frozenset(tuple(x.items())) for x in reduction.motif_history]) for i,reduction in enumerate(self.motif_reduction_list)}
        for u,v in G_reduced_network_based.edges():
                motif_set=self.motif_reduction_list_dictionary[v]-self.motif_reduction_list_dictionary[u]
                motif_list=[list(x) for x in list(motif_set)][0]
                G_reduced_network_based.edges[u, v]['label']=str(u)+"_"+str(v)+" ("+", ".join([x[0]+"="+str(x[1]) for x in motif_list])+")"
        h_dict=dict()
        for i,reduction in enumerate(self.motif_reduction_list):
                newlabel="" if len(reduction.motif_history)==0 else motif_history_text(reduction.motif_history)
                h_dict[i]=str(i)+" ["+newlabel+"]"
        G_reduced_network_based_labeled = nx.relabel_nodes(G_reduced_network_based.copy(), h_dict)
        pos_reduced_network_based = nx.spring_layout(G_reduced_network_based)#nx.nx_pydot.graphviz_layout(G_reduced_network_based, prog='dot')
        xmax=0
        ymax=0
        for node,(x,y) in pos_reduced_network_based.items():
                xmax=max(x,xmax)
                ymax=max(y,ymax)
        for node,(x,y) in pos_reduced_network_based.items():
                G_reduced_network_based_labeled.nodes[h_dict[node]]['x'] = 1.5*(xmax-float(x))
                G_reduced_network_based_labeled.nodes[h_dict[node]]['y'] = 1.5*(ymax-float(y))
                G_reduced_network_based_labeled.nodes[h_dict[node]]['label'] = str(h_dict[node])
        h_dict_edges={(u,v):G_reduced_network_based.edges[u, v]['label'] for u,v in G_reduced_network_based.edges()}
        return(G_reduced_network_based,G_reduced_network_based_labeled,pos_reduced_network_based,h_dict,h_dict_edges)

    def networkx_succession_diagram_motif_based(self,h_dict,h_dict_edges):
        G_motif_based = nx.line_graph(self.G_reduced_network_based)
        pos_motif_based = nx.spring_layout(G_motif_based)#nx.nx_pydot.graphviz_layout(G_motif_based, prog='dot')
        G_motif_based_labeled_temp=G_motif_based.copy()
        xmax=0
        ymax=0
        for node,(x,y) in pos_motif_based.items():
                xmax=max(x,xmax)
                ymax=max(y,ymax)
        for node,(x,y) in pos_motif_based.items():
                G_motif_based_labeled_temp.nodes[node]['x'] = 1.5*(float(x))
                G_motif_based_labeled_temp.nodes[node]['y'] = 1.5*(ymax-float(y))
                G_motif_based_labeled_temp.nodes[node]['label'] = str(h_dict_edges[node])
        for u,v in G_motif_based_labeled_temp.edges():
                G_motif_based_labeled_temp.edges[u,v]['label'] = h_dict[u[1]]
        G_motif_based_labeled=nx.relabel_nodes(G_motif_based_labeled_temp, h_dict_edges)

        return(G_motif_based,G_motif_based_labeled,pos_motif_based)

    def add_attractors_networkx_diagram(self,h_dict,h_dict_edges):
        G=self.G_reduced_network_based
        sink_nodes=[node for node, out_degree in G.out_degree if out_degree == 0]
        sink_node_attractor_dictionary={}
        max_node=len([node for node in G])
        nodes_attractors=list(range(max_node,max_node+len(self.attractor_fixed_nodes_list)))
        edges_attractors=[]
        node_labels_attractors=[(str(max_node+i)+" Attractor_"+str(i),{'label':str(max_node+i)+" Attractor_"+str(i)}) for i,a in enumerate(nodes_attractors)]
        for i,attractor in enumerate(self.attractor_fixed_nodes_list):
            for j,reduction in enumerate(self.motif_reduction_list):
                if j in sink_nodes:
                    if reduction.logically_fixed_nodes.items() <= attractor.items():
                        sink_node_attractor_dictionary[j]=max_node+i
                        edges_attractors.append((j,max_node+i,{'label':""}))
        h_dict.update({max_node+i:str(max_node+i)+" Attractor_"+str(i) for i,a in enumerate(nodes_attractors)})
        h_dict_edges.update({max_node+i:str(max_node+i)+" Attractor_"+str(i) for i,a in enumerate(nodes_attractors)})

        self.G_reduced_network_based.add_nodes_from(nodes_attractors)
        self.G_reduced_network_based.add_edges_from(edges_attractors)
        self.pos_reduced_network_based = nx.spring_layout(self.G_reduced_network_based)#nx.nx_pydot.graphviz_layout(self.G_reduced_network_based, prog='dot')
        self.G_reduced_network_based_labeled.add_nodes_from(node_labels_attractors)
        self.G_reduced_network_based_labeled.add_edges_from([(h_dict[e[0]],h_dict[e[1]],{'label':str(e[0])+"_"+str(e[1])}) for e in edges_attractors])
        G=self.G_reduced_network_based_labeled
        pos=self.pos_reduced_network_based
        xmax=0
        ymax=0
        for node,(x,y) in pos.items():
            xmax=max(x,xmax)
            ymax=max(y,ymax)
        for node,(x,y) in pos.items():
            G.nodes[h_dict[node]]['x'] = 1.5*(float(x))
            G.nodes[h_dict[node]]['y'] = 1.5*(ymax-float(y))

        G = self.G_motif_based
        sink_nodes_motif_based=[node for node, out_degree in G.out_degree if out_degree == 0]
        edges_attractors_motif_based=[(s,sink_node_attractor_dictionary[s[1]],{'label':""}) for s in sink_nodes_motif_based]
        self.G_motif_based.add_nodes_from(nodes_attractors)
        self.G_motif_based.add_edges_from(edges_attractors_motif_based)
        self.pos_motif_based = nx.spring_layout(self.G_motif_based)#nx.nx_pydot.graphviz_layout(self.G_motif_based, prog='dot')
        self.G_motif_based_labeled.add_nodes_from(node_labels_attractors)
        self.G_motif_based_labeled.add_edges_from([(h_dict_edges[e[0]],h_dict_edges[e[1]],{'label':h_dict[e[0][1]]}) for e in edges_attractors_motif_based])
        G=self.G_motif_based_labeled
        pos=self.pos_motif_based
        xmax=0
        ymax=0
        for node,(x,y) in pos.items():
            xmax=max(x,xmax)
            ymax=max(y,ymax)
        for node,(x,y) in pos.items():
            G.nodes[h_dict_edges[node]]['x'] = 1.5*(float(x))
            G.nodes[h_dict_edges[node]]['y'] = 1.5*(ymax-float(y))

    def plot_networkx_succession_diagram_motif_based(self,print_out_labels=False):
        edge_labels_motif_based={(u,v):str(u[1]) for u,v in self.G_motif_based.edges()}
        node_labels_motif_based={n:str(n[0])+"_"+str(n[1]) for n in self.G_motif_based.nodes() if type(n) is tuple}
        node_labels_motif_based.update({n:n for n in self.G_motif_based.nodes() if type(n) is int})
        nx.draw_networkx_nodes(self.G_motif_based, self.pos_motif_based, node_size=200,alpha=0.4)
        nx.draw_networkx_edges(self.G_motif_based, self.pos_motif_based, width=1, arrowsize=20, alpha=1, edge_color='r')
        nx.draw_networkx_edge_labels(self.G_motif_based, self.pos_motif_based, edge_labels = edge_labels_motif_based,font_size=10)
        nx.draw_networkx_labels(self.G_motif_based, self.pos_motif_based,labels = node_labels_motif_based, font_size=12, alpha=1)

        x_values, y_values = zip(*self.pos_motif_based.values())
        x_max = max(x_values)
        x_min = min(x_values)
        x_margin = (x_max - x_min) * 0.5
        plt.xlim(x_min - x_margin, x_max + x_margin)
        if(print_out_labels):
            print("Succession diagram nodes:\n")
            for node in self.G_motif_based_labeled.nodes():
                print(self.G_motif_based_labeled.nodes[node]['label'])
            print("\nSuccession diagram edges:\n")
            for label in sorted(list(set(nx.get_edge_attributes(self.G_motif_based_labeled,'label').values()))):
                print(label.replace("\n"," "))
        plt.show()

    def plot_networkx_succession_diagram_reduced_network_based(self,print_out_labels=False):
        edge_labels={(u,v):str(u)+"_"+str(v) for u,v in self.G_reduced_network_based.edges()}
        nx.draw_networkx_nodes(self.G_reduced_network_based, self.pos_reduced_network_based, node_size=200,alpha=0.4)
        nx.draw_networkx_edges(self.G_reduced_network_based, self.pos_reduced_network_based, width=1, arrowsize=20, alpha=1, edge_color='r')
        nx.draw_networkx_edge_labels(self.G_reduced_network_based, self.pos_reduced_network_based, edge_labels = edge_labels,font_size=10)
        nx.draw_networkx_labels(self.G_reduced_network_based, self.pos_reduced_network_based, font_size=12, alpha=1)

        x_values, y_values = zip(*self.pos_reduced_network_based.values())
        x_max = max(x_values)
        x_min = min(x_values)
        x_margin = (x_max - x_min) * 0.5
        plt.xlim(x_min - x_margin, x_max + x_margin)
        if(print_out_labels):
            print("Succession diagram nodes:\n")
            for node in self.G_reduced_network_based_labeled.nodes():
                print(node.replace("\n"," "))
            print("\nSuccession diagram edges:\n")
            for u,v in self.G_reduced_network_based_labeled.edges():
                print(self.G_reduced_network_based_labeled.edges[u, v]['label'])
        plt.show()

def build_succession_diagram(primes, fixed=None, motif_history=None, diagram=None, merge_equivalent_motifs=True,max_simulate_size=20,prioritize_source_motifs=True):
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
    myMotifReduction=sm_reduction.MotifReduction(motif_history,fixed.copy(),primes,max_simulate_size=max_simulate_size,prioritize_source_motifs=prioritize_source_motifs)
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
            np,fixed2 = sm_reduction.reduce_primes(sm,primes)
            fixed3 = fixed.copy()
            fixed3.update(fixed2)
            diagram = build_succession_diagram(np,fixed3,myMotifReduction.motif_history+[sm],
                diagram, merge_equivalent_motifs=merge_equivalent_motifs,
                max_simulate_size=max_simulate_size,
                prioritize_source_motifs=prioritize_source_motifs)
    return diagram

def motif_history_text(history):
    """
    Obtain the string version of the motif_history of a reduced network
    Given the motif_history of a reduction.from motif_reduction_list, obtain a text version of this motif_history

    Inputs:
    history - motif_history of a reduced network.from motif_reduction_list

    Outputs:
    String of motif history with each motif inside a parenthesis and separated by a line break
    """
    motif_history_str=""
    for motif in history:
        motif_str="("+", ".join([str(k)+"="+str(v) for k,v in motif.items()])+")"
        motif_history_str=motif_history_str+motif_str+"\n"

    return(motif_history_str[:-1])
