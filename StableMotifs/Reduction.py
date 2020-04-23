import PyBoolNet
import itertools as it
import networkx as nx
import re

import StableMotifs.TimeReversal as sm_time
import StableMotifs.RestrictSpace as sm_rspace
import StableMotifs.Format as sm_format
import StableMotifs.Succession as sm_succession
import StableMotifs.DomainOfInfluence as sm_doi

def simplify_primes(primes):
    """
    Simplifies PyBoolNet primes (e.g., A | A & B becomes A)

    Input:
    primes - PyBoolNet primes describing the system update rules

    Output:
    a simplified version of primes
    """
    # reimport to force simplification
    if len(primes) > 0:
        return PyBoolNet.FileExchange.bnet2primes(PyBoolNet.FileExchange.primes2bnet(primes))
    else:
        return primes

def reduce_primes(fixed,primes):
    """
    Simplifies boolean rules when some nodes are held fixed

    Inputs:
    fixed - a dictionary of node states that are held fixed
    primes - PyBoolNet primes describing the system update rules

    Outputs:
    reduced_primes - PyBoolNet primes decribing the simplified update rules
    percolated_states - a dictionary of fixed node states (including the inputs) that were simplified and removed
    """
    reduced_primes = PyBoolNet.PrimeImplicants.create_constants(primes,fixed,Copy=True)
    percolated_states = PyBoolNet.PrimeImplicants._percolation(reduced_primes,True)
    percolated_states.update(fixed)


    return simplify_primes(reduced_primes), percolated_states

def delete_node(primes, node):
    """
    Reduces Boolean rules given by primes by deleting the variable specified by
    node. The deleted node may not appear in its own update function. Any update
    rules depending on the deleted node will have that dependence replaced by
    the update function of the deleted node. The rules are simplified after node
    deletion.

    Inputs:
    primes - PyBoolNet primes describing update rules
    node - the name of the node node to delete

    Outputs:
    new_primes - the reduced primes
    constants - a dictionary of nodes that became logically fixed during the
                simplification process
    """
    G = PyBoolNet.InteractionGraphs.primes2igraph(primes)

    assert not G.has_edge(node,node), ' '.join(["Node",str(node),"has a self-loop and cannot be deleted."])

    new_primes = {k:v for k,v in primes.items() if not k == node}
    constants = {}

    rule1 = sm_format.rule2bnet(primes[node][1])

    for child in G.successors(node):
        crule = sm_format.rule2bnet(primes[child][1])
        crule = re.sub(rf'\b{node}\b',"("+rule1+")",crule)
        crule = PyBoolNet.BooleanLogic.minimize_espresso(crule)
        crule = child + ",\t" + crule

        new_primes[child] = PyBoolNet.FileExchange.bnet2primes(crule)[child]
        nc = PyBoolNet.PrimeImplicants._percolation(new_primes,True)
        constants.update(nc)
    return new_primes, constants

def remove_outdag(primes):
    """
    Removes the terminal directed acyclic part of the regulatory network. This
    part of the network does not influence the attractor repertoire.
    """
    G = PyBoolNet.InteractionGraphs.primes2igraph(primes)
    od = PyBoolNet.InteractionGraphs.find_outdag(G)
    reduced = primes.copy()
    constants = {}
    for node in od:
        if node in reduced:
            reduced, nc = delete_node(reduced, node)
            constants.update(nc)
    return reduced, constants

def deletion_reduction(primes, max_in_degree = float('inf')):
    """
    Implements the reduction method of Veliz-Cuba (2011).
    Deletion order is such that nodes with low in-degree are prioritized for
    removal. Deletion proceeds until all remaining nodes have self-loops.
    """
    reduced, constants = remove_outdag(primes)
    G = PyBoolNet.InteractionGraphs.primes2igraph(reduced)
    cur_order = sorted(reduced,key=lambda x: G.in_degree(x))

    change = True
    while change and len(reduced) > 0:
        change = False
        for node in cur_order:
            retry_node = True
            if not node in reduced or G.in_degree(node) > max_in_degree:
                continue
            elif not any(node in p for p in reduced[node][1]):
                reduced, nc = delete_node(reduced, node)
                constants.update(nc)
                if len(reduced) > 0:
                    G = PyBoolNet.InteractionGraphs.primes2igraph(reduced)
                else:
                    G = nx.DiGraph()
                change = True
                break
        cur_order = sorted(reduced,key=lambda x: G.in_degree(x))

    return reduced, constants

def mediator_reduction(primes):
    """
    Network reduction method of Saadadtpour, Albert, Reluga (2013)
    Preserves number of fixed points and complex attractors, but may change
    qualitative features of complex attractors.
    """

    reduced, constants = remove_outdag(primes)
    cur_order = sorted(reduced)
    G = PyBoolNet.InteractionGraphs.primes2igraph(reduced)
    candidates = [v for v in reduced if G.in_degree(v) == G.out_degree(v) == 1 and not G.has_edge(v,v)]

    for node in candidates:
        u = list(G.predecessors(node))[0]
        w = list(G.successors(node))[0]
        if not w in G.successors(u) and not w in G.predecessors(u):
            reduced, nc = delete_node(reduced, node)
            constants.update(nc)
            G = PyBoolNet.InteractionGraphs.primes2igraph(reduced)
            candidates = [v for v in reduced if G.in_degree(v) == G.out_degree(v) == 1 and not G.has_edge(v,v)]

    return reduced, constants


class MotifReduction:
    """
    A reduced network with additional information stored;
    represents a node in a succession diagram (see Succession.py)

    Variables:
    motif_history - list of stable motifs that can lock in to give the reduced network (in order)
    merged_history_permutations - list of permutations of motif_history (by index) that are also valid
    logically_fixed_nodes - node state dictionary describing nodes that have been
                            fixed and reduced by stable motifs and their
                            logical domain of influence
    reduced_primes - update rules of the reduced network as PyBoolNet primes
    time_reverse_primes - update rules of the time reversed system as PyBoolNet primes
    stable_motifs - list of stable motifs in the reduced network
    time_reverse_stable_motifs - list of stable motifs in the time reversed system
    merged_source_motifs - stable motifs generated by merging the stable motifs corresponding to source nodes
    source_independent_motifs - stable motifs that exist independent of the values of the source nodes

    rspace - the rspace, or "restrict space" of the reduced network, describing a
             necessary condition for the system to avoid activating additional
             stable motifs (see RestrictSpace.py)
    fixed_rspace_nodes - nodes that are fixed in the rspace (stored as a state dictionary)
    rspace_constraint - a Boolean expression that is true in and only in the rspace
    reduced_rspace_constraint - a simplification of the rspace_constraint given
                                the fixed_rspace_nodes states are satisfied
    rspace_update_primes - the update rules obtained from simplifying under the
                           assumption that the fixed_rspace_nodes are fixed
    conserved_functions - a list of Boolean functions that are constant within
                          every attractor, in PyBoolNet update rule format
    rspace_attractor_candidates - attractors in the rspace_update_primes that
                                  satisfy the reduced_rspace_constraint

    partial_STG - subgraph of the state transition graph of the reduced network
                  that contains any and all attractors that do not lie in any
                  of the reduced network's stable motifs
    no_motif_attractors - list of complex attractors that do not "lock in" any additional stable motifs

    Functions:
    __init__(self,motif_history,fixed,reduced_primes,search_partial_STGs=True,prioritize_source_motifs=True)
    merge_source_motifs(self) - merges source node motifs into larger multi-node motifs for efficiency
    test_rspace(self) - for building rspace_attractor_candidates
    build_K0(self) - helper function for build_partial_STG
    build_inspace(self,ss,names) - helper function for build_partial_STG
    build_partial_STG(self) - for building partial_STG
    find_no_motif_attractors(self) - finds no_motif_attractors
    summary(self) - prints a summary of the MotifReduction to screen
    """
    def __init__(self,motif_history,fixed,reduced_primes,max_simulate_size=20,prioritize_source_motifs=True,max_stable_motifs=10000):
        if motif_history is None:
            self.motif_history = []
        else:
            self.motif_history = motif_history.copy()
        self.merged_history_permutations = []
        self.logically_fixed_nodes = fixed
        self.reduced_primes = reduced_primes.copy()
        self.time_reverse_primes =sm_time.time_reverse_primes(self.reduced_primes)
        self.stable_motifs = PyBoolNet.AspSolver.trap_spaces(self.reduced_primes, "max",MaxOutput=max_stable_motifs)
        self.time_reverse_stable_motifs = PyBoolNet.AspSolver.trap_spaces(self.time_reverse_primes, "max",MaxOutput=max_stable_motifs)

        self.merged_source_motifs=None
        self.source_independent_motifs=None
        if self.motif_history == [] and prioritize_source_motifs:
            self.merge_source_motifs()

        self.rspace = sm_rspace.rspace(self.stable_motifs, self.time_reverse_stable_motifs,self.reduced_primes)
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
                print("STG is too large to simulate ("+
                str(simulate_size)+"/"+str(max_simulate_size)+
                "). We will attempt reduction methods. Increase max_simulate_size to force simulation.")
                self.delprimes, self.attractor_constants = mediator_reduction(self.reduced_primes)
                self.delprimes, nc = deletion_reduction(self.delprimes)
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
                        print("The reduction indicates that the branch is not terminal. No need to simulate.")
                        return
                if len(self.delprimes) < max_simulate_size:
                    print("Simulating deletion reduction ("+str(len(self.delprimes))+" nodes)...")
                    self.find_deletion_no_motif_attractors()
                    if len(self.deletion_no_motif_attractors) == 0:
                        self.terminal = "no"
                    else:
                        self.terminal = "yes"
                else:
                    print("The STG is still too large ("+str(len(self.delprimes))+").")
                    print("Further analysis of this branch is needed.")



    def merge_source_motifs(self):
        """
        Merges stable motifs (and time-reversal stable motifs) that correspond to source nodes, e.g. A*=A, into combined motifs to
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

    def test_rspace(self, search_partial_STGs=True):
        STG=PyBoolNet.StateTransitionGraphs.primes2stg(self.rspace_update_primes,"asynchronous")
        steady_states,complex_attractors=PyBoolNet.Attractors.compute_attractors_tarjan(STG)
        names = sorted(self.rspace_update_primes)
        attractors = complex_attractors+[[s] for s in steady_states]
        self.build_rspace_attractor_candidates(attractors)

    def build_rspace_attractor_candidates(self,attractors):
        self.rspace_attractor_candidates = []
        for attractor in attractors:
            possible_rspace_attractor = True
            for state in attractor:
                # state_dict = {** sm_format.statestring2dict(state,names),**self.fixed_rspace_nodes}
                # if PyBoolNet.BooleanLogic.are_mutually_exclusive(self.rspace_constraint,
                #                                                  sm_format.implicant2bnet(state_dict)):
                if sm_rspace.partial_state_contradicts_rspace(sm_format.statestring2dict(state,names),self.rspace):
                    possible_rspace_attractor = False
                    break
            if possible_rspace_attractor:
                self.rspace_attractor_candidates.append(attractor)

        if len(self.rspace_attractor_candidates) == 0:
            self.terminal = "no"

    # Helper function for smart STG building
    def build_K0(self):
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

    # tests whether the (partial) state ss is in any stable motifs
    def in_motif(self,ss,names):
        for sm in self.stable_motifs:
            smin = True
            for i,r in enumerate(names):
                if r in sm and not int(ss[i]) == sm[r]:
                    smin = False
            if smin: return True
        return False

    # Helper function for smart STG building
    # List all tr stable_motifs to which (partial) state ss belongs
    def build_inspace(self,ss,names, tr_stable_motifs = None):
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

    def build_deletion_STG(self):
        names = sorted(self.delprimes)
        name_ind = {n:i for i,n in enumerate(names)}

        trprimes = sm_time.time_reverse_primes(self.delprimes)
        trsms = PyBoolNet.AspSolver.trap_spaces(trprimes,"max")


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

            self.deletion_STG.add_node(ss) # might end up removing later
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

    def find_deletion_no_motif_attractors(self):
        if self.deletion_STG is None:
            self.build_deletion_STG()

        # Note: fixed points of the deletion system are fixed points of the
        # undeleted system, so we ignore these as they must contain stable motifs
        candidates = [x for x in nx.attracting_components(self.deletion_STG) if len(x) > 1]

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

            self.partial_STG.add_node(ss) # might end up removing later
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
        if self.partial_STG is None:
            self.build_partial_STG()
        self.no_motif_attractors = list(nx.attracting_components(self.partial_STG))

    def summary(self,show_original_rules=True,hide_rules=False,show_explicit_permutations=False):
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
