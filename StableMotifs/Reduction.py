import PyBoolNet
import itertools as it
import networkx as nx

from StableMotifs.TimeReversal import time_reverse_primes
from StableMotifs.RestrictSpace import rspace, fixed_rspace_nodes, reduce_rspace_string, attractor_space_candidates
from StableMotifs.Format import pretty_print_rspace, pretty_print_prime_rules, statestring2dict, implicant2bnet

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
    return reduced_primes, percolated_states

# Commented out because this function is no longer used
#def state_in_K(statestring,K,name_inds):
#    for kd in K:
#        in_kd = True
#        for k,v in kd.items():
#            if not int(statestring[name_inds[k]]) == v:
#                in_kd = False
#        if in_kd:
#            return True
#    return False


class MotifReduction:
    """
    A reduced network with additional information stored;
    represents a node in a succession diagram (see Succession.py)

    Variables:
    logically_fixed_nodes - node state dictionary describing nodes that have been
                            fixed and reduced by stable motifs and their
                            logical domain of influence
    reduced_primes - update rules of the reduced network as PyBoolNet primes
    time_reverse_primes - update rules of the time reversed system as PyBoolNet primes
    stable_motifs - list of stable motifs in the reduced network
    tr_stable_motifs - list of stable motifs in the time reversed system

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
    __init__(self,motif_history,fixed,reduced_primes)
    test_rspace(self) - for building rspace_attractor_candidates
    build_K0(self) - helper function for build_partial_STG
    build_inspace(self,ss,names) - helper function for build_partial_STG
    build_partial_STG(self) - for building partial_STG
    find_no_motif_attractors(self) - finds no_motif_attractors
    summary(self) - prints a summary of the MotifReduction to screen
    """
    def __init__(self,motif_history,fixed,reduced_primes,search_partial_STGs=True):
        if motif_history is None:
            self.motif_history = []
        else:
            self.motif_history = motif_history.copy()
        self.logically_fixed_nodes = fixed
        self.reduced_primes = reduced_primes.copy()

        self.time_reverse_primes = time_reverse_primes(self.reduced_primes)
        self.stable_motifs = PyBoolNet.AspSolver.trap_spaces(self.reduced_primes, "max")
        self.tr_stable_motifs = PyBoolNet.AspSolver.trap_spaces(self.time_reverse_primes, "max")
        self.rspace=rspace(self.stable_motifs,self.reduced_primes)

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

        study_possible_oscillation = False

        if self.rspace == [[{'0':1}]] and len(self.stable_motifs) > 0: # a stable motif must lock in
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
            self.fixed_rspace_nodes = fixed_rspace_nodes(self.rspace,self.reduced_primes)
            self.rspace_constraint = pretty_print_rspace(self.rspace)
            self.reduced_rspace_constraint = reduce_rspace_string(self.rspace_constraint,self.fixed_rspace_nodes)
            self.rspace_update_primes = reduce_primes(self.fixed_rspace_nodes,self.reduced_primes)[0]
            self.test_rspace()
            study_possible_oscillation = self.terminal == "possible" # value may be changed by test_rspace
        if study_possible_oscillation:
            self.conserved_functions = attractor_space_candidates(self.stable_motifs,
                                                                  self.tr_stable_motifs)
            if search_partial_STGs:
                self.find_no_motif_attractors()
                if len(self.no_motif_attractors) == 0:
                    self.terminal = "no"
                else:
                    self.terminal = "yes"

    def test_rspace(self):
        STG=PyBoolNet.StateTransitionGraphs.primes2stg(self.rspace_update_primes,"asynchronous")
        steady_states,complex_attractors=PyBoolNet.Attractors.compute_attractors_tarjan(STG)
        names = sorted(self.rspace_update_primes)
        attractors = complex_attractors+[[s] for s in steady_states]
        self.rspace_attractor_candidates = []

        for attractor in attractors:
            possible_rspace_attractor = True
            for state in attractor:
                state_dict = {**statestring2dict(state,names),**self.fixed_rspace_nodes}
                if PyBoolNet.BooleanLogic.are_mutually_exclusive(self.rspace_constraint,
                                                                 implicant2bnet(state_dict)):
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

    # Helper function for smart STG building
    # List all tr stable_motifs to which state ss belongs
    def build_inspace(self,ss,names):
        inspaces = []
        for ts in self.tr_stable_motifs:
            tsin = True
            for i,r in enumerate(names):
                if r in ts and not int(ss[i]) == ts[r]:
                    tsin = False
            if tsin: inspaces.append(ts)
        return inspaces

    def build_partial_STG(self):
        names = sorted(self.reduced_primes)
        name_ind = {n:i for i,n in enumerate(names)}
        N = len(names)
        K = self.build_K0()

        self.partial_STG = nx.DiGraph()

        inspace_dict = {}

        for s in it.product(['0','1'],repeat=N):
            ss = ''.join(s)
            if ss in K: continue
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

                        # Check if landed in K
                        # If not, check if we left a TR stable motif
                        prune = child_state in K
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

    def summary(self):
        print("Motif History:",self.motif_history)
        print()
        print("Logically Fixed Nodes:",self.logically_fixed_nodes)
        print()
        print("Reduced Update Rules:")
        pretty_print_prime_rules(self.reduced_primes)
        print()
        if self.terminal == "no":
            print("At least one additional stable motif must stabilize.")
            print()
            print("Stable motifs:", self.stable_motifs)
        elif self.terminal == "yes":
            if len(self.reduced_primes) > 0:
                print("At least some of the following must oscillate:")
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
                pretty_print_prime_rules(self.rspace_update_primes)

        if not self.conserved_functions is None:
            print()
            if len(self.conserved_functions) > 0:
                print("Found the following functions that are constant on attractors in this branch:")
                for x in self.conserved_functions:
                    if len(x) > 0:
                        pretty_print_rspace([x],silent=False)
                        print()
            else:
                print("Unable to find non-trivial conserved functions for attractors in this branch.")
                print()
            if not self.no_motif_attractors is None:
                if len(self.no_motif_attractors) > 0:
                    print("Found the following complex attractors that do not lock in additional stable motifs:")
                    for x in self.no_motif_attractors:
                        print(x)
