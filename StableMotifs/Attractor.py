import StableMotifs.Reduction as sm_reduction


## NOTE: this code assumes some changes to Reduction.py that are not implemented
# these include:
# 1) Attractor and their STGs are indexed by a new array called attractor_ids

class Attractor:

    def __init__(self,reduction, reduction_attractor_id):
        self.logically_fixed_nodes = reduction.logically_fixed_nodes
        self.representative = (reduction,reduction_attractor_id)
        self.reductions = [reduction]
        self.attractor_dict = reduction.attractor_dict_list[reduction_attractor_id]

        if reduction.no_motif_attractors is not None:
            self.stg = reduction.no_motif_attractors[reduction_attractor_id]
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

    def add_reduction(self,reduction):
        self.reductions.append(reduction)
