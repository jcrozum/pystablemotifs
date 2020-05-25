import StableMotifs.Reduction as sm_reduction


## NOTE: this code assumes some changes to Reduction.py that are not implemented
# these include:
# 1) Attractor and their STGs are indexed by a new array called attractor_ids

class Attractor:

    def __init__(reduction, reduction_attractor_id):
        self.logically_fixed_nodes = reduction.logically_fixed_nodes
        self.representative = (reduction,reduction_attractor_id)
        self.reductions = [reduction]
        self.attractor_dict = reduction.attractor_dicts[reduction_attractor_id]
        self.stg = reduction.attractor_stgs[reduction_attractor_id]

        self.fixed_nodes = {k:v for k,v in attractor_dict.items() if v in ['0','1']}
        if reduction.fixed_rspace_nodes is not None:
            self.fixed_nodes.update(reduction.fixed_rspace_nodes)

        self.oscillation_fixed_nodes = {k:v for k,v in fixed_nodes if k not in self.logically_fixed_nodes}
        self.n_unfixed = len(self.reduced_primes) - len(self.oscilation_fixed_nodes)
        if reduction.rspace_update_primes is not None:
            self.reduced_primes = reduction.rspace_update_primes
        else:
            self.reduced_primes = reduction.reduced_primes

        if self.stg is not None:
            self.size = len(self.stg.nodes())
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
                self.size_upper_bound = 2**(n_unfixed)

        if '?' in attractor_dict.values():
            self.explored = False
            self.guaranteed = True
        elif '!' in attractor_dict.values():
            self.explored = False
            self.guaranteed = False
        else: # X's, 1's, 0's
            self.explored = True
            self.guaranteed = True

        self.reductions = [reduction]

    def add_reduction(reduction):
        self.reductions.append(reduction)
