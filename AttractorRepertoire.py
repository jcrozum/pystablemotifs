import StableMotifs.Succesion as sm_succession
import StableMotifs.Attractor as sm_attractor

class AttractorRepertoire:
    def __init__():
        self.attractors = []
        self.reduction_attractors = {} # dict with values that are lists of matching attractors

    def get_attractors_from_succession_diagram(succession_diagram):
        for ri, reduction in enumerate(succession_diagram.motif_reduction_list):
            if reduction.terminal == "yes": continue

            self.reduction_attractors[ri] = []

            duplicate = False # have we added the attractors for this reduciton yet?
            for attractor in self.attractors:
                if attractor.logically_fixed_nodes == reduction.logically_fixed_nodes:
                    duplicate = True
                    attractor.add_reduction(reduction)
                    self.reduction_attractors[ri].append(attractor)

            if duplicate: continue

            for id in reduction.attractor_ids:
                new_attractor = sm_attractor.Attractor(reduction,id)
                self.attractors.append(new_attractor)
                self.reduction_attractors[ri].append(new_attractor)

    def count_attractors():
        self.fewest_attractors = 0
        self.most_attractors = 0
        for attractor in self.attractors:
            if attractor.guaranteed:
                self.fewest_attractors += 1
            # else: self.fewest_attractors += 0

            if attractor.explored:
                self.most_attractors += 1
            else: # ludicrously conservative upper bound; assumes STG is all 2-cycles
                self.most_attractors = 2**(attractor.n_unfixed - 1)
