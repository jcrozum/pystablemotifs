import StableMotifs.Succession as sm_succession
import StableMotifs.Attractor as sm_attractor

class AttractorRepertoire:
    def __init__(self,primes,max_simulate_size=20,max_stable_motifs=10000):
        self.succession_diagram = None
        self.attractors = []
        self.reduction_attractors = {} # dict with values that are lists of matching attractors
        self.fewest_attractors = None
        self.most_attractors = None

        self.analyze_system(primes,max_simulate_size=max_simulate_size,max_stable_motifs=max_stable_motifs)


    def get_attractors_from_succession_diagram(self):
        for ri, reduction in enumerate(self.succession_diagram.motif_reduction_list):
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
                new_attractor = sm_attractor.Attractor(reduction,id)
                self.attractors.append(new_attractor)
                self.reduction_attractors[ri].append(new_attractor)

    def count_attractors(self):
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

    def analyze_system(self,primes,max_simulate_size=20,max_stable_motifs=10000):
        self.succession_diagram = sm_succession.build_succession_diagram(primes,max_simulate_size=max_simulate_size,max_stable_motifs=max_stable_motifs)
        self.get_attractors_from_succession_diagram()
        self.count_attractors()

    def summary(self):
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
