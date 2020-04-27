import unittest
import StableMotifs as sm
import PyBoolNet

class test_StableMotifs(unittest.TestCase):

    rules='''A*=B
B*=A
C*=A or not D
D*=C
E*=B and F
F*=E
'''
    rules_pbn = sm.Format.booleannet2bnet(rules)
    primes = PyBoolNet.FileExchange.bnet2primes(rules_pbn)
    diag = sm.Succession.build_succession_diagram(primes)

    def test_nr_of_attractors(self):
        self.assertEqual(len(self.diag.attractor_fixed_nodes_list),3)

    def test_booleannet2bnet(self):
        self.assertEqual(self.rules_pbn, 'A,\tB\nB,\tA\nC,\tA | !D\nD,\tC\nE,\tB & F\nF,\tE\n')

    def test_attractor_fixed_nodes_list(self):
        self.assertDictEqual(self.diag.attractor_fixed_nodes_list[0], {'A': 0, 'B': 0, 'E': 0, 'F': 0})
        self.assertDictEqual(self.diag.attractor_fixed_nodes_list[1], {'A': 1, 'B': 1, 'C': 1, 'D': 1, 'E': 1, 'F': 1})
        self.assertDictEqual(self.diag.attractor_fixed_nodes_list[2], {'A': 1, 'B': 1, 'C': 1, 'D': 1, 'E': 0, 'F': 0})

    def test_motif_reduction_list(self):
        motif_history_list=[i.motif_history for i in self.diag.motif_reduction_list]
        validation_motif_history_list=[[],[{'A': 0, 'B': 0}],
        [{'A': 1, 'B': 1}],
        [{'A': 1, 'B': 1}, {'E': 1, 'F': 1}],
        [{'A': 1, 'B': 1}, {'E': 0, 'F': 0}],
        [{'E': 0, 'F': 0}],
        [{'E': 0, 'F': 0}, {'A': 0, 'B': 0}]]
        self.assertListEqual(motif_history_list,validation_motif_history_list)
        #pathological example with a complex attractor on a "ghost" branch

    rules_pathological='''xA*= not xA and not xB or xC
    xB*= not xA and not xB or xC
    xC*= xA and xB'''
    rules_pbn_pathological = sm.Format.booleannet2bnet(rules_pathological)
    primes_pathological = PyBoolNet.FileExchange.bnet2primes(rules_pbn_pathological)
    diag_pathological = sm.Succession.build_succession_diagram(primes_pathological)

    def test_ghost_branch(self):
        '''
        High level test of the correct identification of attractors that are not preceded by stable motif lockins.
        (ghost_branch is not a function in the module)
        '''

        self.assertListEqual(self.diag_pathological.reduced_complex_attractor_list,[[{'000', '010', '100'}], None])
        self.assertListEqual(self.diag_pathological.attractor_fixed_nodes_list,[{}, {'xA': 1, 'xB': 1, 'xC': 1}])

if __name__ == '__main__':
    unittest.main()
