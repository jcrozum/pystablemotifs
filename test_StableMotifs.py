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

if __name__ == '__main__':
    unittest.main()
