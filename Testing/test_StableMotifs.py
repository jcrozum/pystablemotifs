import sys
sys.path.append('../')

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

    def test_motif_reduction_dict(self):
        motif_history_list=[i.motif_history for i in self.diag.motif_reduction_dict.values()]
        validation_motif_history_list=[[],[{'A': 0, 'B': 0}],
        [{'A': 1, 'B': 1}],
        [{'A': 1, 'B': 1}, {'E': 1, 'F': 1}],
        [{'A': 1, 'B': 1}, {'E': 0, 'F': 0}],
        [{'E': 0, 'F': 0}],
        [{'E': 0, 'F': 0}, {'A': 0, 'B': 0}]]
        self.assertListEqual(motif_history_list,validation_motif_history_list)

    def test_AttractorRepertoire_attractor_states(self):
        max_simulate_size=20
        ar = sm.AttractorRepertoire.from_primes(self.primes, max_simulate_size=max_simulate_size)
        attractors_dict_list=[]
        for a in ar.attractors:
            attractors_dict_list.append(a.attractor_dict)
        self.assertListEqual(attractors_dict_list, [{'A': 0, 'B': 0, 'C': 'X', 'D': 'X', 'E': 0, 'F': 0},
                                                    {'A': 1, 'B': 1, 'C': 1, 'D': 1, 'E': 1, 'F': 1},
                                                    {'A': 1, 'B': 1, 'C': 1, 'D': 1, 'E': 0, 'F': 0}])


    #"pathological" example with a complex attractor on a "ghost" branch
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

    def test_ghost_branch_with_max_simulate_size_0(self):
        '''
        If we don't simulate the STG with brute force at all (max_simulate_size=0) the output should give a warning
        on the form of a "!" for nodes that potentially oscillate.
        '''
        max_simulate_size=0
        ar = sm.AttractorRepertoire.from_primes(self.primes_pathological, max_simulate_size=max_simulate_size)
        attractors_dict_list=[]
        for a in ar.attractors:
            attractors_dict_list.append(a.attractor_dict)
        self.assertListEqual(attractors_dict_list, [{'xA': '!', 'xB': '!', 'xC': 0}, {'xA': 1, 'xB': 1, 'xC': 1}])

    def test_two_complex_attractors_on_the_same_branch(self):

        rules='''a* = !b | Z
                b* = c | Z
                c* = a & ~Z
                A* = ~B | Y
                B* = C | Y
                C* = A & ~Y
                Z* = ~(A & B & ~C)
                Y* = ~(a & b & ~c)
                '''
        rules_pbn = sm.Format.booleannet2bnet(rules)
        primes = PyBoolNet.FileExchange.bnet2primes(rules_pbn)

        max_simulate_size=20
        ar = sm.AttractorRepertoire.from_primes(primes, max_simulate_size=max_simulate_size)
        attractors_dict_list=[]
        for a in ar.attractors:
            attractors_dict_list.append(a.attractor_dict)
        self.assertListEqual(attractors_dict_list, [{'A': '1', 'B': '1', 'C': '0', 'Y': '1', 'Z': '0', 'a': 'X', 'b': 'X', 'c': 'X'},
                                                    {'A': 'X', 'B': 'X', 'C': 'X', 'Y': '0', 'Z': '1', 'a': '1', 'b': '1', 'c': '0'}])
        max_simulate_size=0
        ar = sm.AttractorRepertoire.from_primes(primes, max_simulate_size=max_simulate_size)
        attractors_dict_list=[]
        for a in ar.attractors:
            attractors_dict_list.append(a.attractor_dict)
        self.assertListEqual(attractors_dict_list, [{'A': '?', 'B': '?', 'C': '?', 'Y': '?', 'Z': '?', 'a': '?', 'b': '?', 'c': '?'}])

#Testing functions of Export.py


    def test_networkx_succession_diagram_reduced_network_based(self):
        import StableMotifs.Export as ex
        rules='''A*=B
                B*=A
                C*=A or not D
                D*=C
                E*=B and F
                F*=E'''

        rules_pbn = sm.Format.booleannet2bnet(rules)
        primes = PyBoolNet.FileExchange.bnet2primes(rules_pbn)
        max_simulate_size=20
        include_attractors_in_diagram=False
        ar = sm.AttractorRepertoire.from_primes(primes, max_simulate_size=max_simulate_size)
        GR=ex.networkx_succession_diagram_reduced_network_based(ar,include_attractors_in_diagram=include_attractors_in_diagram)
        self.assertDictEqual(dict(GR.nodes(data=True)),{0: {'label': '', 'virtual_nodes': []},
                                                         1: {'label': '{A: 0, B: 0}', 'virtual_nodes': [{'A': 0, 'B': 0}]},
                                                         2: {'label': '{A: 1, B: 1}', 'virtual_nodes': [{'A': 1, 'B': 1}]},
                                                         3: {'label': '{A: 1, B: 1}, {E: 1, F: 1}',
                                                          'virtual_nodes': [{'A': 1, 'B': 1}, {'E': 1, 'F': 1}]},
                                                         4: {'label': '{A: 1, B: 1}, {E: 0, F: 0}',
                                                          'virtual_nodes': [{'A': 1, 'B': 1}, {'E': 0, 'F': 0}]},
                                                         5: {'label': '{E: 0, F: 0}', 'virtual_nodes': [{'E': 0, 'F': 0}]},
                                                         6: {'label': '{E: 0, F: 0}, {A: 0, B: 0}',
                                                          'virtual_nodes': [{'E': 0, 'F': 0}, {'A': 0, 'B': 0}]}})
        self.assertListEqual(list(GR.edges()), [(0, 1), (0, 2), (0, 5), (2, 3), (2, 4), (5, 6)])

    def test_networkx_succession_diagram_reduced_network_based(self):
        import StableMotifs.Export as ex
        rules='''A*=B
                B*=A
                C*=A or not D
                D*=C
                E*=B and F
                F*=E'''

        rules_pbn = sm.Format.booleannet2bnet(rules)
        primes = PyBoolNet.FileExchange.bnet2primes(rules_pbn)
        max_simulate_size=20
        include_attractors_in_diagram=False
        ar = sm.AttractorRepertoire.from_primes(primes, max_simulate_size=max_simulate_size)
        GR=ex.networkx_succession_diagram_reduced_network_based(ar,include_attractors_in_diagram=include_attractors_in_diagram)
        self.assertDictEqual(dict(GR.nodes(data=True)),{0: {'label': '', 'virtual_nodes': []},
                                                         1: {'label': '{A: 0, B: 0}', 'virtual_nodes': [{'A': 0, 'B': 0}]},
                                                         2: {'label': '{A: 1, B: 1}', 'virtual_nodes': [{'A': 1, 'B': 1}]},
                                                         3: {'label': '{A: 1, B: 1}, {E: 1, F: 1}',
                                                          'virtual_nodes': [{'A': 1, 'B': 1}, {'E': 1, 'F': 1}]},
                                                         4: {'label': '{A: 1, B: 1}, {E: 0, F: 0}',
                                                          'virtual_nodes': [{'A': 1, 'B': 1}, {'E': 0, 'F': 0}]},
                                                         5: {'label': '{E: 0, F: 0}', 'virtual_nodes': [{'E': 0, 'F': 0}]},
                                                         6: {'label': '{E: 0, F: 0}, {A: 0, B: 0}',
                                                          'virtual_nodes': [{'E': 0, 'F': 0}, {'A': 0, 'B': 0}]}})
        self.assertListEqual(list(GR.edges()), [(0, 1), (0, 2), (0, 5), (2, 3), (2, 4),(5,4),(5, 6)])

        include_attractors_in_diagram=True
        GR=ex.networkx_succession_diagram_reduced_network_based(ar,include_attractors_in_diagram=include_attractors_in_diagram)
        self.assertDictEqual(dict(GR.nodes(data=True)),{0: {'label': '', 'virtual_nodes': []},
                                                         1: {'label': '{A: 0, B: 0}', 'virtual_nodes': [{'A': 0, 'B': 0}]},
                                                         2: {'label': '{A: 1, B: 1}', 'virtual_nodes': [{'A': 1, 'B': 1}]},
                                                         3: {'label': '{A: 1, B: 1}, {E: 1, F: 1}',
                                                          'virtual_nodes': [{'A': 1, 'B': 1}, {'E': 1, 'F': 1}]},
                                                         4: {'label': '{A: 1, B: 1}, {E: 0, F: 0}',
                                                          'virtual_nodes': [{'A': 1, 'B': 1}, {'E': 0, 'F': 0}]},
                                                         5: {'label': '{E: 0, F: 0}', 'virtual_nodes': [{'E': 0, 'F': 0}]},
                                                         6: {'label': '{E: 0, F: 0}, {A: 0, B: 0}',
                                                          'virtual_nodes': [{'E': 0, 'F': 0}, {'A': 0, 'B': 0}]},
                                                         'A0': {'label': '{A: 0, B: 0, C: X, D: X, E: 0, F: 0}',
                                                          'virtual_nodes': {'A': 0, 'B': 0, 'C': 'X', 'D': 'X', 'E': 0, 'F': 0}},
                                                         'A1': {'label': '{A: 1, B: 1, C: 1, D: 1, E: 1, F: 1}',
                                                          'virtual_nodes': {'A': 1, 'B': 1, 'C': 1, 'D': 1, 'E': 1, 'F': 1}},
                                                         'A2': {'label': '{A: 1, B: 1, C: 1, D: 1, E: 0, F: 0}',
                                                          'virtual_nodes': {'A': 1, 'B': 1, 'C': 1, 'D': 1, 'E': 0, 'F': 0}}})
        self.assertListEqual(list(GR.edges()), [(0, 1),
                                                 (0, 2),
                                                 (0, 5),
                                                 (1, 'A0'),
                                                 (2, 3),
                                                 (2, 4),
                                                 (3, 'A1'),
                                                 (4, 'A2'),
                                                 (5,4),
                                                 (5, 6),
                                                 (6, 'A0')])

    def test_networkx_succession_diagram_motif_based(self):
        import StableMotifs.Export as ex
        rules='''A*=B
                B*=A
                C*=A or not D
                D*=C
                E*=B and F
                F*=E'''

        rules_pbn = sm.Format.booleannet2bnet(rules)
        primes = PyBoolNet.FileExchange.bnet2primes(rules_pbn)
        max_simulate_size=20
        include_attractors_in_diagram=False
        ar = sm.AttractorRepertoire.from_primes(primes, max_simulate_size=max_simulate_size)
        GM=ex.networkx_succession_diagram_motif_based(ar,include_attractors_in_diagram=include_attractors_in_diagram)
        self.assertDictEqual(dict(GM.nodes(data=True)),{(0, 1): {'label': '{A: 0, B: 0}', 'virtual_nodes': {'A': 0, 'B': 0}},
                                                         (0, 2): {'label': '{B: 1, A: 1}', 'virtual_nodes': {'B': 1, 'A': 1}},
                                                         (2, 3): {'label': '{E: 1, F: 1}', 'virtual_nodes': {'E': 1, 'F': 1}},
                                                         (2, 4): {'label': '{E: 0, F: 0}', 'virtual_nodes': {'E': 0, 'F': 0}},
                                                         (0, 5): {'label': '{E: 0, F: 0}', 'virtual_nodes': {'E': 0, 'F': 0}},
                                                         (5, 4): {'label': '{B: 1, A: 1}', 'virtual_nodes': {'B': 1, 'A': 1}},
                                                         (5, 6): {'label': '{A: 0, B: 0}', 'virtual_nodes': {'A': 0, 'B': 0}}})
        self.assertListEqual(list(GM.edges()), [((0, 2), (2, 3)), ((0, 2), (2, 4)), ((0, 5), (5, 4)), ((0, 5), (5, 6))])
        include_attractors_in_diagram=True
        ar = sm.AttractorRepertoire.from_primes(primes, max_simulate_size=max_simulate_size)
        GM=ex.networkx_succession_diagram_motif_based(ar,include_attractors_in_diagram=include_attractors_in_diagram)
        self.assertDictEqual(dict(GM.nodes(data=True)),{(0, 1): {'label': '{A: 0, B: 0}', 'virtual_nodes': {'A': 0, 'B': 0}},
                                                         (0, 2): {'label': '{B: 1, A: 1}', 'virtual_nodes': {'B': 1, 'A': 1}},
                                                         (2, 3): {'label': '{E: 1, F: 1}', 'virtual_nodes': {'E': 1, 'F': 1}},
                                                         (2, 4): {'label': '{E: 0, F: 0}', 'virtual_nodes': {'E': 0, 'F': 0}},
                                                         (0, 5): {'label': '{E: 0, F: 0}', 'virtual_nodes': {'E': 0, 'F': 0}},
                                                         (5, 4): {'label': '{B: 1, A: 1}', 'virtual_nodes': {'B': 1, 'A': 1}},
                                                         (5, 6): {'label': '{A: 0, B: 0}', 'virtual_nodes': {'A': 0, 'B': 0}},
                                                         'A0': {'label': '{A: 0, B: 0, C: X, D: X, E: 0, F: 0}',
                                                          'virtual_nodes': {'A': 0, 'B': 0, 'C': 'X', 'D': 'X', 'E': 0, 'F': 0}},
                                                         'A1': {'label': '{A: 1, B: 1, C: 1, D: 1, E: 1, F: 1}',
                                                          'virtual_nodes': {'A': 1, 'B': 1, 'C': 1, 'D': 1, 'E': 1, 'F': 1}},
                                                         'A2': {'label': '{A: 1, B: 1, C: 1, D: 1, E: 0, F: 0}',
                                                          'virtual_nodes': {'A': 1, 'B': 1, 'C': 1, 'D': 1, 'E': 0, 'F': 0}}})
        self.assertListEqual(list(GM.edges()), [((0, 1), 'A0'),
                                             ((0, 2), (2, 3)),
                                             ((0, 2), (2, 4)),
                                             ((2, 3), 'A1'),
                                             ((2, 4), 'A2'),
                                             ((0, 5), (5, 4)),
                                             ((0, 5), (5, 6)),
                                             ((5, 4), 'A2'),
                                             ((5, 6), 'A0')])

        def test_networkx_succession_diagram_reduced_network_based_pathological_example(self):
            import StableMotifs.Export as ex
            rules='''
            xA*= not xA and not xB or xC
            xB*= not xA and not xB or xC
            xC*= xA and xB
            '''
            rules_pbn = sm.Format.booleannet2bnet(rules)
            primes = PyBoolNet.FileExchange.bnet2primes(rules_pbn)
            max_simulate_size=20
            include_attractors_in_diagram=False
            ar = sm.AttractorRepertoire.from_primes(primes, max_simulate_size=max_simulate_size)
            GR=ex.networkx_succession_diagram_reduced_network_based(ar,include_attractors_in_diagram=include_attractors_in_diagram)
            self.assertDictEqual(dict(GR.nodes(data=True)),{0: {'label': '', 'virtual_nodes': []},
                                                         1: {'label': '{xA: 1, xB: 1, xC: 1}',
                                                          'virtual_nodes': [{'xA': 1, 'xB': 1, 'xC': 1}]}})
            self.assertListEqual(list(GR.edges()), [(0, 1)])

            include_attractors_in_diagram=True
            GR=ex.networkx_succession_diagram_reduced_network_based(ar,include_attractors_in_diagram=include_attractors_in_diagram)
            self.assertDictEqual(dict(GR.nodes(data=True)),{0: {'label': '', 'virtual_nodes': []},
                                                         1: {'label': '{xA: 1, xB: 1, xC: 1}',
                                                          'virtual_nodes': [{'xA': 1, 'xB': 1, 'xC': 1}]},
                                                         'A0': {'label': '{xA: X, xB: X, xC: 0}',
                                                          'virtual_nodes': {'xA': 'X', 'xB': 'X', 'xC': 0}},
                                                         'A1': {'label': '{xA: 1, xB: 1, xC: 1}',
                                                          'virtual_nodes': {'xA': 1, 'xB': 1, 'xC': 1}}})
            self.assertListEqual(list(GR.edges()), [(0, 1), (0, 'A0'), (1, 'A1')])

if __name__ == '__main__':
    unittest.main()
