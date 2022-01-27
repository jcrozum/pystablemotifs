import pystablemotifs as sm
from pyboolnet.external.bnet2primes import bnet_text2primes

class Primes:
    """The class that stores information about pyboolnet formatted primes. Initialize using
    import_rules.

    Attributes
    ----------
    rules : modified pyboolnet primes dictionary
        The pbn rules.
    collapsed_rules : pyboolnet primes dictionary
        The pbn rules collapsed into GAU rules.

    """

    def __init__(self):
        self.probabilities = {}
        # {"node1": [prob1, ...], ...}
        self.rules = {}
        # {"node1": [rule1, ...] , ...}, rule1 in the form [[{}, ...],[{}, ...]]
        self.collapsed_rules = {}
        # {"Node1": rule, ...}
        self.collapsed_bnet = ""

    def import_rules(self, fname, format='BooleanNet', remove_constants=False):
        """Import boolean rules from file and return modified pyboolnet formatted primes list.

        Parameters
        ----------
        fname : str
            Path to (plaintext) file containing Boolean rules in format specified
            by the 'format' option.
            Path to Boolean Expressions folder in case of CellCollective format.
        format : str
             Boolean rule format; options are 'BooleanNet' or 'BNet' or 'CellCollective'
             (the default is 'BooleanNet').
        remove_constants : bool
             If True, variables that are constant are removed and their influence is
             percolated. Otherwise, they remain and we consider initial conditions
             in opposition to their values (the default is False).

        Returns
        -------
        modified pyboolnet primes dictionary
            Update rules in modified pyboolnet format that supports pbn.

        """
        primes = sm.format.import_primes(fname, format='BooleanNet', remove_constants=False)

        probabilities = {}
        rules = {}

        for node in primes:
            if "$" in node:
                x = node.split("$")
                # if the rule has a name
                if len(x) == 3:
                    # if this node is encountered for the first time
                    if x[0] not in rules:
                        probabilities[x[0]] = {}
                        probabilities[x[0]][x[2]] = int(x[1])
                        rules[x[0]] = {}
                        rules[x[0]][x[2]] = primes[node]
                    # if this node was encountered before
                    else:
                        probabilities[x[0]][x[2]] = int(x[1])
                        rules[x[0]][x[2]] = primes[node]
                # if the rule does not have a name
                else:
                    if x[0] not in rules:
                        probabilities[x[0]] = {}
                        probabilities[x[0]]["rule1"] = int(x[1])
                        rules[x[0]] = {}
                        rules[x[0]]["rule1"] = primes[node]
                    else:
                        rule_name = "rule" + str(len(probabilities[x[0]]) + 1)
                        probabilities[x[0]][rule_name] = int(x[1])
                        rules[x[0]][rule_name] = primes[node]
            else:
                continue

        for node in primes:
            if "$" not in node:
                if node in rules:
                    # primes cannot consider A, A$1, A$2 the same, and will store rules for A separately as a source node.
                    # need to avoid adding rules of A as a source node if there are pbn rules for A.
                    continue
                else:
                    probabilities[node] = {}
                    probabilities[node]["rule1"] = 1
                    rules[node] = {}
                    rules[node]["rule1"] = primes[node]
            else:
                continue

        self.probabilities = probabilities
        self.rules = rules
        self._get_collapsed_rules_from_rules(rules)

    def _get_collapsed_rules_from_rules(self, rules):
        """make GAU rules out of pbn rules
        """
        # not A and (rule1 or rule2 ...) or (rume1 and rule2 ...)
        collapsed_bnet = ""
        for node in rules:
            collapsed_bnet += node + ", !" + node + " & ( "
            for rule in rules[node]:
                collapsed_bnet += "(" + sm.format.rule2bnet(rules[node][rule][1]) + ")"
                collapsed_bnet += " | "
            collapsed_bnet = collapsed_bnet.removesuffix(" | ")
            collapsed_bnet += " ) | ( "
            for rule in rules[node]:
                collapsed_bnet += "(" + sm.format.rule2bnet(rules[node][rule][1]) + ")"
                collapsed_bnet += " & "
            collapsed_bnet = collapsed_bnet.removesuffix(" & ")
            collapsed_bnet += " )\n"

        collapsed_rules = bnet_text2primes(collapsed_bnet)

        self.collapsed_bnet = collapsed_bnet
        self.collapsed_rules = collapsed_rules
