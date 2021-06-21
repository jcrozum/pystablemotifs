# PyStableMotifs
A set of tools for attractor and target control of Boolean systems.
Includes stable motif reduction with oscillation checking for attractor identification and control, and Greedy Randomized Adaptive Search Procedure and brute-force methods for target control.

# Installation
Install with pip:
`pip install git+https://github.com/jcrozum/StableMotifs`

# Documentation
See the basic usage example below, or the Tutorial.ipynb notebook for basic instructions. For advanced usage instructions, see Manual.pdf or contact the developers directly.

# Requirements
PyBoolNet (v2.3.0) https://github.com/hklarner/PyBoolNet

Networkx (v2.4+) https://github.com/networkx/networkx/

Sympy (v1.5.1+) https://www.sympy.org/en/index.html

Pandas (v1.0.0+) https://pandas.pydata.org/

NumPy (v1.19.2+) https://numpy.org/

Matplotlib (v3.2.1+) https://matplotlib.org/

# Features
- Import networks in BooleanNet or BNet format

- Integration with PyBoolNet

- Find and explore all attractors of a general asynchronous update Boolean system using the succession diagram method

- Place upper and lower bounds on the number of complex attractors in Boolean networks that are too large to fully analyze with available computational resources

- Identify attractor control strategies by leveraging stable motifs (maximal trap spaces arising from self-sustaining feedback loops)

- Search for drivers of key system behaviors using brute-force of Greedy Randomized Adaptive Search Procedure (GRASP) methods

- Plot succession diagrams, which highlight irreversible decision points in a stochastic system's trajectory

- Apply projection-based network reduction methods

- Generate Kauffman random boolean networks

# Basic usage example
In the example below, we import the Boolean model specified by the file test1.txt provided in the models folder. We then print its rules and finds its attractors, which are displayed in a condensed summary form.

    import PyStableMotifs as sm

    relative_path_to_model = "./models/test1.txt"
    primes = sm.Format.import_primes(relative_path_to_model)

    print("RULES")
    sm.Format.pretty_print_prime_rules({k:primes[k] for k in sorted(primes)})
    print()

    ar = sm.AttractorRepertoire.from_primes(primes)
    ar.summary()

The output is as follows:

    RULES
    xA* = !xA & !xB | xC
    xB* = !xA & !xB | xC
    xC* = xA & xB

    There are 2 attractors.
    {'xA': 'X', 'xB': 'X', 'xC': 0}
    {'xA': 1, 'xB': 1, 'xC': 1}

Alternatively, it is possible to import the Boolean rules from a string, as follows:

    rules="""xA* = !xA & !xB | xC
        xB* = !xA & !xB | xC
        xC* = xA & xB"""
    primes = sm.Format.create_primes(rules)

For further examples, see the IPython notebook Tutorial.ipynb.
