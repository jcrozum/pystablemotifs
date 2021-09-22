# pystablemotifs
A set of tools for attractor and target control of Boolean systems.
Includes stable motif reduction with oscillation checking for attractor identification and control, and Greedy Randomized Adaptive Search Procedure and brute-force methods for target control.

The attractor identification algorithm is described in detail in

> J. C.  Rozum,  J. Gómez  Tejeda  Zañudo,  X. Gan,  D. Deritei,  R. Albert,  Parity  and  time reversal  elucidate  both  decision-making  in  empirical  models  and  attractor  scaling  in  critical  Boolean networks. Sci. Adv. 7 , eabf8124 (2021),

which is freely available here: https://advances.sciencemag.org/content/7/29/eabf8124.

A manuscript detailing the various control algorithms is in preparation. In the meantime, please see the documentation or contact the maintainers of this repository for details about these methods.

# Installation
Install with pip from GitHub (***recommended***):
`pip install git+https://github.com/jcrozum/pystablemotifs`

Install with pip from PyPI (***not recommended***, unless pyboolnet is already installed): `pip install pystablemotifs`

If you install from PyPI, you will need to install pyboolnet separately (instructions at https://github.com/hklarner/pyboolnet). This is because PyPI (apparently) does not support dependencies that are not also on PyPI.

# Documentation
See the basic usage example below, or the Tutorial.ipynb notebook for basic instructions. For advanced usage instructions, see Manual.pdf or contact the developers directly.

# Requirements
pyboolnet (v3.0.2+) https://github.com/hklarner/pyboolnet
(note: pyboolnet requires pyeda, which can be difficult to install in Windows;
    it is recommended to obtain a pyeda Windows wheel from https://www.lfd.uci.edu/~gohlke/pythonlibs/#pyeda)

Networkx (v2.4+) https://github.com/networkx/networkx/

Sympy (v1.5.1+) https://www.sympy.org/en/index.html

Pandas (v1.0.0+) https://pandas.pydata.org/

NumPy (v1.19.2+) https://numpy.org/

Matplotlib (v3.2.1+) https://matplotlib.org/

# Features
- Import networks in BooleanNet or BNet format

- Integration with pyboolnet

- Find and explore all attractors of a general asynchronous update Boolean system using the succession diagram method

- Place upper and lower bounds on the number of complex attractors in Boolean networks that are too large to fully analyze with available computational resources

- Identify attractor control strategies by leveraging stable motifs (maximal trap spaces arising from self-sustaining feedback loops)

- Search for drivers of key system behaviors using brute-force of Greedy Randomized Adaptive Search Procedure (GRASP) methods

- Plot succession diagrams, which highlight irreversible decision points in a stochastic system's trajectory

- Apply projection-based network reduction methods

- Generate Kauffman random boolean networks

# Basic usage example
In the example below, we import the Boolean model specified by the file test1.txt provided in the models folder. We then print its rules and finds its attractors, which are displayed in a condensed summary form.
```python
import pystablemotifs as sm

relative_path_to_model = "./models/simple_model.txt"
primes = sm.format.import_primes(relative_path_to_model)

print("RULES")
sm.format.pretty_print_prime_rules(primes)
print()

ar = sm.AttractorRepertoire.from_primes(primes)
ar.summary()
```
The output is as follows:
```
RULES
xA* = xB
xB* = xA
xC* = !xD | xA
xD* = xC
xE* = xB & xF
xF* = xE

There are 3 attractors.
{'xA': 0, 'xB': 0, 'xC': 'X', 'xD': 'X', 'xE': 0, 'xF': 0}

{'xA': 1, 'xB': 1, 'xC': 1, 'xD': 1, 'xE': 1, 'xF': 1}

{'xA': 1, 'xB': 1, 'xC': 1, 'xD': 1, 'xE': 0, 'xF': 0}
```
Alternatively, it is possible to import the Boolean rules from a string, as follows:
```python
rules="""xA* = !xA & !xB | xC
xB* = !xA & !xB | xC
xC* = xA & xB"""

primes = sm.format.create_primes(rules)
```
It is also possible to compute attractors and control interventions using default parameters and methods using the command line:
```cmd
python -m pystablemotifs "./models/simple_model.txt"
```
In this command, "./models/simple_model.txt" is the relative path to a model file containing Boolean rules.

For further examples, see the IPython notebooks in the "Examples and Tutorials" folder.
