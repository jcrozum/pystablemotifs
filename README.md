# StableMotifs (WORK IN PROGRESS)
A set of tools for attractor and target control of Boolean systems. 
Includes stable motif reduction with oscillation checking for attractor identification and control, and Greedy Randomized Adaptive Search Procedure and brute-force methods for target control.

# Documentation
Extensive documentation is not yet available. See comments in the code for now if detailed usage instructions are needed.

# Requirements
PyBoolNet (v2.2.8) https://github.com/hklarner/PyBoolNet

Networkx (v2.3) https://github.com/networkx/networkx/

For getting PyBoolNet v2.2.8 working with Networkx v2.4, see https://github.com/hklarner/PyBoolNet/issues/33

# Features (Rough Outline)
Import networks in BooleanNet or BNet format.

Integration with PyBoolNet.

Find all attractors, including complex attractors, of a general asynchronous update Boolean system using the succession diagram method. For large networks, some complex attractors might not be identifiable in a reasonable time. In these cases, we place upper and lower bounds on the number of complex attractors.

Identify the stable motif attractor control strategies using one of several methods. Stable motifs can be "locked in" sequentially or collectively during the control search. Drivers for either individual stable motifs or all stable motifs taken together can be searched for using one of three methods: GRASP, brute-force for internal drivers, brute-force for any drivers (internal or external).

Succession diagram plotting.

Network reduction methods: Saadatpour et al. or Veliz-Cuba. Ability to remove the directed acyclic out component.

Generate Kauffman random boolean networks.

