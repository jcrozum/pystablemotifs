def main():
    import pystablemotifs as sm
    from sys import argv
    """
    Example Usage:
    python -m pystablemotifs <model file>

    here <model file> is the path to a model file, e.g. "./models/simple_model.txt".

    This will print the rules, a summary of the system's attractor repertoire, and
    trap space control interventions using default methods and parameters.
    """

    print("""
    This function is for demonstration. It runs the pystablemotifs attractor
    identification routine and identifies control interventions to target each
    identified attractor's trap space, provided it is nontrivial. Default function
    paramters are used. See files in the "Examples and Tutorials" directory for
    instructions regarding typical use cases.
    """)

    primes = sm.format.import_primes(argv[1],remove_constants=True)
    print("RULES")
    sm.format.pretty_print_prime_rules(primes)
    print()
    print("Analyzing network . . .")
    ar = sm.AttractorRepertoire.from_primes(primes,max_simulate_size=20)
    print("Analysis complete.")
    print()
    ar.summary()
    print("-"*10)

    print("\nComputing attractor control . . .")
    for attractor in ar.attractors:
        target = attractor.logically_fixed_nodes
        if len( target ) > 0 and attractor.guaranteed == True:
            print("-"*10)
            print("Target:",target)
            intervention = ar.reprogram_to_trap_spaces(target)
            print("Interventions:")
            for x in intervention: print(x)
    print("Complete.")

if __name__=='__main__':
    main()
