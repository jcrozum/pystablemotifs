import PyBoolNet
import StableMotifs as sm
from timeit import default_timer
import pickle

print("Loading EMT network and creating succssion diagram . . .")
primes = sm.Format.import_primes("models/EMT.txt",remove_constants=True)
reprogramming_target_state = {'EMT':0}


# We do not need the complex attractors for this example, and ruling them
# out in the EMT case is extremely slow
diag = sm.Succession.build_succession_diagram(primes)
#diag.summary()
diag.attractor_candidate_summary()
# try:
#     diag = pickle.load( open( "EMT_SuccessionDiagram.p", "rb" ) )
# except:
#     diag = sm.Succession.build_succession_diagram(primes,search_partial_STGs=True)
#     pickle.dump( diag, open( "EMT_SuccessionDiagram.p", "wb" ) )

print("Computing driver sets (in multiple ways) that reprogram to an attractor with ",reprogramming_target_state," . . .")

# start=default_timer()
# reprogram_sets_minimal = diag.reprogram_to_trap_spaces(reprogramming_target_state,method='minimal')
# end=default_timer()
# print("Time running minimal method:",end-start)
# print("Sets found:")
# for x in reprogram_sets_minimal: print(x)

# start=default_timer()
# reprogram_sets_merge = diag.reprogram_to_trap_spaces(reprogramming_target_state,method='merge')
# end=default_timer()
# print("Time running merge method:",end-start)
# print("Sets found:")
# for x in reprogram_sets_merge: print(x)

# start=default_timer()
# reprogram_sets_history = diag.reprogram_to_trap_spaces(reprogramming_target_state,method='history')
# end=default_timer()
# print("Time running history method:",end-start)
# print("Sets found:")
# for x in reprogram_sets_history: print(x)
#
# start=default_timer()
# reprogram_sets_minimal_history = diag.reprogram_to_trap_spaces(reprogramming_target_state,method='minimal_history')
# end=default_timer()
# print("Time running minimal_history method:",end-start)
# print("Sets found:")
# for x in reprogram_sets_minimal_history: print(x)
