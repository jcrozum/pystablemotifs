import PyBoolNet
import StableMotifs as sm
from timeit import default_timer
import pickle

#print("Loading EMT network . . .")
#primes = sm.Format.import_primes("models/test4.txt",remove_constants=True)
#reprogramming_target = {'E_cadherin':1}

rules = """PDGF*= 0
IL15*=1
Stimuli*=1
Stimuli2*=0
CD45*=0
TAX*=0
CTLA4*=TCR
TCR*= Stimuli and not CTLA4
PDGFR*=S1P or PDGF
FYN*= TCR or IL2RB
Cytoskeleton_signaling*= FYN
LCK*=CD45 or ((TCR or IL2RB) and not ZAP70)
ZAP70*= LCK and not FYN
GRB2*= IL2RB or ZAP70
PLCG1*=GRB2 or PDGFR
RAS*=(GRB2 or PLCG1) and not GAP
GAP*=(RAS or (PDGFR and GAP)) and not (IL15 or IL2)
MEK*=RAS
ERK*=MEK and PI3K
PI3K*=PDGFR or RAS
NFKB*=(TPL2 or PI3K) or (FLIP and TRADD and IAP)
NFAT*=PI3K
RANTES*=NFKB
IL2*=(NFKB or STAT3 or NFAT) and not TBET
IL2RBT*=ERK and TBET
IL2RB*=IL2RBT and (IL2 or IL15)
IL2RAT*=IL2 and (STAT3 or NFKB)
IL2RA*=(IL2 and IL2RAT) and not IL2RA
JAK*=(IL2RA or IL2RB or RANTES or IFNG) and not (SOCS or CD45)
SOCS*=JAK and not (IL2 or IL15)
STAT3*=JAK
P27*=STAT3
Proliferation*=STAT3 and not P27
TBET*=JAK or TBET
CREB*=ERK and IFNG
IFNGT*=TBET or STAT3 or NFAT
IFNG*=((IL2 or IL15 or Stimuli) and IFNGT) and not (SMAD or P2)
P2*=(IFNG) and not Stimuli2
GZMB*=(CREB and IFNG) or TBET
TPL2*=TAX or (PI3K and TNF)
TNF*=NFKB
TRADD*=TNF and not (IAP or A20)
FasL*=STAT3 or NFKB or NFAT or ERK
FasT*=NFKB
Fas*=(FasT and FasL) and not sFas
sFas*=FasT and S1P
Ceramide*=Fas and not S1P
DISC*=FasT and ((Fas and IL2) or Ceramide or (Fas and not FLIP))
Caspase*=(((TRADD or GZMB) and BID) and not IAP) or DISC
FLIP*=(NFKB or (CREB and IFNG)) and not DISC
A20*=NFKB
BID*=(Caspase or GZMB) and not (BclxL or MCL1)
IAP*=NFKB and not BID
BclxL*=(NFKB or STAT3) and not (BID or GZMB or DISC)
MCL1*=(IL2RB and STAT3 and NFKB and PI3K) and not DISC
Apoptosis*=Caspase
GPCR*=S1P
SMAD*=GPCR
SPHK1*=PDGFR
S1P*=SPHK1 and not Ceramide
"""
reprogramming_target = {'Apoptosis':1}
rules = sm.Format.booleannet2bnet(rules)
primes = PyBoolNet.FileExchange.bnet2primes(rules)
PyBoolNet.PrimeImplicants._percolation(primes,True)

# print("Brute-force search for knockout/knockins that achieve",reprogramming_target,". . .")
# start=default_timer()
# koki = sm.DomainOfInfluence.knock_to_partial_state(reprogramming_target,primes,max_drivers=2)
# end=default_timer()
# print("Time running brute-force search method:",end-start)
# print("Sets found:")
# for x in koki: print(x)


print("Building succession diagram . . .")
# We do not need the complex attractors for this example, and ruling them
# out in the EMT case is extremely slow
diag = sm.Succession.build_succession_diagram(primes)

print("Computing driver sets (in multiple ways) that reprogram to an attractor with ",reprogramming_target,". . .")
start=default_timer()
reprogram_sets_minimal = diag.reprogram_to_trap_spaces(reprogramming_target,method='minimal',max_drivers=None)
end=default_timer()
print()
print("Time running minimal method:",end-start)
print("Sets found:")
for x in reprogram_sets_minimal: print(x)

start=default_timer()
reprogram_sets_merge = diag.reprogram_to_trap_spaces(reprogramming_target,method='merge',max_drivers=None)
end=default_timer()
print()
print("Time running merge method:",end-start)
print("Sets found:")
for x in reprogram_sets_merge: print(x)

start=default_timer()
reprogram_sets_history = diag.reprogram_to_trap_spaces(reprogramming_target,method='history',max_drivers=None)
end=default_timer()
print()
print("Time running history method:",end-start)
print("Sets found:")
for x in reprogram_sets_history: print(x)
#
start=default_timer()
reprogram_sets_minimal_history = diag.reprogram_to_trap_spaces(reprogramming_target,method='minimal_history',max_drivers=None)
end=default_timer()
print()
print("Time running minimal_history method:",end-start)
print("Sets found:")
for x in reprogram_sets_minimal_history: print(x,"\n")
