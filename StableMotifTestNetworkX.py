import PyBoolNet
import StableMotifs as sm
from timeit import default_timer
import pickle
import networkx as nx
import pandas as pd

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

print("Loading modified TLGL network . . .")
rules = sm.Format.booleannet2bnet(rules)
primes = PyBoolNet.FileExchange.bnet2primes(rules)
PyBoolNet.PrimeImplicants._percolation(primes,True)

print("Building succession diagram . . .")
diag = sm.Succession.build_succession_diagram(primes)

print("Creating visualization for succession diagram network . . .")
G_reduced_network_based,G_reduced_network_based_labeled,pos_reduced_network_based,G_motif_based,G_motif_based_labeled,pos_motif_based=diag.process_networkx_succession_diagram(include_attractors_in_diagram=True)

print("Writing succession diagram network to file . . .")
nx.write_graphml(G_reduced_network_based_labeled, "SuccessionDiagram_reduced_network_based.graphml")
nx.write_gml(G_reduced_network_based_labeled, "SuccessionDiagram_reduced_network_based.gml")
nx.write_graphml(G_motif_based_labeled, "SuccessionDiagram_motif_based.graphml")
nx.write_gml(G_motif_based_labeled, "SuccessionDiagram_motif_based.gml")

print("Writing attractors to file . . .")
df_attractors=pd.DataFrame(columns=[node for node in diag.unreduced_primes])
df_attractors=df_attractors.merge(pd.DataFrame(diag.attractor_fixed_nodes_list),how="outer")
df_attractors.index=["Attractor_"+str(i) for i,value in enumerate(diag.attractor_fixed_nodes_list)]
df_attractors=df_attractors.fillna("X").astype(str).replace({"0.0": "0", "1.0": "1"}).sort_index(axis=1)
df_attractors.to_csv("Attractors.csv")

print("Visualization for succession diagram network using matplotlib")
sm.Succession.plot_networkx_succession_diagram_reduced_network_based(G_reduced_network_based,G_reduced_network_based_labeled,pos_reduced_network_based,print_out_labels=True)
sm.Succession.plot_networkx_succession_diagram_motif_based(G_motif_based,G_motif_based_labeled,pos_motif_based,print_out_labels=True)
