from pipeline import protein_protein_interaction_network

ppin = protein_protein_interaction_network.ProteinProteinInteractionNetwork()
ppin.add_proteins_from_excel("data/WT_vs_noninvasive_15min.xlsx", "U", 15)
ppin.add_proteins_from_excel("data/WT_vs_noninvasive_1h.xlsx", "U", 60)
ppin.add_proteins_from_excel("data/WT_vs_noninvasive_2h.xlsx", "U", 120)

#ppin.add_interactions_from_STRING()
ppin.add_interactions_from_BioGRID()
#ppin.add_interactions_from_IntAct()