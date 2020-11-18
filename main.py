import math

from pipeline import protein_protein_interaction_network

ppin = protein_protein_interaction_network.ProteinProteinInteractionNetwork()
ppin.add_proteins_from_excel("data/WT_vs_noninvasive_15min.xlsx", "U", 15)
ppin.add_proteins_from_excel("data/WT_vs_noninvasive_1h.xlsx", "U", 60)
ppin.add_proteins_from_excel("data/WT_vs_noninvasive_2h.xlsx", "U", 120)

for time, column in [(15, "wt15minvscontrolratio"),
                     (60, "wt60minvscontrolratio"),
                     (120, "wt120minvscontrolratio")]:
    ppin.add_proteins_from_excel(
        "data/mcp.M113.029918-2.xls",
        "P",
        time,
        header=1,
        protein_id_col="protein",
        protein_id_format=lambda entry: entry.split("|")[1],
        position_col="phosphosites",
        position_format=lambda entry: entry.split("_")[0],
        replicates=[column],
        convert_measurement=lambda measurement: math.log10(measurement
                                                           ) / math.log10(2.0))

ppin.add_interactions_from_BioGRID()
ppin.add_interactions_from_IntAct()
ppin.add_interactions_from_STRING()

ppin.remove_isolates()

ppin.export_as_graphml("test.graphml")