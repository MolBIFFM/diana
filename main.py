import math

from ppi import protein_protein_interaction_network

ppin = protein_protein_interaction_network.ProteinProteinInteractionNetwork()
ppin.add_proteins_from_excel("data/WT_vs_noninvasive_15min.xlsx", "U", 15)
ppin.add_proteins_from_excel("data/WT_vs_noninvasive_1h.xlsx", "U", 60)
ppin.add_proteins_from_excel("data/WT_vs_noninvasive_2h.xlsx", "U", 120)

for t, col in [(15, "wt15minvscontrolratio"), (60, "wt60minvscontrolratio"),
               (120, "wt120minvscontrolratio")]:
    ppin.add_proteins_from_excel(
        file_name="data/mcp.M113.029918-2.xls",
        ptm="P",
        time=t,
        skiprows=1,
        protein_id_col="protein",
        protein_id_format=lambda x: x.split("|")[1],
        position_col="phosphosites",
        position_format=lambda x: x.split("_")[0],
        replicates=[col],
        convert_measurement=lambda x: math.log10(x) / math.log10(2.0))

ppin.add_interactions_from_STRING(combined_score=0.7)