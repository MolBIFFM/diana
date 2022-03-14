from databases import uniprot
from download import download


def get_proteins(neighborhood=0.0,
                 neighborhood_transferred=0.0,
                 fusion=0.0,
                 cooccurence=0.0,
                 homology=0.0,
                 coexpression=0.0,
                 coexpression_transferred=0.0,
                 experiments=0.0,
                 experiments_transferred=0.0,
                 database=0.0,
                 database_transferred=0.0,
                 textmining=0.0,
                 textmining_transferred=0.0,
                 combined_score=0.0,
                 physical=False,
                 taxon_identifier=9606,
                 version=11.5):
    uniprot_id_map = {}
    for row in download.tabular_txt(
            "https://stringdb-static.org/download/protein.aliases.v{version}/{taxon_identifier}.protein.aliases.v{version}.txt.gz"
            .format(taxon_identifier=taxon_identifier, version=version),
            delimiter="\t",
            skiprows=1,
            usecols=[0, 1, 2],
    ):
        if "UniProt_AC" in row[2]:
            if row[0] not in uniprot_id_map:
                uniprot_id_map[row[0]] = set()
            uniprot_id_map[row[0]].add(row[1])

    thresholds = {
        column: threshold
        for column, threshold in {
            "neighborhood": neighborhood,
            "neighborhood_transferred": neighborhood_transferred,
            "fusion": fusion,
            "cooccurence": cooccurence,
            "homology": homology,
            "coexpression": coexpression,
            "coexpression_transferred": coexpression_transferred,
            "experiments": experiments,
            "experiments_transferred": experiments_transferred,
            "database": database,
            "database_transferred": database_transferred,
            "textmining": textmining,
            "textmining_transferred": textmining_transferred,
        }.items() if threshold
    }
    thresholds["combined_score"] = combined_score

    primary_accession = uniprot.get_primary_accession(taxon_identifier)

    for row in download.tabular_txt(
            "https://stringdb-static.org/download/protein.physical.links.full.v{version}/{taxon_identifier}.protein.links.full.v{version}.txt.gz"
            .format(taxon_identifier=taxon_identifier, version=version)
            if physical else
            "https://stringdb-static.org/download/protein.links.full.v{version}/{taxon_identifier}.protein.links.full.v{version}.txt.gz"
            .format(taxon_identifier=taxon_identifier, version=version),
            delimiter=" ",
            header=0,
            usecols=["protein1", "protein2"] + list(thresholds.keys()),
    ):
        if all(row[column] / 1000 >= thresholds[column]
               for column in thresholds):
            for interactor_a in uniprot_id_map.get(row["protein1"], set()):
                for interactor_b in uniprot_id_map.get(row["protein2"], set()):
                    for primary_interactor_a in primary_accession.get(
                            interactor_a, {interactor_a}):
                        for primary_interactor_b in primary_accession.get(
                                interactor_b, {interactor_b}):
                            yield (primary_interactor_a, primary_interactor_b)


def get_protein_protein_interactions(neighborhood=0.0,
                                     neighborhood_transferred=0.0,
                                     fusion=0.0,
                                     cooccurence=0.0,
                                     homology=0.0,
                                     coexpression=0.0,
                                     coexpression_transferred=0.0,
                                     experiments=0.0,
                                     experiments_transferred=0.0,
                                     database=0.0,
                                     database_transferred=0.0,
                                     textmining=0.0,
                                     textmining_transferred=0.0,
                                     combined_score=0.0,
                                     physical=False,
                                     taxon_identifier=9606,
                                     version=11.5):
    uniprot_id_map = {}
    for row in download.tabular_txt(
            "https://stringdb-static.org/download/protein.aliases.v{version}/{taxon_identifier}.protein.aliases.v{version}.txt.gz"
            .format(taxon_identifier=taxon_identifier, version=version),
            delimiter="\t",
            skiprows=1,
            usecols=[0, 1, 2],
    ):
        if "UniProt_AC" in row[2]:
            if row[0] not in uniprot_id_map:
                uniprot_id_map[row[0]] = set()
            uniprot_id_map[row[0]].add(row[1])

    thresholds = {
        column: threshold
        for column, threshold in {
            "neighborhood": neighborhood,
            "neighborhood_transferred": neighborhood_transferred,
            "fusion": fusion,
            "cooccurence": cooccurence,
            "homology": homology,
            "coexpression": coexpression,
            "coexpression_transferred": coexpression_transferred,
            "experiments": experiments,
            "experiments_transferred": experiments_transferred,
            "database": database,
            "database_transferred": database_transferred,
            "textmining": textmining,
            "textmining_transferred": textmining_transferred,
        }.items() if threshold
    }
    thresholds["combined_score"] = combined_score

    primary_accession = uniprot.get_primary_accession(taxon_identifier)

    for row in download.tabular_txt(
            "https://stringdb-static.org/download/protein.physical.links.full.v{version}/{taxon_identifier}.protein.links.full.v{version}.txt.gz"
            .format(taxon_identifier=taxon_identifier, version=version)
            if physical else
            "https://stringdb-static.org/download/protein.links.full.v{version}/{taxon_identifier}.protein.links.full.v{version}.txt.gz"
            .format(taxon_identifier=taxon_identifier, version=version),
            delimiter=" ",
            header=0,
            usecols=["protein1", "protein2"] + list(thresholds.keys()),
    ):
        if all(row[column] / 1000 >= thresholds[column]
               for column in thresholds):
            for interactor_a in uniprot_id_map.get(row["protein1"], set()):
                for interactor_b in uniprot_id_map.get(row["protein2"], set()):
                    for primary_interactor_a in primary_accession.get(
                            interactor_a, {interactor_a}):
                        for primary_interactor_b in primary_accession.get(
                                interactor_b, {interactor_b}):
                            yield (primary_interactor_a, primary_interactor_b,
                                   row["combined_score"] / 1000)
