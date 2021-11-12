from download import download

UNIPROT_SWISSPROT = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz"
ID_MAP = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/{organism}_idmapping.dat.gz"

ORGANISM = {
    9606: "HUMAN_9606",
}


def get_swissprot_entries():
    accessions, entry_gene_name, entry_protein_name = [], {}, {}
    rec_name, taxon_identifier = False, 0
    for line in download.txt(UNIPROT_SWISSPROT):
        if not line.strip():
            continue

        if line.split(maxsplit=1)[0] == "AC":
            if len(line.split(maxsplit=1)) == 1:
                continue

            accessions.extend(
                line.split(maxsplit=1)[1].rstrip(";").split("; "))

        elif line.split(maxsplit=1)[0] == "GN":
            if len(line.split(maxsplit=1)) == 1:
                continue

            for entry in line.split(maxsplit=1)[1].rstrip(";").split("; "):
                if "=" in entry:
                    entry_gene_name[entry.split("=")[0]] = (
                        entry.split("=")[1].split("{")[0].rstrip())

        elif line.split(maxsplit=1)[0] == "DE":
            if len(line.split(maxsplit=1)) == 1:
                continue

            if line.split(maxsplit=1)[1].split(":", 1)[0] == "RecName":
                entries = (line.split(maxsplit=1)[1].split(
                    ":", 1)[1].lstrip().rstrip(";").split("; "))
                rec_name = True
            elif line.split(maxsplit=1)[1].split(":", 1)[0] == "AltName":
                entries = []
                rec_name = False

            elif line.split(maxsplit=1)[1].split(":", 1)[0] == "Contains":
                entries = []
                rec_name = False

            elif line.split(maxsplit=1)[1].split(":", 1)[0] == "Flags":
                entries = []
                rec_name = False

            elif rec_name:
                entries = line.split(maxsplit=1)[1].rstrip(";").split("; ")

            for entry in entries:
                if "=" in entry:
                    if entry.split("=")[0] not in entry_protein_name:
                        entry_protein_name[entry.split("=")[0]] = (
                            entry.split("=")[1].split("{")[0].rstrip())

        elif line.split(maxsplit=1)[0] == "OX":
            if len(line.split(maxsplit=1)) == 1:
                continue

            if line.split(
                    maxsplit=1)[1].split(";")[0].split("=")[0] == "NCBI_TaxID":
                if (line.split(maxsplit=1)[1].split(";")[0].split("=")
                    [1].split("{")[0].isnumeric()):
                    taxon_identifier = int(
                        line.split(maxsplit=1)[1].split(";")[0].split("=")
                        [1].split("{")[0])

        elif line == "//":
            yield (
                accessions,
                entry_gene_name.get("Name", "NA"),
                entry_protein_name.get("Full", "NA"),
                taxon_identifier,
            )

            accessions.clear()
            entry_gene_name.clear()
            entry_protein_name.clear()
            rec_name, taxon_identifier = False, 0


def get_primary_accession(taxon_identifier=9606, proteins=set()):
    primary_accession = {}
    for accessions, _, _, taxon in get_swissprot_entries():
        if taxon == taxon_identifier and (not proteins
                                          or accessions[0] in proteins):
            for accession in accessions:
                if accession not in primary_accession:
                    primary_accession[accession] = set()
                primary_accession[accession].add(accessions[0])

    return primary_accession


def get_id_map(database,
               taxon_identifier=9606,
               proteins=set(),
               identifier_type=str):
    if taxon_identifier not in ORGANISM:
        return {}

    uniprot = {}
    for row in download.tabular_txt(
            ID_MAP.format(organism=ORGANISM[taxon_identifier]),
            delimiter="\t",
            usecols=[0, 1, 2]):
        if row[1] == database and (not proteins or row[0] in proteins):
            if identifier_type(row[2]) not in uniprot:
                uniprot[identifier_type(row[2])] = set()

            uniprot[identifier_type(row[2])].add(row[0])

    return uniprot
