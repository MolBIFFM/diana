BIOGRID_ID_MAP_ZIP_ARCHIVE = "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-IDENTIFIERS-LATEST.tab.zip"
BIOGRID_ID_MAP = r"BIOGRID-IDENTIFIERS-[0-9]\.[0-9]\.[0-9]{3}\.tab\.txt"
BIOGRID_ZIP_ARCHIVE = "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip"
BIOGRID = r"BIOGRID-ORGANISM-{organism}-[0-9]\.[0-9]\.[0-9][0-9][0-9]\.tab3\.txt"
BIOGRID_MV_PHYSICAL_ZIP_ARCHIVE = "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-MV-Physical-LATEST.tab3.zip"
BIOGRID_MV_PHYSICAL = r"BIOGRID-MV-Physical-[0-9]\.[0-9]\.[0-9]{3}\.tab3\.txt"

CORUM_ZIP_ARCHIVE = (
    "http://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip"
)
CORUM = "allComplexes.txt"

INTACT_ZIP_ARCHIVE = (
    "ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip"
)
INTACT = "intact.txt"

REACTOME = "https://reactome.org/download/current/interactors/reactome.{organism}.interactions.tab-delimited.txt"

STRING_ID_MAP = "https://stringdb-static.org/download/protein.aliases.v11.0/{taxon_identifier}.protein.aliases.v11.0.txt.gz"
STRING = "https://stringdb-static.org/download/protein.links.full.v11.0/{taxon_identifier}.protein.links.full.v11.0.txt.gz"
STRING_PHYSICAL = "https://stringdb-static.org/download/protein.physical.links.full.v11.0/{taxon_identifier}.protein.links.full.v11.0.txt.gz"

UNIPROT_ID_MAP = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/{organism}_{taxon_identifier}_idmapping.dat.gz"
UNIPROT_SWISSPROT = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz"

ORGANISM = {
    9606: {BIOGRID: "Homo_sapiens", REACTOME: "homo_sapiens", UNIPROT_ID_MAP: "HUMAN"},
}
