import networkx as nx


class ProteinProteinInteractionNetwork(nx.Graph):
    def __init__(self):
        super(ProteinProteinInteractionNetwork, self).__init__()

    from pipeline.data.excel import add_proteins_from_excel

    from pipeline.database.IntAct import add_interactions_from_IntAct
    from pipeline.database.BioGRID import add_interactions_from_BioGRID
    from pipeline.database.STRING import add_interactions_from_STRING
