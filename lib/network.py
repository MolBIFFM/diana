import csv
import gzip
import itertools
import logging
import math
import os
import pathlib
import random
import sys
import urllib.request
import zipfile
from concurrent import futures
from copy import deepcopy

import networkx as nx
import numpy as np
import pandas as pd
import progressbar
import scipy.cluster.hierarchy
from lxml import etree as ET
from scipy import stats

from lib.interaction_database import InteractionDatabase
from lib.posttranslational_modification_data import PosttranslationalModificationData
from lib.style import (CHART_POSITION, COLOR_MAP, PTM_TO_HEX, SHAPE_MAP, STYLE)

DPI = 300


def shortest_paths(source):
    """
    compute the shortest paths from a particular node "source" to all other nodes using breadth first traversal and
    return (partial) edge betweenness centrality of all edges in the network derived from these shortest paths
    """
    # dictionary associating edges with number of shortest paths from source that they are part of
    edge_betweeness_from_source = {}

    # dictionary associating nodes with shortest path length from source
    # access the network in global scope to reduce multiprocessing overhead
    dis = {node: float("inf") for node in network.nodes}

    # dictionary associating nodes all shortest paths from source
    paths_to_node = {node: [] for node in network.nodes}

    #  associate all neighbors of source with a single shortest paths between them an source
    for v in network[source]:
        dis[v] = 1
        paths_to_node[v] = [[source, v]]

    # u, v are the nodes joined by the edge traversed
    for u, v in nx.algorithms.traversal.edgebfs.edge_bfs(network, source):
        if dis[u] + 1 < dis[v]:
            # a shorter path from source to v exists via u
            dis[v] = dis[u] + 1
            # paths to v updated based on paths to u
            paths_to_node[v] = [path + [v] for path in paths_to_node[u]]
        elif dis[u] + 1 == dis[v]:
            # equally short paths from source to v exist via u
            # paths via u are added to the shortest paths from source to v
            paths_to_node[v].extend([path + [v] for path in paths_to_node[u]])

    # derive edge betweeness centrality based on all shortest paths from source
    for paths in paths_to_node.values():
        for path in paths:
            # exclude first and last edge in a path,
            # i.e. the edges incident to the nodes the shortest path connects
            for i in range(1, len(path) - 2):
                if (path[i], path[i + 1]) not in edge_betweeness_from_source:
                    edge_betweeness_from_source[(path[i], path[i + 1])] = 0.0
                    edge_betweeness_from_source[(path[i + 1], path[i])] = 0.0
                # consider that there may be multiple shortest paths between two nodes
                edge_betweeness_from_source[(path[i],
                                             path[i + 1])] += 1 / len(paths)
                edge_betweeness_from_source[(path[i + 1],
                                             path[i])] += 1 / len(paths)
    return edge_betweeness_from_source


def xml_from_dict(node, parent_node=None):
    """Create an ElementTree node form an applicable dictionary """
    if parent_node is None:
        node_in_xml = ET.Element(node["tag"], node["attrib"])
    else:
        node_in_xml = ET.SubElement(
            parent_node,
            node["tag"],
            node["attrib"],
        )
    for child in node["children"]:
        # recursively proceed with children,
        # attaching them to the current node as parent
        xml_from_dict(child, node_in_xml)
    if parent_node is None:
        return node_in_xml
    else:
        return node["tag"]


def column_data_availability(p, u):
    """return encoding for availability of posttranslational modification data"""
    if not math.isnan(p) and not math.isnan(u):
        return "B"
    elif not math.isnan(p):
        return "P"
    elif not math.isnan(u):
        return "U"
    else:
        return "N"


def ptm_type(p_mean, u_mean, ptm_sig=1.0):
    """
    return encoding for the type of deviation in posttranslational modification
    """

    # interpret lacking data at this time of measurement to represent
    # that there is no deviation
    p_mean = 0.0 if math.isnan(p_mean) else p_mean
    u_mean = 0.0 if math.isnan(u_mean) else u_mean

    if p_mean > 0.0 and u_mean < 0.0:
        return "P_UP_U_DOWN"
    elif p_mean < 0.0 and u_mean > 0.0:
        return "P_DOWN_U_UP"
    elif p_mean == 0.0 and u_mean == 0.0:
        return "MID"
    else:
        regulation = [p_mean, u_mean]
        if min(regulation) >= 0.0:
            if any([regulation[i] >= ptm_sig[i] for i in range(2)]):
                return "UP"
            else:
                return "MID_UP"
        elif max(regulation) <= 0.0:
            if any([regulation[i] <= -ptm_sig[i] for i in range(2)]):
                return "DOWN"
            else:
                return "MID_DOWN"


class PPINetworkStyle(object):
    """
    Class to assemble style specifications for Cytoscape based on user configuration

    """

    def __init__(self, configuration):
        # obtain applicable portion of user configuration
        self.configuration = configuration.get(
            "POST-TRANSLATIONAL MODIFICATION DATA")
        # Assemble an ElementTree form a dictionary stored separately
        # STYLE contains a default specification whose attributes are independent from user-configuration
        self.tree = ET.ElementTree(xml_from_dict(STYLE))

        # determine the uninon of all times of measurement set in the user configuration
        # set for individual types of posttranslational modification
        time_steps = []
        for modification in self.configuration.values():
            time_steps.extend(modification.get("TIME STEPS"))
        time_steps = sorted(set(time_steps))

        # duplicate the default style
        # such that one is available for each time of measurement
        root = self.tree.getroot()
        for _ in range(len(time_steps) - 1):
            root.append(deepcopy(root[0]))

        # add additional attributes
        # specific to each time of measurement and type of modification
        for i, visualStyle in enumerate(root):
            # set style name
            visualStyle.attrib["name"] = "PTM_PPI_Network_{}_min".format(
                time_steps[i])

            visualStyleNode = visualStyle.find("node")
            for j, modification in enumerate(self.configuration, 1):
                if time_steps[i] in self.configuration[modification].get(
                        "TIME STEPS"):
                    # add specification of node data represented by a bar chart
                    ET.SubElement(
                        visualStyleNode, "visualProperty", {
                            "default":
                                self.chart_configuration(
                                    modification, time_steps[i]),
                            "name":
                                "NODE_CUSTOMGRAPHICS_{}".format(j)
                        })
                    # add specification of bar chart position relative to node
                    ET.SubElement(
                        visualStyleNode, "visualProperty", {
                            "default": self.chart_position(modification),
                            "name": "NODE_CUSTOMGRAPHICS_POSITION_{}".format(j)
                        })

            # add specification for data availability representation by node shape
            # default value followed by discrete mapping between identifier and Cytoscape attribute
            map_parent_shape = ET.SubElement(visualStyleNode, "visualProperty",
                                             {
                                                 "default": "HEXAGON",
                                                 "name": "NODE_SHAPE"
                                             })
            map_node_shape = ET.SubElement(
                map_parent_shape, "discreteMapping", {
                    "attributeName": "DATA_{}".format(time_steps[i]),
                    "attributeType": "string"
                })
            for child in SHAPE_MAP:
                ET.SubElement(map_node_shape, child["tag"], child["attrib"])

            # add specification for representation of the type of deviation by node color
            # default value followed by discrete mapping between identifier and Cytoscape attribute
            map_parent_fill_color = ET.SubElement(visualStyleNode,
                                                  "visualProperty", {
                                                      "default": "#CCCCCC",
                                                      "name": "NODE_FILL_COLOR"
                                                  })
            map_node_fill_color = ET.SubElement(
                map_parent_fill_color, "discreteMapping", {
                    "attributeName": "REGULATION_{}".format(time_steps[i]),
                    "attributeType": "string"
                })
            for child in COLOR_MAP:
                ET.SubElement(map_node_fill_color, child["tag"],
                              child["attrib"])

    def export(self, file_name):
        """export the assembled ElementTree to an XML file"""
        # create directory implied by file_name
        if not os.path.exists(os.path.dirname(file_name)):
            pathlib.Path(os.path.dirname(file_name)).mkdir(parents=True,
                                                           exist_ok=True)
        self.tree.write("{}.styles.xml".format(os.path.splitext(file_name)[0]),
                        pretty_print=True)

    def chart_configuration(self, ptm, time_step):
        """
        assemble a Cytoscape bar chart specification 
        representing deviations at individual sites and specific times of measurement
        the bars scale in the range between twice ratios marking significant deviation (default: [-2.0:2.0])
        """
        return "org.cytoscape.BarChart:{{\"cy_range\":[{},{}],\"cy_autoRange\":false, \"cy_showRangeZeroBaseline\":true, \"cy_colors\":[{}],\"cy_dataColumns\":[{}], \"cy_showItemLabels\":true}}".format(
            2 * float(self.configuration[ptm].get("Significance", "-1.0")),
            2 * float(self.configuration[ptm].get("Significance", "1.0")),
            ",".join([
                "\"{}\"".format(self.configuration[ptm].get("Color", "#FFFFFF"))
                for j in range(self.configuration[ptm].get(
                    "SITES", 5))
            ]), ",".join([
                "\"{}_{}_{}\"".format(self.configuration[ptm].get("ID", "NA"),
                                      j + 1, time_step)
                for j in range(self.configuration[ptm].get(
                    "SITES", 5))
            ]))

    def chart_position(self, ptm):
        """
        obtain Cytoscape specification of bar chart position 
        according to configuration of a specific type of posttranslational modification
        """
        return CHART_POSITION.get(self.configuration[ptm].get(
            "Chart position", "R"))


class PPINetwork(object):
    """
    class to assemble and partition a protein-protein interaction network
    based on networkx.graph()
    """

    def __init__(self,
                 configuration=None,
                 graphml=None,
                 subgraph=None,
                 cluster_id=0,
                 score_histograms=None):
        self.logger = logging.getLogger("main")
        # import network from GraphML file
        if graphml:
            self.graph = nx.read_graphml(graphml)

            # index of cluster,
            # positive integer if network resembles a module obtained from modular decomposition,
            # 0 otherwise
            self.cluster_id = self.graph.graph["CLUSTER_ID"]

            # ID for types of posttranslational modification used for node annotation
            self.ptm_ids = self.graph.graph["PTM"].split(";")

            # ratios marking significant deviation
            self.ptm_sig = [
                float(s)
                for s in self.graph.graph["PTM_SIGNIFICANCE_LEVEL"].split(";")
            ]
            # union of times of measurement of all types of posttranslational modification
            self.time_steps = [
                int(ts) for ts in self.graph.graph["TIME_STEPS"].split(";")
            ]

        # network based on subgraph PPINetwork instance
        elif subgraph:
            self.graph = subgraph

            # analogous to GraphML import
            self.cluster_id = cluster_id
            self.graph.graph["CLUSTER_ID"] = self.cluster_id
            self.ptm_ids = self.graph.graph["PTM"].split(";")
            self.ptm_sig = [
                float(s)
                for s in self.graph.graph["PTM_SIGNIFICANCE_LEVEL"].split(";")
            ]
            self.time_steps = [
                int(ts) for ts in self.graph.graph["TIME_STEPS"].split(";")
            ]
            for protein in self.graph.nodes:
                self.graph.nodes[protein]["CLUSTER_ID"] = 0

        # assembly based on user configuration
        elif configuration:
            self.cluster_id = 0
            self.graph = nx.Graph()

            # analogous to GraphML import
            time_steps = set()
            attributes = set()
            self.ptm_ids = []
            self.ptm_sig = []

            for i, ptm in enumerate(
                    configuration.get("POST-TRANSLATIONAL MODIFICATION DATA",
                                      {})):
                # employ specialized class retrieve data on
                # an individual type of posttranslational modification
                # supply applicable portion of configuration
                modification = PosttranslationalModificationData(
                    ptm,
                    configuration["POST-TRANSLATIONAL MODIFICATION DATA"][ptm])
                self.ptm_ids.append(modification.id)
                self.ptm_sig.append(modification.significance)
                # add proteins retrieved from data set to the network as nodes
                for protein in modification.proteins:
                    self.graph.add_node(protein)
                    for time_step in modification.proteins[protein]:
                        time_steps.add(time_step)
                        # annotate node with ratios of deviation
                        # single largest ratio of deviation by absolute value
                        self.graph.nodes[protein]["{}_M_{}".format(
                            modification.id, time_step
                        )] = modification.proteins[protein][time_step]["MAX"]
                        attributes.add("{}_M_{}".format(modification.id,
                                                        time_step))

                        # sites exhibiting largest deviation
                        for i in range(
                                len(modification.proteins[protein][time_step]
                                    ["SITES"])):
                            self.graph.nodes[protein]["{}_{}_{}".format(
                                modification.id, i + 1,
                                time_step)] = modification.proteins[protein][
                                    time_step]["SITES"][i]
                            attributes.add("{}_{}_{}".format(
                                modification.id, i + 1, time_step))

            for protein in self.graph.nodes:
                # ID signifies ID of module that a particular protein was assigned to
                # by modular decomposition
                self.graph.nodes[protein]["CLUSTER_ID"] = 0
                # assign NaN to all combinations of modification ratio and time steps
                # for which no data is available
                for attribute in attributes:
                    if attribute not in self.graph.nodes[protein]:
                        self.graph.nodes[protein][attribute] = float("nan")

            self.time_steps = sorted(time_steps)

            # assign identifiers for data availability and type of deviation
            # used for visual representation in Cytoscape
            for protein in self.graph.nodes:
                for time_step in self.time_steps:
                    self.graph.nodes[protein]["DATA_{}".format(
                        time_step
                    )] = column_data_availability(
                        self.graph.nodes[protein]["P_M_{}".format(time_step)],
                        self.graph.nodes[protein]["U_M_{}".format(time_step)])
                    self.graph.nodes[protein]["REGULATION_{}".format(
                        time_step
                    )] = ptm_type(
                        self.graph.nodes[protein]["P_M_{}".format(time_step)],
                        self.graph.nodes[protein]["U_M_{}".format(time_step)],
                        self.ptm_sig)

            # store attributes in graph,
            # such that they are exportable to and importable from a GraphML file
            self.graph.graph["PTM"] = ";".join(self.ptm_ids)
            self.graph.graph["PTM_SIGNIFICANCE_LEVEL"] = ";".join(
                [str(s) for s in self.ptm_sig])
            self.graph.graph["CLUSTER_ID"] = self.cluster_id
            self.graph.graph["TIME_STEPS"] = ";".join(
                [str(ts) for ts in self.time_steps])

            # employ class to retrieve protein-protein interactions
            # from individual databases
            # supply applicable portion of configuration
            databases = []
            interactions = {}
            for database_name in configuration.get("PPI DATABASES", {}):
                databases.append(database_name)
                database = InteractionDatabase(
                    database_name, configuration["PPI DATABASES"][database_name])

                # transfer interactions from database-specific dictionary
                # to general dictionary
                for interactants, score in database.interactions.items():
                    if interactants not in interactions:
                        interactions[interactants] = {}
                    interactions[interactants][database_name] = score

            if interactions:
                # count interactions including proteins not included in
                # the posttranslational modification data
                num_uncovered = 0
                # integrate interactions to graph
                for interactants, scores in interactions.items():
                    interactants = sorted(interactants)
                    if interactants[0] in self.graph and interactants[
                            1] in self.graph:
                        sc = [score for score in scores.values() if score]
                        # assign largest confidence among all databases
                        # in which an interaction is reported
                        conf = max(sc) if sc else float("nan")
                        self.graph.add_edge(interactants[0],
                                            interactants[1],
                                            CONFIDENCE=conf)
                    else:
                        num_uncovered += 1

            if self.logger:
                # report to log file the number of interactions discarded
                # because at least one of its participants was not found
                # in posttranslational modification data
                self.logger.info(
                    "INTERACTIONS DISCARDED: NO POSTTRANSLATIONAL MODIFICATION DATA: {}/{} ({}%)"
                    .format(
                        num_uncovered,
                        num_uncovered + self.graph.number_of_edges(),
                        (num_uncovered /
                        (num_uncovered + self.graph.number_of_edges()) *100)))

            # remove isolate nodes, i.e. proteins not interacting with any others
            num_proteins = self.graph.number_of_nodes()
            isolates = list(nx.isolates(self.graph))
            if self.logger:
                self.logger.info(
                    "PROTEINS DISCARDED: NO INTERACTION DATA: {}/{} ({}%)".
                    format(len(isolates), num_proteins,
                            (len(isolates) / num_proteins) * 100))
            self.graph.remove_nodes_from(isolates)

    def export(self, file_name):
        """export the network to a GraphML file"""
        # if the network represents a module
        # obtained from modular decomposition of a network
        # add suffix designating its module ID to the file name
        if self.cluster_id:
            file_name = "{}.{}.graphml".format(
                os.path.splitext(file_name)[0], self.cluster_id)
        else:
            file_name = "{}.graphml".format(os.path.splitext(file_name)[0])

        # create the target directory implied by file_name
        if not os.path.exists(os.path.dirname(file_name)):
            pathlib.Path(os.path.dirname(file_name)).mkdir(parents=True,
                                                           exist_ok=True)
        nx.write_graphml(self.graph, file_name)

    def topological_clustering(self,
                  search_bound=None,
                  num_communities=None,
                  output=None,
                  weight=1.0,
                  use_density=False,
                  num_processes=1):
        """cluster the network using the Girvan-Newman algorithm"""
        if not self.graph.edges:
            return

        # tuples of node sets representing successive decompositions
        # convert to hash value of node set for memory efficiency
        initial_components = [component for component in nx.connected_components(self.graph)]
        opt_decomposition = initial_components

        # modularities of each partition
        # implied by the separate connected components of the modified graph
        if use_density:
            modularities = [self.modularity_density(initial_components)]
        else:
            modularities = [self.modularity(initial_components)]


        # largest modularity marking the optimal decomposition
        max_modularity = modularities[0]

        # number of divisions of connected components
        # without improving/increasing modularity
        iter_no_imp_mod = 0

        if not (search_bound or num_communities):
            bar = progressbar.ProgressBar(max_value=len(self.graph.nodes),
                                          min_value=len(initial_components))
        if num_communities:
            bar = progressbar.ProgressBar(max_value=num_communities,
                                          min_value=len(initial_components))
        elif search_bound:
            bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength,
                                          min_value=len(initial_components))
        if not search_bound:
            search_bound = len(self.graph.nodes)

        print("COMMUNITY OPTIMIZATION:")
        bar.update(len(initial_components))
        # self.girvan_newman() yields tuples of node sets
        # once a connected component is divided by the removal of edges
        for i, communities in enumerate(self.girvan_newman(num_processes), len(initial_components)):
            bar.update(i+1)

            # use either modularity or modularity density
            if use_density:
                modularities.append(self.modularity_density(
                    communities, weight))
            else:
                modularities.append(self.modularity(communities, weight))

            # update largest modularity/index of optimal decomposition
            if modularities[-1] > max_modularity:
                max_modularity = modularities[-1]
                max_index = i
                opt_decomposition = communities

            # assess whether modularity has improved with respect to the previous decomposition
            if len(modularities) > 1 and modularities[-1] <= modularities[-2]:
                iter_no_imp_mod += 1
            else:
                iter_no_imp_mod = 0

            if self.logger:
                # report to log file the  modularity/modularity density of the current partition
                if use_density:
                    self.logger.info(
                        "\tMODULARITY DENSITY OF {} COMMUNITIES: {} ({}/{} ITERATIONS WITHOUT IMPROVEMENT)"
                        .format(len(communities), modularities[-1],
                                iter_no_imp_mod, search_bound))
                else:
                    self.logger.info(
                        "\tMODULARITY OF {} COMMUNITIES: {} ({}/{} ITERATIONS WITHOUT IMPROVEMENT)"
                        .format(len(communities), modularities[-1],
                                iter_no_imp_mod, search_bound))

            # assess whether search bounds have been reached
            # regarding the number of divisions of connected components
            # without improvement of modularity or
            # the number of separate connected components
            if (iter_no_imp_mod == search_bound) or (
                    num_communities and len(communities) == num_communities):
                bar.finish()
                break

        if not num_communities and iter_no_imp_mod < search_bound:
            bar.finish()

        if self.logger:
            # report to log file the number of modules obtained by clustering
            self.logger.info("NETWORK PARTITIONED INTO {} COMMUNITIES".format(
                len(opt_decomposition)))

        # yield individual modules as PPINetwork objects
        for i in range(len(opt_decomposition)):
            for protein in opt_decomposition[i]:
                self.graph.nodes[protein]["CLUSTER_ID"] = i+1
            if output:
                self.induced_subgraph(opt_decomposition[i], i+1).export(output)

    def induced_subgraph(self, nodes, cluster_id=None):
        """return a PPINetwork resembling a subgraph induced by nodes"""
        if cluster_id is None:
            cluster_id = self.cluster_id
        return PPINetwork(subgraph=self.graph.subgraph(nodes).copy(),
                          cluster_id=cluster_id)

    def girvan_newman(self, num_processes=1):
        """
        Girvan-Newman algorithm yielding node partitions 
        upon each division of a connected component
        """
        # avoid multiprocessing overhead
        # that passing an object to a process repeatedly would cause
        global network

        # copy graph to remove edges without altering network
        # omit annotations of nodes and edges (memory efficiency)
        network = nx.Graph()
        network.add_nodes_from(self.graph.nodes)
        network.add_edges_from(self.graph.edges)

        # dictionary associating an edges
        # with their respective edge betweeness centrality
        edge_betweeness = {}

        def update_edge_betweeness(nodes):
            # reset edge betweeness of all edges
            # using nodes this update is restricted to the connected component
            # in which the last edge removal occurred
            # the edge betweeness centrality of edges within other
            # connected components is not affected
            for u, v in nx.subgraph(network, nodes).edges:
                edge_betweeness[(u, v)] = 0.0
                edge_betweeness[(v, u)] = 0.0

            # multiprocessing to search all shortest paths from each node
            with futures.ProcessPoolExecutor(num_processes) as executor:
                for partial_eb in executor.map(shortest_paths, nodes):
                    # accumulate edge betweeness centrality computed for all
                    # source nodes separately
                    for edge in partial_eb:
                        edge_betweeness[edge] += partial_eb[edge]

        # initialize edge betweeness centrality
        update_edge_betweeness(network.nodes)

        num_connected_components = nx.algorithms.components.number_connected_components(network)
        # proceed removing edges as long as there is an edge
        # in the network with an edge betweeness centrality > 0.0
        while network.edges and max(edge_betweeness.values()):

            # remove edge maximizing edge betweeness centrality
            u, v = max(edge_betweeness, key=edge_betweeness.get)
            network.remove_edge(u, v)

            # effectively exclude edge
            edge_betweeness[(u, v)] = 0.0
            edge_betweeness[(v, u)] = 0.0

            # recalculate edge betweeness folling the removal of an edge
            # restricted to the connected component that the edge belonged to
            update_edge_betweeness(
                nx.algorithms.components.node_connected_component(network, u))
            if nx.algorithms.components.number_connected_components(
                    network) > num_connected_components:
                num_connected_components += 1
                update_edge_betweeness(
                    nx.algorithms.components.node_connected_component(
                        network, v))
                # produce a tuple of separate node sets of all connected components
                yield tuple([
                    component for component in
                    nx.algorithms.components.connected_components(network)
                ])

    def modularity(self, communities, w=1.0):
        """Compute the weighted modularity of a node partition"""

        def sigma(u, v):
            for community in communities:
                if u in community:
                    if v in community:
                        return True
                    else:
                        return False
                elif v in community:
                    if u in community:
                        return True
                    else:
                        return False

        Q = 0.0
        m = self.graph.number_of_edges()
        adj = nx.adjacency_matrix(self.graph)
        nodes = list(self.graph.nodes)
        for i, u in enumerate(nodes):
            for j, v in enumerate(nodes[i + 1:], i + 1):
                if sigma(u, v):
                    Q += 2 * (adj[i, j] -
                              w * self.graph.degree[u] * self.graph.degree[v] /
                              (2.0 * m))
        return Q / (2.0 * m)

    def modularity_density(self, communities, w=1.0):
        """Compute the weighted modularity density of a node partition"""
        D = 0.0
        for community in communities:
            deg = 0.0
            for node in community:
                for adj_node in self.graph[node]:
                    if adj_node in community:
                        deg += 2 * w
                    else:
                        deg -= 2 * (1 - w)
            D += deg / len(community)
        return D
