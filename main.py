#!/usr/bin/env python3
import argparse
import logging
import os
import pathlib

import yaml


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help="GraphML file to import PPINetwork object from")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        help="file destination to export PPINetwork object to")
    parser.add_argument(
        "-c",
        "--configuration",
        type=str,
        help="configuration file specifying the data sets to build a network from")
    parser.add_argument("-l",
                        "--log",
                        type=str,
                        default=None,
                        help="file destination of log file")
    parser.add_argument(
        "-g",
        "--girvan-newman",
        default=False,
        action="store_true",
        help=
        "perform topological clustering of the network (Girvan-Newman algorithm)"
    )
    parser.add_argument("-n",
                        "--communities",
                        type=int,
                        default=None,
                        help="Topological clustering: limit on number of separate connected components to produce")
    parser.add_argument(
        "-b",
        "--bound",
        type=int,
        default=0,
        help=
        "Topological clustering: limit on number of cluster divisions without partition quality improvement")
    parser.add_argument("-w",
                        "--weight",
                        type=float,
                        default=1.0,
                        help="weight parameter for modularity and modularity density")
    parser.add_argument(
        "-d",
        "--density",
        default=False,
        action="store_true",
        help=
        "Topological clustering: modularity density will be used instead of modularity to assess partition quality"
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=1,
        help=
        "Topological clustering: number of processes to use for edge betweenness centrality computation")
    args = parser.parse_args()

    if args.log:
        if not os.path.exists(os.path.dirname(args.log)):
            pathlib.Path(os.path.dirname(args.log)).mkdir(parents=True,
                                                          exist_ok=True)
        logging.basicConfig(filename=args.log,
                            filemode="w",
                            format="%(asctime)s %(levelname)s\t%(message)s",
                            level=logging.INFO)

    logger = logging.getLogger("main")
    logger.propegate = False

    from lib.network import PPINetwork
    from lib.network import PPINetworkStyle

    if args.configuration:
        config = yaml.safe_load(open(args.configuration, "r"))
        if args.output:
            PPINetworkStyle(configuration=config).export(args.output)
    else:
        config = {}

    if args.input:
        network = PPINetwork(graphml=args.input)
    elif args.configuration:
        network = PPINetwork(configuration=config)

    if args.output:
        network.export(args.output)

    if args.girvan_newman:
        network.topological_clustering(args.bound, args.communities, args.output, 
            args.weight, args.density, args.processes)

        if args.output:
            network.export(args.output)


if __name__ == "__main__":
    main()
