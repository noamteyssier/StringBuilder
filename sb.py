#!/usr/bin/env python3


import argparse
import requests
import pandas as pd
from time import sleep
import sys
import io


class StringBuilder:

    def __init__(self, genes, prefix="out"):
        self.genes = genes
        self.prefix = prefix
        self.api_url = "https://string-db.org/api"
        self.request_url = "{}/{}/{}"

    def call(self, output_format, method, params, name=None):
        """
        request information through API
        """
        if name:
            print(
                "Calling String:\n\tmethod: {}\n\tname: {}".format(
                    method, name
                    ),
                file=sys.stderr
            )

        response = requests.post(
            self.request_url.format(self.api_url, output_format, method),
            data=params
        )
        sleep(0.2)

        return response

    def write_image(self, response):
        """
        writes a network png to file
        """
        file_name = "{}_network.png".format(self.prefix)
        print("Writing : {}".format(file_name), file=sys.stderr)

        with open(file_name, 'wb') as fh:
            fh.write(response.content)

    def write_go(self, go_frame):
        """
        writes go frame to file
        """
        file_name = "{}_GO.tsv".format(self.prefix)
        print("Writing : {}".format(file_name), file=sys.stderr)

        go_frame.to_csv(
            file_name, sep="\t", index=False
        )

    def write_net(self, net_frame):
        """
        writes network to tab-delim file
        """
        file_name = "{}_net.tsv".format(self.prefix)
        print("Writing : {}".format(file_name), file=sys.stderr)

        net_frame.to_csv(
            file_name, sep="\t", index=False
        )

    def write_map(self, map_frame):
        """
        writes gene-name map to tab-delim file
        """

        file_name = "{}_map.tsv".format(self.prefix)
        print("Writing : {}".format(file_name), file=sys.stderr)

        map_frame.to_csv(
            file_name, sep="\t", index=False
        )

    def get_image(self, genes=None, save=False, n_nodes=10,
                  flavor="confidence", resolution='low'):
        """
        prepares API for network image call and passes request
        """

        if isinstance(genes, type(None)):
            genes = self.genes

        if resolution == 'low':
            output_format = "image"
        else:
            output_format = "highres_image"

        method = "network"

        params = {
            "identifiers": "%0d".join(genes),
            "species": 9606,
            "add_white_nodes": n_nodes,
            "network_flavor": flavor,
            "caller_identity": "Kampmann Lab"
        }

        response = self.call(
            output_format, method, params,
            name="Network Image"
            )

        if save:
            self.write_image(response)

        return response

    def get_extended_nodes(self, genes=None, n_nodes=10):
        """
        requests all nodes found + extended network
        """

        if isinstance(genes, type(None)):
            genes = self.genes

        output_format = "tsv-no-header"
        method = "network"

        params = {
            "identifiers": "%0d".join(genes),
            "species": 9606,
            "add_white_nodes": 10,
            "caller_identity": "Kampmann Lab"
        }

        response = self.call(
            output_format, method, params
            )

        all_nodes = set()
        for line in response.text.strip().split("\n"):
            [all_nodes.add(i) for i in line.split("\t")[2:4]]

        return list(all_nodes)

    def get_functional_enrichment(self, genes=None, save=False):
        """
        requests functional enrichment of extended network
        """

        if isinstance(genes, type(None)):
            genes = self.genes

        output_format = "json"
        method = "enrichment"

        params = {
            "identifiers": "%0d".join(genes),
            "species": 9606,
            "called_identity": "Kampmann Lab"
        }

        response = self.call(
            output_format, method, params,
            name="GO Enrichment"
            )
        go_frame = pd.read_json(response.text)

        if save:
            self.write_go(go_frame)

        return go_frame

    def get_network(self, genes=None, save=False, n_nodes=0):
        """
        requests inter-node network without network extension
        """

        if isinstance(genes, type(None)):
            genes = self.genes

        method = "network"
        output_format = "tsv"
        params = {
            "identifiers": "%0d".join(genes),
            "species": 9606,
            "caller_identity": "Kampmann Lab",
            "add_white_nodes": n_nodes
        }

        response = self.call(
            output_format, method, params,
            name="Network TSV"
            )

        # print(response.text)
        net_frame = pd.read_csv(
            io.StringIO(response.content.decode('utf-8')),
            sep="\t"
            )

        if save:
            self.write_net(net_frame)

        return net_frame

    def get_identifiers(self, genes=None, save=False):
        """
        maps gene names to identifiers
        """

        if isinstance(genes, type(None)):
            genes = self.genes

        method = "get_string_ids"
        output_format = "tsv"
        params = {
            "identifiers": "%0d".join(genes),
            "species": 9606,
            "limit": 1,
            "echo_query": 1,
            "caller_identity": "Kampmann Lab"
        }

        response = self.call(
            output_format, method, params,
            name="Network TSV"
            )

        map_frame = pd.read_csv(
            io.StringIO(response.content.decode('utf-8')),
            sep="\t"
            )

        if save:
            self.write_map(map_frame)

        return map_frame


def read_genes(txt):
    with open(txt, "r") as f:
        return set([i.strip() for i in f])


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "-i", "--input", required=True,
        help="input text file of genes to process"
        )
    p.add_argument(
        "-o", "--output_prefix", required=False, default="out",
        help="Output prefix to prepend to recovered data"
    )
    p.add_argument(
        "-n", "--nodes", required=False, default=10,
        help="Number of nodes to extend the network (default=10)"
    )
    p.add_argument(
        "-f", "--flavor", required=False, default="confidence",
        help="flavor to color interactions - choices=confidence/evidence"
    )
    p.add_argument(
        "-r", "--resolution", required=False, default="low",
        help="resolution of image to request from string (low/high)"
    )
    p.add_argument(
        "--network", required=False, action='store_true',
        help="only process network to tsv"
    )
    args = p.parse_args()
    return args


def main():
    args = get_args()
    genes = read_genes(args.input)

    sb = StringBuilder(genes, prefix=args.output_prefix)

    if args.network:
        sb.get_network(save=True, n_nodes=args.nodes)
        sb.get_identifiers(save=True)

    else:
        sb.get_image(
            n_nodes=args.nodes,
            flavor=args.flavor,
            resolution=args.resolution,
            save=True
            )

        extended_network = sb.get_extended_nodes(
            n_nodes=args.nodes
        )

        sb.get_functional_enrichment(extended_network, save=True)
        sb.get_identifiers(genes, save=True)


if __name__ == '__main__':
    main()
