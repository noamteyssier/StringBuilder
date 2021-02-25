#!/usr/bin/env python3


import argparse
import requests
import pandas as pd
from time import sleep


class StringBuilder:

    def __init__(self, genes, prefix="out"):
        self.genes = genes
        self.prefix = prefix
        self.api_url = "https://string-db.org/api"
        self.request_url = "{}/{}/{}"

    def call(self, output_format, method, params):
        """
        request information through API
        """

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
        with open(file_name, 'wb') as fh:
            fh.write(response.content)

    def write_go(self, go_frame):
        """
        writes go frame to file
        """
        file_name = "{}_GO.tsv".format(self.prefix)
        go_frame.to_csv(
            file_name, sep="\t", index=False
        )

    def get_image(self, genes=None, save=False, n_nodes=10,
                  flavor="confidence"):
        """
        prepares API for network image call and passes request
        """

        if isinstance(genes, type(None)):
            genes = self.genes

        output_format = "highres_image"
        method = "network"

        params = {
            "identifiers": "%0d".join(genes),
            "species": 9606,
            "add_white_nodes": n_nodes,
            "network_flavor": flavor,
            "caller_identity": "Kampmann Lab"
        }

        response = self.call(output_format, method, params)

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

        response = self.call(output_format, method, params)

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

        response = self.call(output_format, method, params)
        go_frame = pd.read_json(response.text)

        if save:
            self.write_go(go_frame)

        return go_frame


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
    args = p.parse_args()
    return args


def main():
    args = get_args()
    genes = read_genes(args.input)

    sb = StringBuilder(genes, prefix=args.output_prefix)

    sb.get_image(save=True)
    extended_network = sb.get_extended_nodes()
    sb.get_functional_enrichment(extended_network, save=True)


if __name__ == '__main__':
    main()
