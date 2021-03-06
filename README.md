# StringBuilder

A wrapper of STRING-DB's web api to build networks and GO tables from a given input of genes.


# Usage

## From the terminal
```bash

# can be called simply with a plain text file of genes
./sb.py -i <name>.txt -o <prefix>

# with additional arguments (extend network by 30 nodes, color by evidence, highres image)
./sb.py -i <name>.txt -o <prefix> \
    -n 30 \
    -f evidence \
    -r high

# requests network as tsv (without extension)
./sb.py -i <name>.txt -o <prefix> --network
```

## From python
```python3

# assemble gene list
genes = set([i.strip() for i in open("genes.txt", 'r')])

# initialize for a given set of genes
sb = StringBuilder(genes, prefix="out")

# request the network image and save it
sb.get_image(save=True)

# request the GO enrichment and save the dataframe
sb.get_functional_enrichment(save=True)

# request the network with additional arguments
sb.get_image(n_nodes=30, flavor="evidence", save=True)

```
