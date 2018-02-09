# SIMLR
# Language: R
# Dependency: Requires SIMLR package, and specifically references projsplx_R.so
#             SIMLR is available at: https://github.com/BatzoglouLabSU/SIMLR
#             This file can be placed in the PluMA root directory
# Input: prefix (for multiple network files)
# Output: prefix (for CSV and NOA file for clusters) 
# Tested with: PluMA 1.0, R 3.2.5

PluMA plugin to cluster a network using the Single-Cell Interpretation
via Multiple Kernel Learning (SIMLR, Wang et al. 2017).

The plugin accepts input as a CSV file where rows and columns each
represent nodes, and entry (i, j) corresponds to the weight of the edge
from node i to node j.

The plugin will generate clusters in two formats: NOA and CSV.  The user
provides the prefix for the output files in the configuration file, which
will then generate files prefix.csv and prefix.noa.

The cluster CSV file will take the following format: 

"",     "x"
"1",Family.Lachnospiraceae.0091
"2",Class.Clostridia.0003
"3",Family.Lachnospiraceae.0093
"4",Family.Lachnospiraceae.0094
"",     "x"
"1",Family.Lachnospiraceae.0006
"2",Lactobacillus.0002
"3",Family.Porphyromonadaceae.0009
"4",Family.Lachnospiraceae.0008
"5",Family.Lachnospiraceae.0009
"6",Lactobacillus.0003
"7",Oscillibacter.0002
"8",Family.Ruminococcaceae.0001
"9",Family.Porphyromonadaceae.0016
"10",Atopostipes.0001
"11",TM7_genus_incertae_sedis.0004
"12",Clostridium_XlVa.0003
"13",Family.Lachnospiraceae.0106
"",     "x"
"1",Kingdom.Bacteria.0007
"2",Order.Clostridiales.0017
"3",Phylum.Firmicutes.0003
"4",Family.Ruminococcaceae.0024
"5",Family.Lachnospiraceae.0092

Note that each cluster is separated by a line "", "x"

The other output file is a NOde Attribute (NOA) file for Cytoscape, which
is just a table with each node name and the cluster to which it belongs:

Name    Cluster
Family.Porphyromonadaceae.0001  23
Family.Porphyromonadaceae.0002  23
Family.Porphyromonadaceae.0003  18
Bacteroides.0001        6
Alistipes.0001  6
Family.Porphyromonadaceae.0004  6
Barnesiella.0001        9
Family.Porphyromonadaceae.0005  6
Lactobacillus.0001      18
Family.Lachnospiraceae.0001     18
Family.Porphyromonadaceae.0006  16
Family.Lachnospiraceae.0002     18
Kingdom.Bacteria.0001   21

This cluster then becomes an attribute of each node in Cytoscape,
and can be used for further downstream analysis and visualization.
