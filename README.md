# Olave_etal2018-JEB
DESCRIPTION:
These functions are written to calculate the Extra Lineage Contribution (XLC) statistics by using the PhyloNet package (Than and Nakhleh 2009) .
The original article describing the method:
Olave, M., Avila, L.J., Sites, J.W. Jr and Morando, M (2018). Hybridization could be a common phenomenon within the highly diverse lizard genus Liolaemus. Journal of Evolutionary Biology.

REQUIREMENTS
functions.R --> github.com/melisaolave/Olave_etal2018-JEB
ape, foreach, doMC installed in R
Phylonet v2.3 (Than and Nakhleh 2009): http://old-bioinfo.cs.rice.edu/phylonet/index.html # NOTE: it might not work with a different PhyloNet version!

Files you should have ready
1. a species tree in newick format. See example sp-tree.tre
2. a set of gene trees in newick format. See examples genetree1.tre - genetree10.tre
3. an Imap file containing the individual-species associations. See Imap.txt as example. Columns separated by tab, and columns names: traits and species.

EXECUTION
See the tutorial R script for an example and explanation.

Any question, please refer to:
Email: melisa.olave@uni-konstanz.de
