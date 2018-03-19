# When using, please cite: Olave, M., Avila, L.J., Sites, J.W. Jr and Morando, M (2018). Hybridization could be a common phenomenon within the highly diverse lizard genus Liolaemus. Journal of Evolutionary Biology.

# This is a tutorial to obtain a matrix of Extra Lineage Contribution statistic (XLC) as decribed by Olave et al. (2018; JEB)
# any question: melisa.olave@uni-konstanz.de

# REQUIREMENTS
## functions.R --> github.com/melisaolave/Olave_etal2018-JEB
## ape, foreach, doMC installed in R
## Phylonet v2.3 (Than and Nakhleh 2009): http://old-bioinfo.cs.rice.edu/phylonet/index.html # NOTE: it might not work with a different PhyloNet version!

# Files you should have ready
## a species tree in newick format. See example sp-tree.tre
## A set of gene trees in newick format. See examples genetree1.tre - genetree10.tre
## an Imap file containing the individual-species associations. See Imap.txt as example. Columns separated by tab, and columns names: traits and species.

###### Below starts the tutorial  #########
# make sure you have installed the ape, foreach, doMC libraries

# download the functions from github.com/melisaolave/Olave_etal2018-JEB. Place them into your working directory folder. 
# In this tutorial we assume that all the files are in:
wd <- "/Users/Desktop/tutorial/examples"

#thus, set working directory
setwd(wd);

# and read the functions:
source("functions.R");

# The folder examples has a species tree (sp-tree.tre), a set of gene trees (genetree1.tre - genetree10.tre) and an Imap (Imap.txt).
# use the getXLC function to obtain the XLC matrix:
getXLC(wd=wd, genetreePattern=".tre$", sp.treeName="sp-tree.tre", ImapName="Imap.txt", 
       PhylonetPath="/Users/melisaolave/Downloads/Programas/Trabajo/phylonet_v2_3.jar", outputName="XLC.csv", cores=1);

# the attributes of this function:
  #wd : working directory. By default it will take the current working directory.
  #genetreePattern : a pattern to match all gene tree files using regular expression. The species tree name is removed from the list if it is also reached by the regular expression.
    # In the example, the $ symbol at the end of .tre$ is the regular expression form of saying that the .tre has to be found at the end of the file name.
    # By default it will match all files ending in .tre.
    # For example, here all the gene tree files end with .tre, thus we use this pattern. 
  # sp.treeName : this is the full name of the file containing the species tree in newick format. By default: sp-tree.tre
  # ImapName : this is the full name of the map file contaning the individual - species associations. By default: Imap.txt
  # PhylonetPath : a path to the PhyloNet program executable. E.g.: if the java executable is at the example folder = Users/Desktop/tutorial/examples/phylonet_v2_3.jar
    # If PhylonetPath = NULL, then by default it will try to find the executable phylonet_v2_3.jar
    # IMPORTANT only PhyloNet version 2.3. has been tested, other versions might not work.
  # outputName : name for the output files.
  # cores : it is possible to run the calculations in parallel. Change the number of cores you want to use. 
    # Calculations for each column (i.g. gene tree) in the XLC matrix can be done in a different core. Thus, the maximum efficiency can be reached by setting: number of cores = number of gene trees.
    # Larger number of cores than gene trees will not make sense.
    # Note that there is not progress report when running in parallel.
  

# OUTPUTS
# The outputs in this example are 
  # 1. absolute-XLC.csv : it summarized the absolute numbers of extra lineage (XL) counted per allele (or individual) per gene tree. The last row is the final count per gene tree (= XLi in our ecuation in section 2.3)
  # 2. allelesSummary.XLC.csv : some summary statistics of the extra lineages counted per allele (from values shown in absolute-XLC.csv). For example, for ind1 the minimum value of XL is 0 and the maximum is 1, with a mean = 0.2387.
  # 3. genetreeSummary.XLC.csv : some summary statistics of the extra lineages counted per gene tree (from values shown in absolute-XLC.csv).
  # 4. XLC.csv : the extra lineage contribution (XLC) matrix. The last column is the XLCj (as shown in our ecuation in section 2.3).



# Plot the XLC matrix
# We have also included the plotXLCmatrix to procude a similar plot the our Figure 4a. I recommend you to play with the values to understand what they do.
plotXLCmatrix(wd=wd, tableName="XLC.csv", reportInd=NULL, cex.txt=1,
                          oma=c(1, 0, 1, 0), mar=c(1, 1, 1, 1), min.xlim=-5, min.ylim=-3, max.ylim=4, 
                          pch=15, square.cex=0.75, x.square.sep=1, y.square.sep=1,
                          spBar=T, sp.color=c("firebrick3", "dodgerblue4", "darkolivegreen4", "darkorchid4"), x.spBar=2, 
                          height=4, width=2.5,
                          gradientRate=40, colorGradient=c("gold","orange","red"), XLCij.label=1,
                          colorGradient.XLCj=c("gray95","black"), XLCj.label=1,
                          max.XLC.pos=1.75, min.XLC.pos=1, x.scaleBar=2, y.scaleBar=-1, scaleBar.adj=0.2,
                          y.label.line = -1, x.label.line=0,
                          savePlot=T, plotName="XLCplot.pdf");

# the attributes of this function:
  #wd : working directory. By default it will take the current working directory.
  # tableName : name of the XLC matrix exported using the getXLC() function.
  # reportInd : a character vector can be provided with the individuals that want to be highlighthed in the matrix. Then an asteriks will be plotted on the right.
    # For example try reportInd = c("ind1", "ind5") to report the position of ind1 and ind5.
  # cex.txt : cex of text.
  # oma & mar : used in par() function.
  # min.xlim, min.ylim & max.ylim : they can be used to control the x and y limits.
  # pch : type of point for plotting the matrix.
  # square.cex : cex for the pch.
  # x.square.sep & y.square.sep : they are used to increase or decrease the separation of the points pch.
  # spBar : if TRUE, then a bar on the left will be drawn to highlight the different species.
  # sp.color : if spBAR is TRUE, then provide a character vector of color names to be used to each species bar
  # x.spBar : it controls the position in the x axis of the species bars
  # height  &  width : heoght and width of plot to be saved.
  # gradientRate : it is used for the gradienteRate for the colorRampPalette. Can be between 1 and 40.
  # colorGradient : a character vector with color names for the XLCij bar 
  # XLCij.label : it is used to adjust the position of the XLCij label
  # colorGradient.XLCj: a character vector with color names for the XLCj bar
  # XLCj.label : it is used to adjust the position of the XLCj label
  # max.XLC.pos & min.XLC.pos : they control the position in the x axis of the maximum and minimum values in both XLC scale bars
  # x.scaleBar & y.scaleBar : they control the position of both XLC scale bars
  # scaleBar.adj : it controls the position of both XLC scale bars
  # y.label.line & x.label.line : they control the position of the y and x axis labels (Alleles and Gene trees, respectively)
  # savePlot : if TRUE, then the plot is exported in pdf format
  # plotName : of savePlot is TRUE, then provide a name for the output plot