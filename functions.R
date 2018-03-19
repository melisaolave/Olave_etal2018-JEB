getXLC <- function(wd=getwd(), genetreePattern=".tre$", sp.treeName="sp-tree.tre", ImapName="Imap.txt", PhylonetPath=NULL, outputName="XLC.csv", cores=1){
  require(ape);
  require(foreach);
  require(doMC);
  registerDoMC(cores);
  setwd(wd);
  Imap <- read.table(ImapName, sep="\t", header=T);
  gene.treeNames <- list.files(pattern=genetreePattern);
  gene.treeNames <- gene.treeNames[gene.treeNames != sp.treeName]; #removing species tree in case it was matched
  species <- unique(as.character(Imap$species));
  cat(length(species), "species and ", nrow(Imap)," individuals found in Imap file\n");
  sptree.newick <- scan(sp.treeName, what="character");
  sp <- as.character(unique(Imap$species));
  summary.table <- NULL
  cat("Gene trees found:", gene.treeNames, "\n")
  cat("Calculating extra lineages, be pacient\n");
  if(cores > 1){
    cat("Not progress report when running in parallel\n");
  }
  absolute.table <- foreach(i=1:length(gene.treeNames), .combine=cbind) %dopar%{
    xl.all <- sptree.vs.genetrees(wd=wd, PhylonetDir=PhylonetPath, tree.newick=sptree.newick, 
                                  genetreeNames=gene.treeNames[i], species=species, Imap=Imap, sp.treeName=sp.treeName);
    tree <- read.tree(gene.treeNames[i]);
    xlVec <- NULL;
    table <- data.frame(as.character(Imap$traits), as.character(Imap$species),stringsAsFactors=FALSE);
    table <- cbind(table, rep(NA, nrow(table)), stringsAsFactors=FALSE);
    table <- rbind(table, c("XLt", "-", xl.all), stringsAsFactors=FALSE);
    report.xl <- NULL;
    xlVec.plot <- NULL;
    for(k in 1:length(sp)){
      sp.tip <-  as.character(Imap[Imap$species == sp[k],1]);
      new.tree <- drop.tip(tree, sp.tip);
      tempName <- paste(gene.treeNames[i],"-dropTip", i, ".temp", sep="");
      write.tree(new.tree, tempName);
      xl.sp.allremoved <- sptree.vs.genetrees(wd=wd, PhylonetDir=PhylonetPath, tree.newick=sptree.newick, 
                                              genetreeNames=tempName, species=species, Imap=Imap, sp.treeName=sp.treeName);
      for(j in 1:length(sp.tip)){
        if(sp.tip[j] %in% tree$tip.label){
          tip.rm <-  setdiff(sp.tip, sp.tip[j]);
          new.tree <- drop.tip(tree, tip.rm);
          write.tree(new.tree, tempName);
          PhyloImapName <- paste("PhyloImap-",tempName,sep="");
          xl <- sptree.vs.genetrees(wd=wd, PhylonetDir=PhylonetPath, PhyloImapName=PhyloImapName,tree.newick=sptree.newick, 
                                    genetreeNames=tempName, species=species, Imap=Imap, sp.treeName=sp.treeName);
          XLCij <- xl-xl.sp.allremoved;
          table[table[,1] == sp.tip[j],3] <- XLCij;
        }
      }
    }
    file.remove(tempName);
    cat("Gene tree",i, "/", length(gene.treeNames), "done!\n");
    return(table[,3]);
  }
  table <- data.frame(as.character(Imap$traits), as.character(Imap$species),stringsAsFactors=FALSE);
  table <- rbind(table, c("XLt", "-"), stringsAsFactors=FALSE);
  table <- cbind(table, absolute.table);
  colnames(table) <- c("individuals", "species", gene.treeNames);
  absTableName <- paste("absolute-", outputName, sep="");
  write.table(table,absTableName, col.names=T, row.names=T, quote=F, sep="\t");
  table<- read.table(absTableName, sep="\t", header=T)
  ### getting XLC per allele
  new.table <- table[1:(nrow(table)-1),];
  new.table <- cbind(new.table, rep(NA, nrow(new.table)));
  ind.summaryTable <- NULL;
  for(i in 1:(nrow(table)-1)){
    XLCij.vec <- NULL;
    subtable <- table[i,3:ncol(table)];
    XLt.vec <- as.integer(table[nrow(table),3:ncol(table)])
    XLt.vec[XLt.vec == 0] <- 0.000000000001 # this is to prevent values = 0 and infinit results in the division below (note that it will anyway rounded to 4, and all small values will be = 0 again);
    XLCij.vec <- round(as.numeric(subtable/XLt.vec),4);
    no.missing <- is.na(subtable) == F;
    XLCj <- sum(as.integer(subtable[no.missing]));
    XLt <- sum(XLt.vec[no.missing]);
    min <- min(XLCij.vec[no.missing])
    ind.SummaryStats <- summary(XLCij.vec[no.missing]);
    ind.SummaryStats <- c(as.character(new.table[i,2], stringsAsFactors=FALSE),ind.SummaryStats, sd(XLCij.vec[no.missing]))
    ind.summaryTable <- rbind(ind.summaryTable, ind.SummaryStats);
    colnames(ind.summaryTable) <- c("species","Min", "1st Qu.", "Median", "Mean", "3rd Qu.","Max.", "sd");
    row.names(ind.summaryTable) <- as.character(new.table[1:i,1]);
    write.table(ind.summaryTable, paste("allelesSummary.",outputName, sep=""), col.names=T, row.names=T, quote=F, sep="\t");
    XLCij.vec <- c(XLCij.vec, round(XLCj/XLt, 4));
    new.table[i,3:ncol(new.table)] <- XLCij.vec;
  }
  colnames(new.table) <- c("individuals", "species", gene.treeNames, "XLCj");
  write.table(new.table,outputName, col.names=T, row.names=F, quote=F, sep="\t");
  
  ## getting summary stats per gene tree
  summary.table <- NULL
  for(i in 3:(ncol(new.table)-1)){
    summaryStats <- summary(as.numeric(new.table[is.na(new.table[,i]) == F,i]));
    summaryStats <- c(summaryStats, sd(as.numeric(new.table[is.na(new.table[,i]) == F,i])));
    summary.table <- rbind(summary.table, summaryStats);
    colnames(summary.table) <- c("Min", "1st Qu.", "Median", "Mean", "3rd Qu.","Max.", "sd");
    row.names(summary.table) <- gene.treeNames[1:(i-2)];
    write.table(summary.table, paste("genetreeSummary.",outputName, sep=""), col.names=T, row.names=T, quote=F, sep="\t");
  }
}

plotXLCmatrix <- function(wd=getwd(), tableName="XLC.csv", reportInd=NULL, cex.txt=1,
                          oma=c(1, 0, 1, 0), mar=c(1, 1, 1, 1), min.xlim=-5, min.ylim=-3, max.ylim=4, 
                          pch=15, square.cex=0.75, x.square.sep=1, y.square.sep=1,
                          spBar=F, sp.color=NULL, x.spBar=2, 
                          height=8.27, width=8.27/2.5,
                          gradientRate=40, colorGradient=c("white","orange","red"), XLCij.label=3,
                          colorGradient.XLCj=c("gray95","black"), XLCj.label=3,
                          max.XLC.pos=1.5, min.XLC.pos=1, x.scaleBar=3, y.scaleBar=-3, scaleBar.adj=0.2,
                          y.label.line = -1, x.label.line=0,
                          savePlot=T, plotName="XLCplot.pdf"){
  setwd(wd);
  table <- read.table(tableName, sep="\t", header=T);
  max.XLC <- max(table[,3:(ncol(table)-1)], na.rm=T);
  colfunction <- colorRampPalette(colorGradient);
  colorRange <- colfunction(gradientRate+1);
  max.XLCj <- max(table[,ncol(table)], na.rm=T);
  colfunction.XLCj <- colorRampPalette(colorGradient.XLCj);
  colorRange.XLCj <- colfunction.XLCj(gradientRate+1);
  dev.new(width=width, height=height);
  plot.new(); #clear all previous parameters
  quartz.options(height=7, width=10, dpi=72);
  par(oma = oma); #outer margins
  par(mar = mar);
  plot.window(xlim = c(min.xlim,(ncol(table)+1)), ylim = c(min.ylim,nrow(table)+max.ylim));
  y.coor <- rev(1:nrow(table));
  x.coor <- 1:(ncol(table)-2);
  x.coor[length(x.coor)] <- x.coor[length(x.coor)]+1;
  x.coor <- x.coor*x.square.sep
  y.coor <- y.coor*y.square.sep
  for(j in 3:ncol(table)){
    for(k in 1:nrow(table)){
      index <- as.numeric(table[k,j]);
      if(is.na(index) == F){
        if(j == ncol(table)){
          color <- colorRange.XLCj[round(index*gradientRate/max.XLCj, 0)+1];
          points(x=x.coor[j-2],y=y.coor[k],pch=pch,col=color, cex=square.cex);
        }else{
          color <- colorRange[round(index*gradientRate/max.XLC, 0)+1];
          points(x=x.coor[j-2],y=y.coor[k],pch=pch,col=color, cex=square.cex);
        }
      }
    }
  }
  if(length(reportInd) > 0){
    report.matched.rows <- grep(pattern=paste("^",reportInd, "$",collapse="|", sep=""), table[,1]);
    points(x=rep(ncol(table), length(report.matched.rows))*x.square.sep, y=((nrow(table)-report.matched.rows)+1)*y.square.sep, pch="*", cex=square.cex);
  }
  if(spBar){
    if(length(sp.color) < length(sp)){
      stop("If spBar = TRUE, then provide a vector of colors with length equal to species names\n");
    }
    sp <- as.character(unique(table$species));
    for(i in 1:length(sp)){
      rows <- nrow(table) - grep(pattern=sp[i], table[,2]) + 1;
      polygon(x=c(-1,-1,-0.7,-0.7), y=c(min(rows)*y.square.sep, max(rows)*y.square.sep, max(rows)*y.square.sep, min(rows)*y.square.sep), col=sp.color[i], border=F)
      text(y=mean(rows)*y.square.sep, x = min.xlim+x.spBar, label=sp[i], xpd = F, cex = 1, srt = 90);
    }
  }
  # scale bar XLCij
  x.coor <- (1:length(colorRange)*scaleBar.adj)+x.scaleBar;
  y.coor <- y.scaleBar;
  points(x=x.coor,y=rep(y.coor, length(colorRange)),pch=15,col=colorRange, cex=square.cex+0.5);
  text(x= x.coor[length(colorRange)/2], y=y.coor-XLCij.label, label="XLCij");
  text(x= x.coor[1]-min.XLC.pos, y=y.coor, label=as.character(0))
  text(x= x.coor[length(x.coor)]+max.XLC.pos, y=y.coor, label=as.character(round(max.XLC, 1)));
  # scale bar XLCj
  #x.coor <- (1:length(colorRange)*0.2)+3 #(1:length(colorRange)*0.1)+(max(x.coor)+5)
  y.coor <- nrow(table)+3;
  points(x=x.coor,y=rep(y.coor, length(colorRange)),pch=15,col=colorRange.XLCj, cex=square.cex+0.5);
  text(x= x.coor[length(colorRange)/2], y=y.coor+XLCj.label, label="XLCj");
  text(x= x.coor[1]-min.XLC.pos, y=y.coor, label=as.character(0));
  text(x= x.coor[length(x.coor)]+max.XLC.pos, y=y.coor, label=as.character(round(max.XLCj, 1)))
  mtext(side = 2, text = "Alleles", line = y.label.line, cex = square.cex+0.25, outer = T);
  mtext(side = 1, text = "Gene trees", line = x.label.line, cex = square.cex+0.25, outer = T);
  if(savePlot){
    dev.copy(pdf,plotName, height=height, width=width); #height=5.83, width=8.27 for A5 size
    dev.off();
  }
}

### write PhyloImap -- input Imap can either be a data.frame or name of an external Imap file
writePhyloImap <- function(Imap, species = NULL, msTaxaNames = F, taxaVec = NULL, PhyloImapName="PhyloImap.txt"){
  PhyloImap <- character();
  if(is.character(Imap)){
    Imap  <- read.table(Imap, sep="\t", header=T);
  }
  if(is.null(species)){
    species <- as.character(unique(Imap$species));
  }
  # writting Imap required by Phylonet
  if(msTaxaNames){
    if(length(taxaVec) == 0){
      stop("To write PhyloImap using ms taxa names, you need to provide a taxaVec\n");
    }
    if(is.integer(taxaVec) == FALSE){
      taxaVec <- as.integer(unlist(strsplit(taxaVec, split=" ")));
    }
    count <- 1;
    start <- 1;
    for(x in taxaVec){
      PhyloImap[count] <- paste(species[count], ":", paste(start:(start-1+x), collapse=","), ";", sep="");
      start <- start + x;
      count <- count +1;
    }
  }else{
    for(z in 1:length(species)){
      taxa <- Imap[Imap$species == species[z],];
      taxa <- taxa$traits;
      taxa <- paste(taxa, collapse=",");
      PhyloImap[z] <- paste(species[z], ":", taxa, ";", sep="");
    }
  }
  write(PhyloImap[1:length(species)], PhyloImapName);
}

sptree.vs.genetrees <- function(wd= getwd(), PhylonetDir = NULL, PhyloImapName="PhyloImap.txt",
                                tree.newick, genetreeNames, species, Imap, msTaxaNames = F, taxaVec = NULL, sp.treeName, sptree.dir = getwd(), geneTreeTempName="genetree.temp"){
  setwd(wd);
  if(is.null(PhylonetDir)){
    if("phylonet_v2_3.jar" %in% list.files()){
      PhylonetDir <- "phylonet_v2_3.jar";
    }else{
      stop("Phylonet executable not found, add phylonet_v2_3.jar to working directory or provide a path\nNote that sptree.vs.genetrees() might not work with a different version.\n
           To download Phylonet v2 visit: http://old-bioinfo.cs.rice.edu/phylonet/index.html\n");
    }
  }
  writePhyloImap(Imap=Imap, species=species, msTaxaNames = msTaxaNames, taxaVec = taxaVec, PhyloImapName = PhyloImapName);
  distVec <- integer();
  genetrees <- read.tree(genetreeNames);
  if(class(genetrees) == "multiPhylo"){
    cat("Estimating expected distribution of extra lineages, be patient...\n");
    count <- 1;
    for(l in 1:length(genetrees)){
      write.tree(genetrees[[l]], geneTreeTempName);
      command <- paste("java -jar ", PhylonetDir, " deep_coal_count ", sptree.dir, "/", sp.treeName," ", wd, "/", geneTreeTempName, " -a ", wd, "/", PhyloImapName, sep="");
      deepcoal <- system(command, intern = T);
      deepcoal <- as.integer(gsub(pattern=".*lineages: (\\d+)", replacement="\\1", deepcoal[2]));
      distVec <- c(distVec, deepcoal);
      cat("progress: ",l , " / ", length(genetrees), "\n");
      count <- count +1;
      file.remove(geneTreeTempName);
    }
  }else{
    command <- paste("java -jar ", PhylonetDir, " deep_coal_count ", wd, "/", sp.treeName," ", wd, "/", genetreeNames, " -a ", wd, "/", PhyloImapName, sep="");
    deepcoal <- system(command, intern = T);
    deepcoal <- as.integer(gsub(pattern=".*lineages: (\\d+)", replacement="\\1", deepcoal[2]));
    distVec <- c(distVec, deepcoal);
  }
  file.remove(PhyloImapName);
  return(distVec)
}