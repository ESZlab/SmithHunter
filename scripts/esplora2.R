#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("Correct syntax is: Rscript esplora.R /path/COV1.txt /path/COV2.txt /path/folderconsoglie /path/fastacandidati.fas")
}
fileprima <- args[1]
filedopo <- args[2]
foldersoglie <- args[3]
fastacandidati <- args[4]
# Following from the main script: 
# if the nuclear mapping WAS performed COV1 is overall coverage, COV2 is coverage of uniquely mitochondrial reads
# if the nuclear mapping was NOT performed, COV2 is a copy of COV1


# reads in thresholds from files T1 e T2
T1 <- as.numeric(readLines(paste(foldersoglie, '/T1', sep = ""))) # global threshold as number
T2files <- list.files(path = foldersoglie, pattern = 'T2', full.names = TRUE) # files with replicate-specific thresholds, namesorted
library(gtools)
T2files <- mixedsort(T2files) # numerical sort
T2 <- vector() # will take T2 of replicates
for (T2file in T2files) {
    T2 <- append(T2, as.numeric(readLines(T2file)))
} # T2 of replicates as one vector


# output file names
plot_nuc_mt <- gsub('.txt$' ,'.pdf', filedopo) # coverage plot (including global and uniquely mitochondrial reads, if nuclear mapping is performed)
plot_replicates <- gsub('.txt$' ,'.replicates.pdf', filedopo) # replicates plot
plot_candidates <- gsub('.txt$' ,'.clusters.pdf', filedopo) # clusters plot
thresholds <- gsub('.txt$' ,'.stats', filedopo) # stats as text


# READS IN COVERAGE DATA
# reads prima, coverage after remapping of all reads
# reads dopo, coverage after remapping of uniquely mitochondrial reads (copy of prima if nuclear remapping was NOT performed)
# columns are: base, coverage in each replicate, total coverage
prima <- read.csv(file=fileprima, header=F, sep='\t')
dopo <- read.csv(file=filedopo, header=F, sep='\t')
n_replicates <- ncol(prima)-2

# gives informative column names to both dataframes
names <- colnames(prima)
names[1] <- "base"
for (r in 1:n_replicates) {
    names[r+1] <- paste('R', r, sep='')
}
names[n_replicates+2] <- 'tot'
colnames(prima) <- names
colnames(dopo) <- names



# PLOT COVERAGE BEFORE AND (if available) AFTER NUCLEAR REMAPPING 
pdf(plot_nuc_mt, onefile=T, paper='A4', width = 21/2.54, height = 29.7/2.54)
par(mfrow=c(2,1), mar=c(4, 4.1, 1, 1), mgp=c(2.5,1,0))

# coverage from prima and dopo, global threshold is indicated
# difference from red and green are peaks that were removed as of nuclear origin (i.e. reads can remap on nuclear geome)
plot(prima$base, prima$tot, type='l', xlab='position over genome', ylab='coverage (all replicates)', col='red')
points(dopo$base, dopo$tot, type='l', col='darkgreen')
abline(h=T1, col='blue')

# coverage from dopo, global threshold is indicated, detail of low coverage area
# note: threshold is here indicated over coverage, but will (in fact) be applied to cluster depth
plot(dopo$base, dopo$tot, type='l', xlab='position over genome', ylab='coverage (all replicates, detail)', col='darkgreen', ylim=c(0, 5*T1))
abline(h=T1, col='blue')

trash <- dev.off()



# PLOT COVERAGE OF INDIVIDUAL REPLICATES
# coverage from individual replicates, replicate-specific thresholds are indicated, detail of low coverage area
# note: threshold is here indicated over coverage, but will (in fact) be applied to cluster depth
pdf(plot_replicates, onefile=T, paper='A4', width = 21/2.54, height = 29.7/2.54)
par(mfrow=c(n_replicates,1), mar=c(4, 4.1, 1, 1), mgp=c(2.5,1,0))

# plot cycle
for (r in 1:n_replicates) {
    rname <- paste('R', r, sep='')
    T2_replica <- T2[r] # threshold for this specific replicates, from T2
    plot(dopo$base, dopo[,r+1], type='l', xlab='position over genome', ylab=paste('coverage (', rname, '), detail', sep=''), col='darkgreen', ylim=c(0, 5*T2_replica))
    abline(h=T2_replica, col='blue')
}
trash <- dev.off()



# STATS TEXTUAL OUTPUT
# note: stats are different depending of whether nuclear remapping has been performed or not
if (!identical(prima, dopo)){ # prima and dopo differ, there has been a nuclear remapping, more info is produced 
    
    # thresholds
    cat(paste('Global_threshold: ',T1, '\n', sep=''),file=thresholds)
    cat(paste('Replicates_thresholds: ', paste(T2, collapse=" "), '\n', sep=''),file=thresholds,append=TRUE)
    # total coverage
    cat(paste('All_replicates_meancoverage_first_remapping: ', round(mean(prima$tot)), '\n', sep=''),file=thresholds,append=TRUE)
    cat(paste('All_replicates_meancoverage_second_remapping: ', round(mean(dopo$tot)), '\n', sep=''),file=thresholds,append=TRUE)
    cat(paste('Coverage_decreased_to: ', round((mean(dopo$tot)/mean(prima$tot))*100), '%\n', sep=''),file=thresholds,append=TRUE)
    # replicate coverage
    for (r in 1:n_replicates) {
        rname <- paste('R', r, sep='')
        cat(paste(rname, '_meancoverage_first_remapping: ', round(mean(prima[,r+1])), '\n', sep=''),file=thresholds,append=TRUE)
        cat(paste(rname, '_meancoverage_second_remapping: ', round(mean(dopo[,r+1])), '\n', sep=''),file=thresholds,append=TRUE)
    }
    # compare replicates
    cat(paste('Relative_coverage: '),file=thresholds,append=TRUE)
    for (r in 1:n_replicates) {
        cat(paste(round((mean(dopo[,r+1])/mean(dopo$tot))*100), '% ', sep=''),file=thresholds,append=TRUE)
    }
    cat('\n',file=thresholds,append=TRUE)
    
} else { # prima and dopo are equal, there has NOT been a nuclear remapping, less info is produced. Based on 'dopo'.
    
    # thresholds
    cat(paste('Global_threshold: ',T1, '\n', sep=''),file=thresholds)
    cat(paste('Replicates_thresholds: ', paste(T2, collapse=" "), '\n', sep=''),file=thresholds,append=TRUE)
    # total coverage
    cat(paste('All_replicates_meancoverage: ', round(mean(dopo$tot)), '\n', sep=''),file=thresholds,append=TRUE)
    # replicate coverage
    for (r in 1:n_replicates) {
        rname <- paste('R', r, sep='')
        cat(paste(rname, '_meancoverage: ', round(mean(dopo[,r+1])), '\n', sep=''),file=thresholds,append=TRUE)
    }
    # compare replicates
    cat(paste('Relative_coverage: '),file=thresholds,append=TRUE)
    for (r in 1:n_replicates) {
        cat(paste(round((mean(dopo[,r+1])/mean(dopo$tot))*100), '% ', sep=''),file=thresholds,append=TRUE)
    }
    cat('\n',file=thresholds,append=TRUE)
    
}



#PLOT CANDIDATES (CLUSTERS)
# read candidates from fasta
fasta <- readLines(fastacandidati)
headers <- fasta[grepl("^>", fasta)]

# headers to dataframe
candidates = data.frame(from=numeric(0),to=numeric(0),strand=character(0),depth=numeric(0),clusterid=numeric(0))
for (c in 1:length(headers)) { # for each candidate
        candidate <- unlist(strsplit(gsub('>clusterid','',headers[c]), "_size|_pos|_|_strand")) # from, to, strand, depth, clusterid
        candidater <- list(as.numeric(candidate[3]), as.numeric(candidate[4]), candidate[5] , as.numeric(candidate[2]), as.numeric(candidate[1])) # ordine, tipo
        candidates <- rbind(candidates, candidater) # adds candidates to dataframe
}
colnames(candidates) <- c('from', 'to', 'strand', 'depth', 'clusterid')
candidates <-candidates[order(candidates$from),] # order based on start

# finds appropriate vertical space between cluster and its label based on y scale
vs <- ((max(candidates[,4])*1.1)/30)

# adds cluster mids (as column 6) and depths (as column 7) to help reposition labels avoiding overlaps
candidates$center <- (candidates[,1]+candidates[,2])/2
candidates$ylabel <- candidates[,4]+vs # right over cluster, but there may be overalp with nearby clusters

# moves overlapping labels upwards, in 5 successive approximations
for (rounds in 1:5){ # successive approximations
    for (r in 2:nrow(candidates)){
        for (rp in 1:(r-1)){ 
            # if one of the previous clusters is too close on both x and y (i.e. there would be overlap)...
            if (((candidates[r,6] - candidates[rp,6]) < 400) & (abs(candidates[r,7] - candidates[rp,7]) < vs)){
                candidates[r, 7] = candidates[r, 7]+vs
            } # moves upwards the labels of clusters nearby
        }
    }
}

# appropriate ylim for the plot
mp <- (max(candidates[,7])*1.2) # highest label +20%

# splits the genome in 4 segments to foster readability
genomelen <- dopo[nrow(dopo),1]
primo <- c(-50, floor(genomelen/4)+50)
secondo <- c(ceiling(genomelen/4)-50, floor(genomelen/2)+50)
terzo <- c(ceiling(genomelen/2)-50, floor(genomelen*3/4)+50)
quarto <- c(ceiling(genomelen*3/4)-50, genomelen+50)
campi <- list(primo, secondo, terzo, quarto)

# the plot
pdf(plot_candidates, onefile=T, paper='A4', width = 21/2.54, height = 29.7/2.54)

par(mfrow=c(4,1), mar=c(4, 4.1, 0.5, 1), mgp=c(2.5,1,0), xaxs="i", yaxs="i") # xaxs='i' rende i margini precisi
for (campo in campi) { # for each of the four plotting fields
    # plot coverage and global threshold
    plot(dopo$base, dopo$tot, type='l', xlim=campo, ylim=c(0, mp), xlab='position over genome', ylab='cluster depth', col='grey')
        abline(h=T1, col='blue')
    # and then plot clusters as rectangles
    for (c in 1:nrow(candidates)) { # for each candidate
        if (candidates[c,3] == '+'){ # in red if on plus strand (genome as submitted)
            rect(candidates[c,1], 0, candidates[c,2], candidates[c,4], border='red')
        }
        if (candidates[c,3] == '-'){ # in green if on minus strand
            rect(candidates[c,1], 0, candidates[c,2], candidates[c,4], border='green')
        }
    }
    text(candidates[,6], candidates[,7], candidates[,5], cex=0.7) # plot labels
}

trash <- dev.off()
