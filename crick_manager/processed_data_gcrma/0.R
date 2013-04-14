scatter.with.hist <-
    function(x, y, xlab = 'x', ylab = 'y') {
    xy <- data.frame(x,y)  
    scatter <- qplot(x,y, data=xy, xlab = xlab, ylab = ylab)  + 
         scale_x_continuous(limits=c(min(x),max(x))) + 
         scale_y_continuous(limits=c(min(y),max(y))) + 
         geom_rug(col=rgb(.5,0,0,alpha=.2)) +
         opts(title = paste(xlab, "VS", ylab))
    scatter

}

source('/home/baldig/projects/proteomics/Huang/maratio/R/bayesreg.R')

fname = "0.csv"
numC = 7
numE = 7
conf = 1
winSize = 101
bayesIntC = T
bayesIntE = T
doMulttest = T
ppde = F

raw.data <- read.table(fname, sep = '\t', header = T)
d <- data.frame(raw.data[,c(2:dim(raw.data)[2])])
#here we start at row 2, assuming only one probe id later used for conversion
dim(d)

#combined_ans <- bayesT(d, numC, numE, winSize = winSize, conf = conf, bayesIntC=T, bayesIntE=T, doMulttest=T, ppde=F)

ans <- bayesT(d, numC, numE, winSize = winSize, conf = conf, bayesIntC=bayesIntC, bayesIntE=bayesIntE, doMulttest=doMulttest, ppde=ppde)

ans$diffMean <- ans$meanE - ans$meanC

control <- data.frame(cbind(ans$meanC, ans$stdC, ans$bayesSDC))
experiment <- data.frame(cbind(ans$meanE, ans$stdE, ans$bayesSDE))

control <- subset(control, X1 > 0)
experiment <- subset(experiment, X1 > 0)


## MEAN
mc <- ans$meanC
minmeanC <- min(mc[mc > 0]) / 1000.0
mc[mc==0] <- minmeanC
ratioMean <- ans$meanE / mc
ans <- cbind(ans, ratioMean)


#proteinname <- raw.data[,1]
#accession <- raw.data[,2]
#wd <- cbind(proteinname, accession, ans[,c('meanC', 'meanE', 'stdC', 'stdE', 'bayesSDC', 'bayesSDE', 'pVal', 'diffMean', 'ratioMean')])
write.table(ans, paste0(fname,"cybert.out"), sep='\t', row.names=F, quote=F)


if(F){
library(ggplot2)

all.pval <- data.frame(ans)
above.pval <- data.frame(subset(ans, pVal > 0.05)) #, select=c('meanE', 'bayesSDE', 'stdE', 'diffMean')))
below.pval <- data.frame(subset(ans, pVal <= 0.05)) #, select=c('meanE', 'bayesSDE', 'stdE', 'diffMean')))

all.pval$th <- paste('ALL-MEANS', dim(all.pval)[1])
above.pval$th <- paste('FALSE', dim(above.pval)[1])
below.pval$th <- paste('TRUE', dim(below.pval)[1])
all.data <- rbind(all.pval, above.pval, below.pval)

control <- data.frame(cbind(ans$meanC, ans$stdC, ans$bayesSDC))
experiment <- data.frame(cbind(ans$meanE, ans$stdE, ans$bayesSDE))
control <- subset(control, X1 > 0)
experiment <- subset(experiment, X1 > 0)
pdf(paste('Mean_vs_SD_', fname, '.pdf', sep=''), 8, 6)
ggplot(all.data, aes(log(meanE), fill = th)) + geom_density(alpha = 0.2) + opts(title = expression("MEAN OF THE EXPERIMENT"))
ggplot(all.data, aes(log(diffMean), fill = th)) + geom_density(alpha = 0.2) + opts(title = expression("DIFF MEAN OF THE EXPERIMENT"))
ggplot(all.data, aes(log(ratioMean), fill = th)) + geom_density(alpha = 0.2) + opts(title = expression("RATIO MEAN OF THE EXPERIMENT"))
ggplot(all.data, aes(log(bayesSDE), fill = th)) + geom_density(alpha = 0.2) + opts(title = expression("BAYES SD OF THE EXPERIMENT"))
ggplot(all.data, aes(log(stdE), fill = th)) + geom_density(alpha = 0.2) + opts(title = expression("SD OF THE EXPERIMENT"))
scatter.with.hist(log(experiment[,1]), log(experiment[,2]), xlab='Experiment Mean', ylab='SD')
scatter.with.hist(log(experiment[,1]), log(experiment[,3]), xlab='Experiment Mean', ylab='POOLED-bayesSD')
scatter.with.hist(log(experiment[,2]), log(experiment[,3]), xlab='Experiment SD', ylab='POOLED-bayesSD')
scatter.with.hist(log(experiment[,3]), log(experiment2[,3]), xlab='POOLED-bayesSD', ylab='AVG-bayesSD')
dev.off()
}
