#R code used for running CMR analysis portion of the simulations in

#Benjamin, Griffin, and Douglas, "A nonparametric technique for analysis of state-trace functions

#This code uses data generated from runsims.m, in the folder simulated data. It requires the 
#Java Source and STACMR-R folders to run.

#3 parameters, shape, sep, and nhalf determine which dataset will be used for the analysis.
#For the analysis itself, nsamples remained at 10000 except for the 50 point case, which used 1000.



library("R.matlab")
library(xlsx)

datfolder = 'simulated data/'
shape = 'cdown'   #Shape can be linear, cdown, or cup
sep = 'large'     #sep can be small or large
nhalf = 4         #nhalf can be 4,10 or 50.
npoints = nhalf*2
readname = c('allsim_', sep, as.character(nhalf), '_', shape, '.mat')
readname = paste0(c(datfolder, readname), collapse = '')
rawdat = readMat(readname)
fulldat = rawdat$fulldata
fulldat = lapply(fulldat, '[[', 1)
dat = as.data.frame(fulldat[[120]])
rm(rawdat) #clear up workspace.

dirbase = getwd()
dirsim = paste0(dirbase, '/STACMR-R')
setwd(dirsim)
source("staCMRsetup.R") 

starttime = Sys.time()
writefile = 1
nsamples = 10000 #use 1000 if nhalf = 50.


iterations = 10 
noisevals = c(.1, .2, .4) #These don't get used aside from labels
interactvals = c(0, .1, .15, .2)

cmrstats = data.frame(matrix(vector(), 0, 4,
                             dimnames=list(c(), c("pval", "datafit", "avgfits", "stdfit"))),
                      stringsAsFactors=F)

cmrsummary = data.frame(matrix(vector(), 12, 4,
                             dimnames=list(c(), c("pval", "datafit", "avgfits", "stdfit"))),
                      stringsAsFactors=F)
allcmrstats = rep(list(cmrsummary), 12)


for (h in 1:length(noisevals)){
  for (i in 1:length(interactvals)){
    drow = i+(h-1)*length(interactvals); 
    
    for (j in 1:iterations){
      dcol = 12*(j-1)
      elem = drow+dcol
      dat = fulldat[[elem]]
      
      
      #Format data 
      nsubs = dim(dat)[1]
      datnew = data.frame(matrix(vector(), nsubs*4,3+nhalf))
      
      crow = 1
      for (k in 1:nsubs){
        for (depend in 1:2){
          for (att in 1:2){
            ccols = (depend-1)*(npoints) + (att-1)*nhalf + 1
            ccols = ccols:(ccols+nhalf-1)
            
            datnew[crow,1:3] = c(k, depend, att)
            datnew[crow,4:dim(datnew)[2]] = dat[k,ccols]
            
            crow = crow+1
          }
        }
      }
      
      
      if (nhalf > 10){
        output = staCMRFIT(datnew, nsample=nsamples, approx = 1)
      } else {
        output = staCMRFIT(datnew, nsample=nsamples)
      }
      
      cmrstats[j,"pval"] = output$p
      cmrstats[j,"datafit"] = output$datafit
      cmrstats[j,"avgfits"] = mean(output$fits)
      cmrstats[j,"stdfit"] = sd(output$fits)
    }
    allcmrstats[[drow]] = cmrstats
    
    for (n in 1:4){
      cmrsummary[drow,n] = mean(cmrstats[,n])
    }
    midtime = Sys.time()
    print(midtime - starttime)
  }
}


cmrall = data.frame(matrix(nrow = 12, ncol = 4*iterations+1))
for (m in 1:4){
  colstart = (m-1)*iterations+1
  colfinish = m*iterations
  for (n in 1:12){
    cmrall[n,colstart:colfinish] = allcmrstats[[n]][,m]
  }
}
cmrall = cbind(cmrall, cmrsummary)
clabel = c(rep('pval', 10), rep('datafit', 10), rep('avgfits', 10), rep('stdfit', 10), '', 
           'mean p-val', 'mean datafit', 'mean avgfits', 'mean stdfit')
colnames(cmrall) = clabel


writename = paste0('Rsimulations_cmr_', shape, sep, '.xlsx')
writenameall = paste0('Rsimulations_cmr_full', shape, sep, '.xlsx')

sheetname = 0
if (sep == 'small'){
  sheetname = sheetname+3
}
if (nhalf == 4) {
  sheetname = sheetname+1
}else if (nhalf == 10){
  sheetname = sheetname+2
}else {
  sheetname = sheetname+3
}
sheetname = paste0('Sheet', as.character(sheetname))


setwd(dirbase)
if (writefile == 1){
  write.xlsx(cmrsummary, writename, sheetName = sheetname, row.names=FALSE, append=TRUE)
  write.xlsx(cmrall, writenameall, sheetName = sheetname, row.names=FALSE, append=TRUE)
}
