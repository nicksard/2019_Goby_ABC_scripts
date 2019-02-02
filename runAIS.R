#Running AIS simulations in R via strataG
#Be sure to install strataG from github
#install.packages("devtools")
#library(devtools)
#install_github("EricArcher/strataG")
#Also install KScorrect 

#INFINITE SITES MODEL!
#NO MISSING DATA!

library(strataG)
library(swfscMisc)

args = commandArgs(TRUE)
model = as.character(args[1])  
node = as.numeric(args[2])
nrep = as.numeric(args[3])
#Usage -- Rscript runAIS.R Natural 1 100
#in bash: Rscript runAIS.R '${m[$i]}' ${PBS_ARRAYID} 100


codedir = "/mnt/research/ABC/RG_Flint/SRC"
simdir = system("echo $TMPDIR", intern = TRUE)
outdir = "/mnt/research/ABC/RG_Flint/OUT"

### DELETE ME ###
#codedir = "~/Google Drive/MSU/AquaticInvasives/FINAL/DEPLOY/NoMissing"
#simdir = "~/Desktop/AIS_SIMS"
#outdir = "~/Desktop/AIS_SIMS"
#node = 1
#nrep = 10
#r=1
#model = "Local"
### END DELETE ###

setwd(codedir)
source("AIS_strataG_fxns.R")
source("fast_stats_v2.R")
source("fast_mask_v2.R")   

bounds = read.csv("FlintPriors.csv", header = TRUE, row.names = 1)
popDF = read.csv("Flint_popDF.csv",  header = TRUE)

nSNP = 2921

#Run simulations for each replicate
for(r in 1:nrep) {
	setwd(simdir)

	timeneeds = c()
	maskout = c()
	if(exists("repout")) {
		rm(repout)
	}
		
	parmdraw = draw.parms(bounds = bounds, nrep = 1)
	fscmod = build.model(parmvec = parmdraw, model = model, npop = 9, nSNP = 100000, SNP.mu = 1e-7, seq.length = 10, Tnomig=0)
	attr(fscmod[[3]], "ploidy") = 2
	attr(fscmod[[3]], "opts") = paste("-I -s", nSNP)

	start = Sys.time()
	ntries = 1

	while(!exists("repout")) {
		if(ntries <= 5) {
			if(is.gtypes(maskout) == FALSE) {
				tmpout = myfsc(pop.info = fscmod[[1]], locus.params = fscmod[[3]], hist.ev = fscmod[[2]], mig.rates = fscmod[[4]], label = paste0(model, "_",node, "-", r), delete.files = TRUE, exec = "fsc25", label.haplotypes = FALSE)
				maskout = nomask(tmpout = tmpout, popDF = popDF, nSNP = nSNP)
				rm(tmpout)
				ntries = ntries + 1
				if(is.gtypes(maskout) == FALSE) {
					attr(fscmod[[3]], "num.chrom") = attr(fscmod[[3]], "num.chrom")*2
				}
			} else {
				print(paste("Calculating summary statistics for", model, "replicate", r))
				stats = getStats2(maskout)	  
				end = Sys.time()
				timeneeds = end-start
				parms = as.vector(parmdraw)
				names(parms) <- colnames(parmdraw)
				repout = c(mod=model, rep = r, parms, stats, time = timeneeds)
				rm(parms, parmdraw, stats, timeneeds, start, end)

				setwd(outdir)
				if(!file.exists(paste0(model,"_",node, ".dat"))) {
					write(names(repout), paste0(model,"_",node, ".dat"), ncolumns = length(repout))
					write(repout, paste0(model,"_",node, ".dat"), ncolumns = length(repout), append = TRUE)
				} else {
					write(repout, paste0(model,"_",node, ".dat"), ncolumns = length(repout), append = TRUE)
				}
			}
		} else {
			print(paste("BUMMER!  Not enough SNPs after", ntries-1,"attempts, moving to next replicate"))
			repout = c()
			rm(start)
			ntries = c()
		}
	}
}
