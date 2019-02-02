mask.data = function(tmpout, popDF, nSNP, mask) {

	#Rename populations if dealing with simulated data...
	strat_id = c()
	for(popu in 1:length(popDF$id)) {
		if(popu == 1) {
			strat_id = rep(as.character(popDF$id[popu]), popDF$sample.size[popu])
		} else {
			strat_id = append(strat_id, rep(as.character(popDF$id[popu]), popDF$sample.size[popu]))
		}
	}

	afreqs = alleleFreqs(tmpout)
	freqmat = sapply(afreqs, FUN=function(x){a=x[3]})
	variable = which(!is.na(freqmat) & freqmat > 0 & freqmat < 1)

	if(length(variable) <= nSNP) {
		return(print(paste("try again -- too few SNPs:", length(variable))))
	} else {
		print(paste(length(variable), "SNPs initially"))
	}

	#prelimgt = tmpout[,variable]
	
	tmpmat = as.matrix(tmpout)
	ordervec = sapply(strsplit(tmpmat[,2], " "), FUN=function(x) as.numeric(x[2]))
	tmpmat = tmpmat[order(ordervec),]
	tmpmat[,2] = strat_id
	tmpgt = as.data.frame(tmpmat)
	tmpgt = tmpgt[order(tmpgt$strata),]
	tmpgt$ids = paste0(tmpgt$strata, "_", c(1:length(tmpmat[,1])))
	
	#print("df2gtypes 1/2")
	#tmpgt = df2gtypes(tmpgt, ploidy = 2)

	onemat = as.matrix(tmpgt)
	
#	#a = tmpgt[,1:(2*nSNP+2)]
#	#masked = as.matrix(a)
#	#masked[mask == TRUE] = NA

	all1 = seq(3,dim(onemat)[2],2)
	starters = sample(all1, nSNP, replace = FALSE)   #
	a = tmpgt[,c(1:2,sort(c(starters, starters+1)))]    #

	masked = as.matrix(a)
	masked[mask == TRUE] = NA
	
	checknum = 1
	#print(paste0("checking for polymorphism #", checknum))
	onecoldata = matrix(data = NA, nrow = 2*dim(masked)[1], ncol = nSNP)
	for(x in 1:nSNP) {
		onecoldata[,x] = c(masked[,all1[x]], masked[,all1[x]+1])
	}
	fixme = apply(onecoldata, 2, FUN=function(x){length(unique(x[is.na(x)==FALSE]))})
	fixme = which(fixme == 1)

	#print(paste0("Need to fix ", length(fixme), " loci"))

	
#	#inuse = sort(c(all1[1:nSNP], all1[1:nSNP]+1))
	inuse = sort(c(starters, starters+1))   #

	if(length(fixme) > 0) {
		for(fl in 1:length(fixme)) {
			remaining = onemat
			remaining[,inuse] = NA
		
			mat.col.leaving = all1[fixme[fl]]  #I'm not sure about this?
			currmask = mask[,all1[fixme[fl]]]  #or this...

			remaining[currmask==TRUE,] = NA

			checknum = checknum+1
			#print(paste0("checking for polymorphism #", checknum))
		
			onecoldata = matrix(data = NA, nrow = 2*dim(remaining)[1], ncol = length(all1))
			for(x in 1:length(all1)) {
				onecoldata[,x] = c(remaining[,all1[x]], remaining[,all1[x]+1])
			}

			choices = apply(onecoldata, 2, FUN=function(x){length(unique(x[is.na(x)==FALSE]))})
			if(length(choices[choices == 2]) > 0) {
				if(length(choices[choices == 2] > 1)) {
					chose = sample(which(choices==2), 1, replace = FALSE)
				} else {
					chose = choices[choices == 2]
				}
			} else {
				return("BUMMER! too few SNPs due to missing data")
				print("BUMMER! too few SNPs due to missing data")
			}

			masked[,mat.col.leaving] = remaining[,all1[chose]]
			masked[,mat.col.leaving+1] = remaining[,all1[chose]+1]

			inuse = c(inuse, all1[chose], all1[chose]+1)
		}
	}

	postmask = as.data.frame(masked)
	#print("df2gtypes 2/2")
	postmask = df2gtypes(postmask, ploidy = 2)
	postmask
}

nomask = function(tmpout, popDF, nSNP) {
	#Rename populations if dealing with simulated data...
	strat_id = c()
	for(popu in 1:length(popDF$id)) {
		if(popu == 1) {
			strat_id = rep(as.character(popDF$id[popu]), popDF$sample.size[popu])
		} else {
			strat_id = append(strat_id, rep(as.character(popDF$id[popu]), popDF$sample.size[popu]))
		}
	}

	afreqs = alleleFreqs(tmpout)
	freqmat = sapply(afreqs, FUN=function(x){a=x[3]})
	variable = which(!is.na(freqmat) & freqmat > 0 & freqmat < 1)

	if(length(variable) < nSNP) {
		return(print(paste("try again -- too few SNPs:", length(variable))))
	} else {
		print(paste(length(variable), "SNPs initially"))
	}

	tmpmat = as.matrix(tmpout)
	ordervec = sapply(strsplit(tmpmat[,2], " "), FUN=function(x) as.numeric(x[2]))
	tmpmat = tmpmat[order(ordervec),]
	tmpmat[,2] = strat_id
	tmpgt = as.data.frame(tmpmat)
	tmpgt = tmpgt[order(tmpgt$strata),]
	tmpgt$ids = paste0(tmpgt$strata, "_", c(1:length(tmpmat[,1])))

	onemat = as.matrix(tmpgt)
	df = as.data.frame(onemat)
	gtype = df2gtypes(df, ploidy = 2)
	gtype
}

