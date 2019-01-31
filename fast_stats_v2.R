getStats2 = function(out) {

	freqs = alleleFreqs(out, TRUE)
	popid = strataNames(out)
	split_out = strataSplit(out)

	freqmat = matrix(data = NA, nrow = length(freqs), ncol = length(popid))
	for(p in 1:length(freqmat[1,])) {
		freqmat[,p] = sapply(freqs, FUN=function(x){a=x[,,popid[p]][3]})
	}
	colnames(freqmat) = popid

	variable_loci = matrix(data = 0, nrow = length(freqs), ncol = length(popid))
	for(p in 1:length(variable_loci[1,])) {
		variable_loci[freqmat[,p] > 0 & freqmat[,p] < 1 & !is.na(freqmat[,p]),p] = 1
	}
	colnames(variable_loci) = popid

	nind_mat = matrix(data = NA, nrow = length(freqs), ncol = length(popid))
	for(p in 1:length(nind_mat[1,])) {
		nind_mat[,p] = sapply(freqs, FUN=function(x){a=x[,,popid[p]][1]+x[,,popid[p]][2]})
	}
	colnames(nind_mat) = popid

	#Local SNP number
	localSNP = colSums(variable_loci)
	names(localSNP) = paste0("S.", popid)

	Smean = mean(localSNP)
	names(Smean) = "mean_S"

	Ssd = sd(localSNP)
	names(Ssd) = "sd_S"

	#Private SNP number
	privateSNP = c()
	for(p in 1:length(variable_loci[1,])) {
		privateSNP[p] = length(which(rowSums(variable_loci) == 1 & variable_loci[,p] == 1))
	}
	names(privateSNP) = paste0("pS.", popid)

	prSmean = mean(privateSNP)
	names(prSmean) = "mean_prS"

	prSsd = sd(privateSNP)
	names(prSsd) = "sd_prS"

	total_priv = sum(privateSNP)
	names(total_priv) = "tot_priv"

	#Pairwise private SNPs and their frequencies
	combos = combn(popid, 2)
	pwPriv = c()
	pwPrivFREQ = c()
	for(pair in 1:length(combos[1,])) {
		subset = which(nind_mat[,combos[1,pair]] > 0 & nind_mat[,combos[2,pair]] > 0)
		tmpfreq = freqmat[subset,c(combos[1,pair], combos[2,pair])]
		tmpvar = variable_loci[subset,c(combos[1,pair], combos[2,pair])]

		tmpfreq2 = tmpfreq[tmpvar[,1] == 1 & tmpvar[,2] == 0,]
		tmpfreq3 = tmpfreq[tmpvar[,1] == 0 & tmpvar[,2] == 1,]

		if(!class(tmpfreq2) == "matrix") {
			tmpfreq2 = matrix(data = tmpfreq2, ncol = 2, nrow = 1, byrow = TRUE)
		}
		
		if(!class(tmpfreq3) == "matrix") {
			tmpfreq3 = matrix(data = tmpfreq3, ncol = 2, nrow = 1, byrow = TRUE)
		}

		tmpfreq2[tmpfreq2[,2] == 1,] = 1-tmpfreq2[tmpfreq2[,2] == 1,]
		tmpfreq3[tmpfreq3[,1] == 1,] = 1-tmpfreq3[tmpfreq3[,1] == 1,]

		tmp = c(dim(tmpfreq2)[1], dim(tmpfreq3)[1])
		#tmp = c(length(which(tmpvar[,1] == 1 & tmpvar[,2] == 0)), length(which(tmpvar[,2] == 1 & tmpvar[,1] == 0)))
		names(tmp) = c(paste0("pwPriv.",combos[1,pair],".",combos[2,pair]), paste0("pwPriv.",combos[2,pair],".",combos[1,pair]))
		pwPriv = c(pwPriv, tmp)

		tmp2 = c(mean(tmpfreq2[,1]), mean(tmpfreq3[,2]))
		names(tmp2) = c(paste0("pwPrivFREQ.",combos[1,pair],".",combos[2,pair]), paste0("pwPrivFREQ.",combos[2,pair],".",combos[1,pair]))
		pwPrivFREQ = c(pwPrivFREQ,tmp2)
		
		rm(tmpvar, tmpfreq, tmpfreq2, tmpfreq3, tmp, tmp2, subset)
	}

	#He
	Hetmp = sapply(split_out, FUN=function(x){exptdHet(x)})
	Hetmp[nind_mat == 0] = NA
	He = colMeans(Hetmp, na.rm = TRUE)
	names(He) = paste0("He.",  colnames(Hetmp))

	HeALL = mean(exptdHet(out))
	names(HeALL) = "HeALL"

	#Simple Fst calculation... FAST  ~ 1 second
	pairFst = c()
	for(pair in 1:length(combos[1,])) {
		subset = which(nind_mat[,combos[1,pair]] > 0 & nind_mat[,combos[2,pair]] > 0)
		Hs_tmp = (nind_mat[subset,combos[1,pair]]*Hetmp[subset,combos[1,pair]] + nind_mat[subset,combos[2,pair]]*Hetmp[subset,combos[2,pair]])/(nind_mat[subset,combos[1,pair]]+nind_mat[subset,combos[2,pair]])
		combo_freq = (nind_mat[subset,combos[1,pair]]*freqmat[subset,combos[1,pair]] + nind_mat[subset,combos[2,pair]]*freqmat[subset,combos[2,pair]])/(nind_mat[subset,combos[1,pair]]+nind_mat[subset,combos[2,pair]])
		Ht_tmp = 1-combo_freq^2-(1-combo_freq)^2
		Fst_tmp = (Ht_tmp-Hs_tmp)/Ht_tmp
		pairFst[pair] = mean(Fst_tmp, na.rm = TRUE)
		names(pairFst)[pair] = paste0(combos[1,pair],".",combos[2,pair])
		rm(subset, Hs_tmp, combo_freq, Ht_tmp, Fst_tmp)
	}

	#Fst calculation in hierfstat -- VERY SLOOOW ~ 15 min
	#genind_out = gtypes2genind(out)
	#FstMat = pairwise.fst(genind_out, pop = genind_out$pop)
	#pwFst2 = as.vector(FstMat)
	#names(pwFst2) = paste0(combos[1,],".", combos[2,])

	#Fst calculation using popStructTest from strataG -- NOT FAST ~ 4 min
#	fst = popStructTest(out, nrep = 0, type = "both", stats = "fst", quietly = TRUE)
#	pairFst = fst$pairwise$pair.mat$Fst[lower.tri(fst$pairwise$pair.mat$Fst) == TRUE]
#	compared = rownames(fst$pairwise$pair.mat$Fst)
#	Fstnames = c()
#	for(pid in 1:(length(compared)-1)) {
#		Fstnames = c(Fstnames, paste0(compared[pid],".",compared[(pid+1):length(compared)]))
#	}
#	names(pairFst) = Fstnames
#	tot_Fst = fst$overall$result[,1]
	tot_Fst = statFst(out)$result["estimate"]
	names(tot_Fst) = "tot_Fst"

	stats = c(localSNP, Smean, Ssd, privateSNP, prSmean, prSsd, total_priv, pwPriv, pwPrivFREQ, He, HeALL, pairFst, tot_Fst)
	stats
}

