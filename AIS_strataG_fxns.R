#Functions for AIS model simulations using strataG instead of fsc25 directly

#draw.parms() 
#function to draw parameter values for the analysis -- will want to write this to a CSV 
#so it can be used for both the Flint and Au Sable systems
draw.parms = function(bounds=NULL, nrep = NULL) {
	#Make sure all columns are present
	if(dim(bounds)[2] != 5) {
		stop("bounds should be a matrix with 4 columns (min, max, distribution, is integer?, hyper) and # parameter rows")
	}

	#Establish a matrix for the parameter values
	parmmat = matrix(data = NA, nrow = nrep, ncol = length(bounds[,1]))
	colnames(parmmat) = rownames(bounds)

	#Group sorted parameters are sampled first
	sorted = which(bounds$dist == "groupsort")
	nsorted = length(sorted)
	min=as.numeric(as.character(bounds$min[sorted[1]]))
	max=as.numeric(as.character(bounds$max[sorted[1]]))
	for(r in 1:nrep) {
		#parmmat[r,sorted] = sort(runif(nsorted, min, max))
		parmmat[r,sorted] = sort(sample(c(min:max),nsorted,replace=FALSE))
	} 
	rm(min)
	rm(max)

	#Loop over the number of replicates to draw
	for(p in 1:length(bounds[,1])) {
		if(bounds$dist[p] == "uniform") {
			min=as.numeric(as.character(bounds$min[p]))
			max=as.numeric(as.character(bounds$max[p]))
			parmmat[,p] = runif(nrep, min, max)
			rm(min)
			rm(max)
		} else if(bounds$dist[p] == "loguniform") {
			require(KScorrect)
			min=as.numeric(as.character(bounds$min[p]))
			max=as.numeric(as.character(bounds$max[p]))
			parmmat[,p] = rlunif(nrep, min, max, base=10)
			rm(min)
			rm(max)
		} else if(bounds$dist[p] == "conditional") {
			for(r in 1:nrep) {
				if(length(which(rownames(bounds)==bounds$min[p])) == 0) {
					min=as.numeric(as.character(bounds$min[p]))
					max=parmmat[r,which(rownames(bounds)==bounds$max[p])]
					parmmat[r,p] = runif(1, min, max)
					rm(min)
					rm(max)
				} else if(length(which(rownames(bounds)==bounds$max[p])) == 0) {
					min=parmmat[r,which(rownames(bounds)==bounds$max[p])]
					max=as.numeric(as.character(bounds$min[p]))
					parmmat[r,p] = runif(1, min, max)
					rm(min)
					rm(max)
				}
			}
		} 
	}

	#Hyper parameters drawn here...
	n.hyper = length(which(bounds$dist == "hyper"))
	mean_col = which(bounds$hyper == "mean")
	cv_col = which(bounds$hyper == "cv")

	for(repl in 1:nrep) {
		hypsd = parmmat[repl,mean_col]*parmmat[repl,cv_col]
		draw = rnorm(10*n.hyper, mean = parmmat[repl,mean_col], sd = hypsd)
		draw = draw[draw > 0]
		parmmat[repl,which(bounds$dist == "hyper")] = sample(draw, n.hyper, replace = FALSE)
		rm(draw, hypsd)
	}



	for(p in 1:length(parmmat[1,])) {
		#Round the drawn values if the parameter is an integer
		if(bounds$isint[p] == 1) {
			parmmat[,p] = round(parmmat[,p], 0)
		}
	}

	#Name columns and output matrix
	parmmat
}


#build.model()
#function to build a model using strataG
#This takes one line from parmmat... keep that in mind
build.model = function(parmvec = NULL, model = NULL, npop = NULL, nSNP = NULL, SNP.mu = NULL, seq.length=NULL, Tnomig=0) {
	require(strataG)
	fscmodel = vector("list", 4)
	#parmvec = as.data.frame(t(parmvec))
	parmvec = as.data.frame(parmvec)

	#1 fscPopInfo()
	npop = npop
	pop.size = c(parmvec$NATIVE, parmvec$STCLAIR, parmvec$ERIE, parmvec$MICH, parmvec$SAGBAY, parmvec$SAGRIV, parmvec$FLRIV, parmvec$MOTT, parmvec$HOLL)
	sample.size = c(0,12,20,0,13,0,20,20,20) #Native, St. Clair, Erie, Michigan, Saginaw Bay, Saginaw River, below Mott, Mott, Holloway
	fscmodel[[1]] = fscPopInfo(pop.size = pop.size, sample.size = sample.size, sample.times = 0, growth.rate = 0)

	#2 fscHistEv()
	fscmodel[[2]] = getEvents(parmvec = parmvec, model = model)
	
	#3 fscLocusParams()
    #createLocusParams <- function(chr, type, num.markers, recomb.rate, 
    #    param.4, param.5, param.6, ploidy, num.chrom) {
    #    if (is.null(chr)) 
    #        chr <- 1
    #    df <- data.frame(chromosome = chr, type = type, num.markers = num.markers, 
    #        recomb.rate = recomb.rate, param.4 = param.4, param.5 = param.5, 
    #        param.6 = param.6, stringsAsFactors = FALSE)
    #    df <- df[order(df$chromosome), ]
    #    attr(df, "num.chrom") <- num.chrom[1]
    #    attr(df, "ploidy") <- ploidy
    #    return(df)
    #}

	#fscmodel[[3]] = createLocusParams(1,"DNA",seq.length,0,SNP.mu,0.33,NA,2,nSNP) 
	#attr(fscmodel[[3]], "opts") = "-I"
	fscmodel[[3]] = fscLocusParams(locus.type = "dna", sequence.length = seq.length, num.chrom = nSNP, mut.rate = SNP.mu)
	
	#fscmodel[[3]] = fscLocusParams(locus.type = "snp", num.loci = nSNP, mut.rate = seq.length*SNP.mu)
	#4 Migration Matrix - elsewhere [getMig() fxn]
	fscmodel = addMig(fscmod = fscmodel, parmvec = parmvec, model = model, Tnomig=Tnomig)
	
	fscmodel

}


addMig = function(fscmod = fscmod, parmvec = parmdraw, model = "Local", Tnomig=0) {
	zerome = c()
	fscmod[[2]] = fscmod[[2]][order(as.numeric(fscmod[[2]][,1]), as.numeric(fscmod[[2]][,4])),]

	eventlist = fscmod[[2]]

	migmat = vector("list", 1)
	migmat[[1]]=matrix(data = 0, nrow = 9, ncol = 9)
	migmat[[1]][2:5,2:5] = parmvec$SHP
	migmat[[1]][5,6] = parmvec$MBR
	migmat[[1]][6,5] = parmvec$MBR
	migmat[[1]][6,7] = parmvec$MDS 
	migmat[[1]][7,8] = parmvec$MDS
	migmat[[1]][8,9] = parmvec$MDS

	if(Tnomig == 0) {
		if(model == "Local") {
			migmat[[1]][7,6] = parmvec$MUS
			migmat[[1]][8,7] = parmvec$MUS
			migmat[[1]][9,8] = parmvec$MUS
		} else if(model == "MixH") {
			migmat[[1]][9,2:5] = parmvec$BBM
		} else if(model == "MixM") {
			migmat[[1]][8,2:5] = parmvec$BBM
		} else if(model == "MixHM") {
			migmat[[1]][8:9,2:5] = parmvec$BBM
		} else if(model == "ErieH") {
			migmat[[1]][9,3] = parmvec$BBM
		} else if(model == "ErieM") {
			migmat[[1]][8,3] = parmvec$BBM
		} else if(model == "ErieHM") {
			migmat[[1]][8:9,3] = parmvec$BBM
		} else if(model == "MichH") {
			migmat[[1]][9,4] = parmvec$BBM
		} else if(model == "MichM") {
			migmat[[1]][8,4] = parmvec$BBM
		} else if(model == "MichHM") {
			migmat[[1]][8:9,4] = parmvec$BBM
		} else if(model == "StCH") {
			migmat[[1]][9,2] = parmvec$BBM
		} else if(model == "StCM") {
			migmat[[1]][8,2] = parmvec$BBM
		} else if(model == "StCHM") {
			migmat[[1]][8:9,2] = parmvec$BBM
		} else if(model == "SagBayH") {
			migmat[[1]][9,5] = parmvec$BBM
		} else if(model == "SagBayM") {
			migmat[[1]][8,5] = parmvec$BBM
		} else if(model == "SagBayHM") {
			migmat[[1]][8:9,5] = parmvec$BBM
		}
	}

	diag(migmat[[1]]) = 0

	m_number = 0
	have_mig_mat = FALSE

	for(y in 1:length(eventlist[,1])) {
		if(eventlist[y,1] < Tnomig) {
			if(eventlist[y,2] == eventlist[y,3]) {
				fscmod[[2]][y,7] = m_number
			} else if(eventlist[y,2] != eventlist[y,3]) {
				m_number = m_number+1

				fscmod[[2]][y,7] = m_number

				migmat[[m_number+1]] = migmat[[m_number]]
				zerome = append(zerome, as.numeric(eventlist[y,2])+1)
				migmat[[m_number+1]][zerome,] = 0
				migmat[[m_number+1]][,zerome] = 0
			}
		} else if(eventlist[y,1] >= Tnomig) {
			if(eventlist[y,1] == min(eventlist[,1][eventlist[,1] >= Tnomig])) {
				if(eventlist[y,2] == eventlist[y,3]) {
					if(have_mig_mat == FALSE) {
						m_number = m_number+1

						fscmod[[2]][y,7] = m_number

						migmat[[m_number+1]] = migmat[[m_number]]

						if(Tnomig > 0) {
							if(model == "Local") {
								migmat[[m_number+1]][7,6] = parmvec$MUS
								migmat[[m_number+1]][8,7] = parmvec$MUS
								migmat[[m_number+1]][9,8] = parmvec$MUS
							} else if(model == "MixH") {
								migmat[[m_number+1]][9,2:5] = parmvec$BBM
							} else if(model == "MixM") {
								migmat[[m_number+1]][8,2:5] = parmvec$BBM
							} else if(model == "MixHM") {
								migmat[[m_number+1]][8:9,2:5] = parmvec$BBM
							} else if(model == "ErieH") {
								migmat[[m_number+1]][9,3] = parmvec$BBM
							} else if(model == "ErieM") {
								migmat[[m_number+1]][8,3] = parmvec$BBM
							} else if(model == "ErieHM") {
								migmat[[m_number+1]][8:9,3] = parmvec$BBM
							} else if(model == "MichH") {
								migmat[[m_number+1]][9,4] = parmvec$BBM
							} else if(model == "MichM") {
								migmat[[m_number+1]][8,4] = parmvec$BBM
							} else if(model == "MichHM") {
								migmat[[m_number+1]][8:9,4] = parmvec$BBM
							} else if(model == "StCH") {
								migmat[[m_number+1]][9,2] = parmvec$BBM
							} else if(model == "StCM") {
								migmat[[m_number+1]][8,2] = parmvec$BBM
							} else if(model == "StCHM") {
								migmat[[m_number+1]][8:9,2] = parmvec$BBM
							} else if(model == "SagBayH") {
								migmat[[m_number+1]][9,5] = parmvec$BBM
							} else if(model == "SagBayM") {
								migmat[[m_number+1]][8,5] = parmvec$BBM
							} else if(model == "SagBayHM") {
								migmat[[m_number+1]][8:9,5] = parmvec$BBM
							}
						}

						have_mig_mat = TRUE
						
						migmat[[m_number+1]][zerome,] = 0
						migmat[[m_number+1]][,zerome] = 0

						diag(migmat[[m_number+1]]) = 0
					} else {
						fscmod[[2]][y,7] = m_number
					}

				} else if(eventlist[y,2] != eventlist[y,3]) {
					m_number = m_number+1

					fscmod[[2]][y,7] = m_number

					migmat[[m_number+1]] = migmat[[m_number]]

					if(model == "Local") {
						migmat[[m_number+1]][7,6] = parmvec$MUS
						migmat[[m_number+1]][8,7] = parmvec$MUS
						migmat[[m_number+1]][9,8] = parmvec$MUS
					} else if(model == "MixH") {
						migmat[[m_number+1]][9,2:5] = parmvec$BBM
					} else if(model == "MixM") {
						migmat[[m_number+1]][8,2:5] = parmvec$BBM
					} else if(model == "MixHM") {
						migmat[[m_number+1]][8:9,2:5] = parmvec$BBM
					} else if(model == "ErieH") {
						migmat[[m_number+1]][9,3] = parmvec$BBM
					} else if(model == "ErieM") {
						migmat[[m_number+1]][8,3] = parmvec$BBM
					} else if(model == "ErieHM") {
						migmat[[m_number+1]][8:9,3] = parmvec$BBM
					} else if(model == "MichH") {
						migmat[[m_number+1]][9,4] = parmvec$BBM
					} else if(model == "MichM") {
						migmat[[m_number+1]][8,4] = parmvec$BBM
					} else if(model == "MichHM") {
						migmat[[m_number+1]][8:9,4] = parmvec$BBM
					} else if(model == "StCH") {
						migmat[[m_number+1]][9,2] = parmvec$BBM
					} else if(model == "StCM") {
						migmat[[m_number+1]][8,2] = parmvec$BBM
					} else if(model == "StCHM") {
						migmat[[m_number+1]][8:9,2] = parmvec$BBM
					} else if(model == "SagBayH") {
						migmat[[m_number+1]][9,5] = parmvec$BBM
					} else if(model == "SagBayM") {
						migmat[[m_number+1]][8,5] = parmvec$BBM
					} else if(model == "SagBayHM") {
						migmat[[m_number+1]][8:9,5] = parmvec$BBM
					}		

					zerome = append(zerome, as.numeric(eventlist[y,2])+1)
					migmat[[m_number+1]][zerome,] = 0
					migmat[[m_number+1]][,zerome] = 0

					diag(migmat[[m_number+1]]) = 0
				}
			} else {
				if(eventlist[y,2] == eventlist[y,3]) {
					fscmod[[2]][y,7] = m_number
				} else if(eventlist[y,2] != eventlist[y,3]) {
					m_number = m_number+1

					fscmod[[2]][y,7] = m_number

					migmat[[m_number+1]] = migmat[[m_number]]

					zerome = append(zerome, as.numeric(eventlist[y,2])+1)
					migmat[[m_number+1]][zerome,] = 0
					migmat[[m_number+1]][,zerome] = 0

					diag(migmat[[m_number+1]]) = 0
				}
			}
		}
	}

	matsums = sapply(migmat, sum)
	migs = which(matsums > 0)
	newmat = vector("list", length(migs))
	for(m in 1:length(migs)) {
		newmat[[m]] = migmat[[m]]
	}

	nomigs = which(matsums == 0)
	fscmod[[2]][fscmod[[2]][,7] %in% (nomigs-1),7] = "nomig"

	fscmod[[4]] = newmat

	fscmod

}	


#getEvents()
#Function to build the HistEv list for the particular model... complicated!!
getEvents = function(parmvec = NULL, model = NULL) {

	Tnomig = 6

	if(model == "Local") {
		div.times=as.numeric(parmvec[1,17:23])
		div.times = c(div.times[1:6], div.times[6:7])
		bott.times = div.times-1
		event.times = sort(c(bott.times, div.times))
		events = matrix(data = NA, nrow = length(event.times), ncol = 7)
		events[,1] = event.times #time
		events[,2] = c(8,8,7,7,6,6,5,5,4,4,2,2,3,3,1,1) #source (from)
		events[,3] = c(8,7,7,6,6,5,5,4,4,1,2,1,3,1,1,0) #sink (to)
		events[,4] = rep(c(0,1),8) #proportion migrants
		events[,5] = c(rep(c(parmvec$RF,1),4), rep(c(parmvec$SF,1),3), c(parmvec$IF,1)) #new size
		events[,6] = 0 #new growth
		events[,7] = c(0,1,1,2,2,3,3,4,4,5,5,"nomig",5,rep("nomig",3)) #migration matrix
		events[events[,2] == "2" & events[,3] == "1",1] = as.numeric(events[events[,2] == "2" & events[,3] == "1",1])+1
		events[events[,2] == "3" & events[,3] == "3",1] = as.numeric(events[events[,2] == "3" & events[,3] == "3",1])-1

		#events = rbind(events, c(6, 0, 0, 0, 1, 0, 6))

	} else if(model == "SagBayH") {
		div.times=as.numeric(parmvec[1,17:23])  
		div.times = c(div.times[1:6], div.times[6:7])
		bott.times = div.times-1
		event.times = sort(c(bott.times, div.times))
		events = matrix(data = NA, nrow = length(event.times), ncol = 7)
		events[,1] = event.times #time
		events[,2] = c(6,6,7,7,8,8,5,5,4,4,2,2,3,3,1,1) #source (from)
		events[,3] = c(6,7,7,8,8,4,5,4,4,1,2,1,3,1,1,0) #sink (to)
		events[,4] = rep(c(0,1),8) #proportion migrants
		events[,5] = c(rep(c(parmvec$RF,1),4), rep(c(parmvec$SF,1),3), c(parmvec$IF,1)) #new size
		events[,6] = 0 #new growth
		events[,7] = c(0,1,1,2,2,3,3,4,4,5,5,"nomig",5,rep("nomig",3))
		events[events[,2] == "2" & events[,3] == "1",1] = as.numeric(events[events[,2] == "2" & events[,3] == "1",1])+1
		events[events[,2] == "3" & events[,3] == "3",1] = as.numeric(events[events[,2] == "3" & events[,3] == "3",1])-1

		#events = rbind(events, c(6, 0, 0, 0, 1, 0, 6))

	} else if(model == "SagBayHM") {
		div.times=as.numeric(parmvec[1,17:23])  
		div.times = c(div.times[1:6], div.times[6:7])
		bott.times = div.times-1
		event.times = sort(c(bott.times, div.times))
		events = matrix(data = NA, nrow = length(event.times), ncol = 7)
		events[,1] = event.times #time
		events[,2] = c(6,6,8,8,7,7,5,5,4,4,2,2,3,3,1,1) #source (from)
		events[,3] = c(6,7,8,4,7,4,5,4,4,1,2,1,3,1,1,0) #sink (to)
		events[,4] = rep(c(0,1),8) #proportion migrants
		events[,5] = c(rep(c(parmvec$RF,1),4), rep(c(parmvec$SF,1),3), c(parmvec$IF,1)) #new size
		events[,6] = 0 #new growth
		events[,7] = c(0,1,1,2,2,3,3,4,4,5,5,"nomig",5,rep("nomig",3))
		events[events[,2] == "2" & events[,3] == "1",1] = as.numeric(events[events[,2] == "2" & events[,3] == "1",1])+1
		events[events[,2] == "3" & events[,3] == "3",1] = as.numeric(events[events[,2] == "3" & events[,3] == "3",1])-1

		#events = rbind(events, c(6, 0, 0, 0, 1, 0, 6))

	} else if(model == "SagBayM") {
		div.times=as.numeric(parmvec[1,17:23])  
		div.times = c(div.times[1:6], div.times[6:7])
		bott.times = div.times-1
		event.times = sort(c(bott.times, div.times))
		events = matrix(data = NA, nrow = length(event.times), ncol = 7)
		events[,1] = event.times #time
		events[,2] = c(8,8,6,6,7,7,5,5,4,4,2,2,3,3,1,1) #source (from)
		events[,3] = c(8,7,6,7,7,4,5,4,4,1,2,1,3,1,1,0) #sink (to)
		events[,4] = rep(c(0,1),8) #proportion migrants
		events[,5] = c(rep(c(parmvec$RF,1),4), rep(c(parmvec$SF,1),3), c(parmvec$IF,1)) #new size
		events[,6] = 0 #new growth
		events[,7] = c(0,1,1,2,2,3,3,4,4,5,5,"nomig",5,rep("nomig",3))
		events[events[,2] == "2" & events[,3] == "1",1] = as.numeric(events[events[,2] == "2" & events[,3] == "1",1])+1
		events[events[,2] == "3" & events[,3] == "3",1] = as.numeric(events[events[,2] == "3" & events[,3] == "3",1])-1

		#events = rbind(events, c(6, 0, 0, 0, 1, 0, 6))

	} else if(model == "ErieH") {
		div.times=as.numeric(parmvec[1,17:23])  
		div.times = c(div.times[1:6], div.times[6:7])
		bott.times = div.times-1
		event.times = sort(c(bott.times, div.times))
		events = matrix(data = NA, nrow = length(event.times), ncol = 7)
		events[,1] = event.times #time
		events[,2] = c(6,6,7,7,8,8,5,5,4,4,2,2,3,3,1,1) #source (from)
		events[,3] = c(6,7,7,8,8,2,5,4,4,1,2,1,3,1,1,0) #sink (to)
		events[,4] = rep(c(0,1),8) #proportion migrants
		events[,5] = c(rep(c(parmvec$RF,1),4), rep(c(parmvec$SF,1),3), c(parmvec$IF,1)) #new size
		events[,6] = 0 #new growth
		events[,7] = c(0,1,1,2,2,3,3,4,4,5,5,"nomig",5,rep("nomig",3))
		events[events[,2] == "2" & events[,3] == "1",1] = as.numeric(events[events[,2] == "2" & events[,3] == "1",1])+1
		events[events[,2] == "3" & events[,3] == "3",1] = as.numeric(events[events[,2] == "3" & events[,3] == "3",1])-1

		#events = rbind(events, c(6, 0, 0, 0, 1, 0, 6))


	} else if(model == "ErieHM") {
		div.times=as.numeric(parmvec[1,17:23])  
		div.times = c(div.times[1:6], div.times[6:7])
		bott.times = div.times-1
		event.times = sort(c(bott.times, div.times))
		events = matrix(data = NA, nrow = length(event.times), ncol = 7)
		events[,1] = event.times #time
		events[,2] = c(6,6,8,8,7,7,5,5,4,4,2,2,3,3,1,1) #source (from)
		events[,3] = c(6,7,8,2,7,2,5,4,4,1,2,1,3,1,1,0) #sink (to)
		events[,4] = rep(c(0,1),8) #proportion migrants
		events[,5] = c(rep(c(parmvec$RF,1),4), rep(c(parmvec$SF,1),3), c(parmvec$IF,1)) #new size
		events[,6] = 0 #new growth
		events[,7] = c(0,1,1,2,2,3,3,4,4,5,5,"nomig",5,rep("nomig",3))
		events[events[,2] == "2" & events[,3] == "1",1] = as.numeric(events[events[,2] == "2" & events[,3] == "1",1])+1
		events[events[,2] == "3" & events[,3] == "3",1] = as.numeric(events[events[,2] == "3" & events[,3] == "3",1])-1

		#events = rbind(events, c(6, 0, 0, 0, 1, 0, 6))

	} else if(model == "ErieM") {
		div.times=as.numeric(parmvec[1,17:23])  
		div.times = c(div.times[1:6], div.times[6:7])
		bott.times = div.times-1
		event.times = sort(c(bott.times, div.times))
		events = matrix(data = NA, nrow = length(event.times), ncol = 7)
		events[,1] = event.times #time
		events[,2] = c(8,8,6,6,7,7,5,5,4,4,2,2,3,3,1,1) #source (from)
		events[,3] = c(8,7,6,7,7,2,5,4,4,1,2,1,3,1,1,0) #sink (to)
		events[,4] = rep(c(0,1),8) #proportion migrants
		events[,5] = c(rep(c(parmvec$RF,1),4), rep(c(parmvec$SF,1),3), c(parmvec$IF,1)) #new size
		events[,6] = 0 #new growth
		events[,7] = c(0,1,1,2,2,3,3,4,4,5,5,"nomig",5,rep("nomig",3))
		events[events[,2] == "2" & events[,3] == "1",1] = as.numeric(events[events[,2] == "2" & events[,3] == "1",1])+1
		events[events[,2] == "3" & events[,3] == "3",1] = as.numeric(events[events[,2] == "3" & events[,3] == "3",1])-1

		#events = rbind(events, c(6, 0, 0, 0, 1, 0, 6))
	
	} else if(model == "MichH") {
		div.times=as.numeric(parmvec[1,17:23])  
		div.times = c(div.times[1:6], div.times[6:7])
		bott.times = div.times-1
		event.times = sort(c(bott.times, div.times))
		events = matrix(data = NA, nrow = length(event.times), ncol = 7)
		events[,1] = event.times #time
		events[,2] = c(6,6,7,7,8,8,5,5,4,4,2,2,3,3,1,1) #source (from)
		events[,3] = c(6,7,7,8,8,3,5,4,4,1,2,1,3,1,1,0) #sink (to)
		events[,4] = rep(c(0,1),8) #proportion migrants
		events[,5] = c(rep(c(parmvec$RF,1),4), rep(c(parmvec$SF,1),3), c(parmvec$IF,1)) #new size
		events[,6] = 0 #new growth
		events[,7] = c(0,1,1,2,2,3,3,4,4,5,5,"nomig",5,rep("nomig",3))
		events[events[,2] == "2" & events[,3] == "1",1] = as.numeric(events[events[,2] == "2" & events[,3] == "1",1])+1
		events[events[,2] == "3" & events[,3] == "3",1] = as.numeric(events[events[,2] == "3" & events[,3] == "3",1])-1

		#events = rbind(events, c(6, 0, 0, 0, 1, 0, 6))


	} else if(model == "MichHM") {
		div.times=as.numeric(parmvec[1,17:23])  
		div.times = c(div.times[1:6], div.times[6:7])
		bott.times = div.times-1
		event.times = sort(c(bott.times, div.times))
		events = matrix(data = NA, nrow = length(event.times), ncol = 7)
		events[,1] = event.times #time
		events[,2] = c(6,6,8,8,7,7,5,5,4,4,2,2,3,3,1,1) #source (from)
		events[,3] = c(6,7,8,3,7,3,5,4,4,1,2,1,3,1,1,0) #sink (to)
		events[,4] = rep(c(0,1),8) #proportion migrants
		events[,5] = c(rep(c(parmvec$RF,1),4), rep(c(parmvec$SF,1),3), c(parmvec$IF,1)) #new size
		events[,6] = 0 #new growth
		events[,7] = c(0,1,1,2,2,3,3,4,4,5,5,"nomig",5,rep("nomig",3))
		events[events[,2] == "2" & events[,3] == "1",1] = as.numeric(events[events[,2] == "2" & events[,3] == "1",1])+1
		events[events[,2] == "3" & events[,3] == "3",1] = as.numeric(events[events[,2] == "3" & events[,3] == "3",1])-1

		#events = rbind(events, c(6, 0, 0, 0, 1, 0, 6))

	} else if(model == "MichM") {
		div.times=as.numeric(parmvec[1,17:23])  
		div.times = c(div.times[1:6], div.times[6:7])
		bott.times = div.times-1
		event.times = sort(c(bott.times, div.times))
		events = matrix(data = NA, nrow = length(event.times), ncol = 7)
		events[,1] = event.times #time
		events[,2] = c(8,8,6,6,7,7,5,5,4,4,2,2,3,3,1,1) #source (from)
		events[,3] = c(8,7,6,7,7,3,5,4,4,1,2,1,3,1,1,0) #sink (to)
		events[,4] = rep(c(0,1),8) #proportion migrants
		events[,5] = c(rep(c(parmvec$RF,1),4), rep(c(parmvec$SF,1),3), c(parmvec$IF,1)) #new size
		events[,6] = 0 #new growth
		events[,7] = c(0,1,1,2,2,3,3,4,4,5,5,"nomig",5,rep("nomig",3))
		events[events[,2] == "2" & events[,3] == "1",1] = as.numeric(events[events[,2] == "2" & events[,3] == "1",1])+1
		events[events[,2] == "3" & events[,3] == "3",1] = as.numeric(events[events[,2] == "3" & events[,3] == "3",1])-1

		#events = rbind(events, c(6, 0, 0, 0, 1, 0, 6))

	} else if(model == "StCH") {
		div.times=as.numeric(parmvec[1,17:23])  
		div.times = c(div.times[1:6], div.times[6:7])
		bott.times = div.times-1
		event.times = sort(c(bott.times, div.times))
		events = matrix(data = NA, nrow = length(event.times), ncol = 7)
		events[,1] = event.times #time
		events[,2] = c(6,6,7,7,8,8,5,5,4,4,2,2,3,3,1,1) #source (from)
		events[,3] = c(6,7,7,8,8,1,5,4,4,1,2,1,3,1,1,0) #sink (to)
		events[,4] = rep(c(0,1),8) #proportion migrants
		events[,5] = c(rep(c(parmvec$RF,1),4), rep(c(parmvec$SF,1),3), c(parmvec$IF,1)) #new size
		events[,6] = 0 #new growth
		events[,7] = c(0,1,1,2,2,3,3,4,4,5,5,"nomig",5,rep("nomig",3))
		events[events[,2] == "2" & events[,3] == "1",1] = as.numeric(events[events[,2] == "2" & events[,3] == "1",1])+1
		events[events[,2] == "3" & events[,3] == "3",1] = as.numeric(events[events[,2] == "3" & events[,3] == "3",1])-1

		#events = rbind(events, c(6, 0, 0, 0, 1, 0, 6))


	} else if(model == "StCHM") {
		div.times=as.numeric(parmvec[1,17:23])  
		div.times = c(div.times[1:6], div.times[6:7])
		bott.times = div.times-1
		event.times = sort(c(bott.times, div.times))
		events = matrix(data = NA, nrow = length(event.times), ncol = 7)
		events[,1] = event.times #time
		events[,2] = c(6,6,8,8,7,7,5,5,4,4,2,2,3,3,1,1) #source (from)
		events[,3] = c(6,7,8,1,7,1,5,4,4,1,2,1,3,1,1,0) #sink (to)
		events[,4] = rep(c(0,1),8) #proportion migrants
		events[,5] = c(rep(c(parmvec$RF,1),4), rep(c(parmvec$SF,1),3), c(parmvec$IF,1)) #new size
		events[,6] = 0 #new growth
		events[,7] = c(0,1,1,2,2,3,3,4,4,5,5,"nomig",5,rep("nomig",3))
		events[events[,2] == "2" & events[,3] == "1",1] = as.numeric(events[events[,2] == "2" & events[,3] == "1",1])+1
		events[events[,2] == "3" & events[,3] == "3",1] = as.numeric(events[events[,2] == "3" & events[,3] == "3",1])-1

		#events = rbind(events, c(6, 0, 0, 0, 1, 0, 6))

	} else if(model == "StCM") {
		div.times=as.numeric(parmvec[1,17:23])  
		div.times = c(div.times[1:6], div.times[6:7])
		bott.times = div.times-1
		event.times = sort(c(bott.times, div.times))
		events = matrix(data = NA, nrow = length(event.times), ncol = 7)
		events[,1] = event.times #time
		events[,2] = c(8,8,6,6,7,7,5,5,4,4,2,2,3,3,1,1) #source (from)
		events[,3] = c(8,7,6,7,7,1,5,4,4,1,2,1,3,1,1,0) #sink (to)
		events[,4] = rep(c(0,1),8) #proportion migrants
		events[,5] = c(rep(c(parmvec$RF,1),4), rep(c(parmvec$SF,1),3), c(parmvec$IF,1)) #new size
		events[,6] = 0 #new growth
		events[,7] = c(0,1,1,2,2,3,3,4,4,5,5,"nomig",5,rep("nomig",3))
		events[events[,2] == "2" & events[,3] == "1",1] = as.numeric(events[events[,2] == "2" & events[,3] == "1",1])+1
		events[events[,2] == "3" & events[,3] == "3",1] = as.numeric(events[events[,2] == "3" & events[,3] == "3",1])-1

		#events = rbind(events, c(6, 0, 0, 0, 1, 0, 6))

	} else if(model == "MixH") {
		div.times=as.numeric(parmvec[1,17:23])  
		div.times = c(div.times[1:6], div.times[6:7])
		bott.times = div.times-1
		div.times = c(div.times[1:2], rep(div.times[3], 4), div.times[4:8])
		event.times = sort(c(bott.times, div.times))
		events = matrix(data = NA, nrow = length(event.times), ncol = 7)
		events[,1] = event.times #time
		events[,2] = c(6,6,7,7,8,8,8,8,8,5,5,4,4,2,2,3,3,1,1) #source (from)
		events[,3] = c(6,7,7,8,8,2,3,1,4,5,4,4,1,2,1,3,1,1,0) #sink (to)
		events[,4] = c(rep(c(0,1),2), 0, 0.25, 0.33, 0.5, 1, rep(c(0,1),5)) #proportion migrants
		events[,5] = c(rep(c(parmvec$RF,1),3), rep(1,3), rep(c(parmvec$RF,1),1), rep(c(parmvec$SF,1),3), c(parmvec$IF,1)) #new size
		events[,6] = 0 #new growth
		events[,7] = c(0,1,1,2,2,3,3,3,3,3,4,4,5,5,"nomig",5,rep("nomig",3))
		events[events[,2] == "2" & events[,3] == "1",1] = as.numeric(events[events[,2] == "2" & events[,3] == "1",1])+1
		events[events[,2] == "3" & events[,3] == "3",1] = as.numeric(events[events[,2] == "3" & events[,3] == "3",1])-1

		#events = rbind(events, c(6, 0, 0, 0, 1, 0, 6))


	} else if(model == "MixHM") {
		div.times=as.numeric(parmvec[1,17:23])  
		div.times = c(div.times[1:6], div.times[6:7])
		bott.times = div.times-1
		div.times = c(div.times[1], rep(div.times[2], 4), rep(div.times[3], 4), div.times[4:8])
		event.times = sort(c(bott.times, div.times))
		events = matrix(data = NA, nrow = length(event.times), ncol = 7)
		events[,1] = event.times #time
		events[,2] = c(6,6,8,8,8,8,8,7,7,7,7,7,5,5,4,4,2,2,3,3,1,1) #source (from)
		events[,3] = c(6,7,8,2,3,1,4,7,2,3,1,4,5,4,4,1,2,1,3,1,1,0) #sink (to)
		events[,4] = c(rep(c(0,1),1), 0, 0.25, 0.33, 0.5, 1, 0, 0.25, 0.33, 0.5, 1, rep(c(0,1),5)) #proportion migrants
		events[,5] = c(rep(c(parmvec$RF,1),2), rep(1,3), rep(c(parmvec$RF,1),1), rep(1,3), rep(c(parmvec$RF,1),1), rep(c(parmvec$SF,1),3), c(parmvec$IF,1)) #new size
		events[,6] = 0 #new growth
		events[,7] = c(0,1,1,2,2,2,2,2,3,3,3,3,3,4,4,5,5,"nomig",5,rep("nomig",3))
		events[events[,2] == "2" & events[,3] == "1",1] = as.numeric(events[events[,2] == "2" & events[,3] == "1",1])+1
		events[events[,2] == "3" & events[,3] == "3",1] = as.numeric(events[events[,2] == "3" & events[,3] == "3",1])-1

		#events = rbind(events, c(6, 0, 0, 0, 1, 0, 6))

	} else if(model == "MixM") {
		div.times=as.numeric(parmvec[1,17:23])  
		div.times = c(div.times[1:6], div.times[6:7])
		bott.times = div.times-1
		div.times = c(div.times[1:2], rep(div.times[3], 4), div.times[4:8])
		event.times = sort(c(bott.times, div.times))
		events = matrix(data = NA, nrow = length(event.times), ncol = 7)
		events[,1] = event.times #time
		events[,2] = c(8,8,6,6,7,7,7,7,7,5,5,4,4,2,2,3,3,1,1) #source (from)
		events[,3] = c(8,7,6,7,7,2,3,1,4,5,4,4,1,2,1,3,1,1,0) #sink (to)
		events[,4] = c(rep(c(0,1),2), 0, 0.25, 0.33, 0.5, 1, rep(c(0,1),5)) #proportion migrants
		events[,5] = c(rep(c(parmvec$RF,1),3), rep(1,3), rep(c(parmvec$RF,1),1), rep(c(parmvec$SF,1),3), c(parmvec$IF,1)) #new size
		events[,6] = 0 #new growth
		events[,7] = c(0,1,1,2,2,3,3,3,3,3,4,4,5,5,"nomig",5,rep("nomig",3))
		events[events[,2] == "2" & events[,3] == "1",1] = as.numeric(events[events[,2] == "2" & events[,3] == "1",1])+1
		events[events[,2] == "3" & events[,3] == "3",1] = as.numeric(events[events[,2] == "3" & events[,3] == "3",1])-1

		#events = rbind(events, c(6, 0, 0, 0, 1, 0, 6))
	}

	events

}

#myfsc()
#Replacement fastsimcoal() fxn (turning off -S option)

myfsc = function (pop.info, locus.params, mig.rates = NULL, hist.ev = NULL, 
    label = NULL, quiet = TRUE, delete.files = TRUE, exec = "fsc252", 
    num.cores = NULL, label.haplotypes = TRUE) 
{
    if (is.null(label)) 
        label <- "fsc.run"
    label <- make.names(label)
    if (file.exists(label)) 
        for (f in dir(label, full.names = T)) file.remove(f)
    if (!quiet) 
        cat("fastsimcoal: writing input file\n")
    infile <- fscWrite(pop.info = pop.info, locus.params = locus.params, 
        mig.rates = mig.rates, hist.ev = hist.ev, label = label)
    if (!quiet) 
        cat("fastsimcoal: running\n")
    cores.spec <- if (!is.null(num.cores)) {
        num.cores <- max(1, num.cores)
        num.cores <- min(num.cores, min(detectCores(), 12))
        if (num.cores == 1) 
            ""
        else paste(c("-c", "-B"), num.cores, collapse = " ")
    }
    else ""
    cmd.line <- paste(exec, "-i", infile, "-n 1", ifelse(quiet, 
        "-q", ""), cores.spec, attr(locus.params, "opts"))
    err <- if (.Platform$OS.type == "unix") {
        system(cmd.line, intern = F)
    }
    else {
        shell(cmd.line, intern = F)
    }
    if (err == 0) {
        if (!quiet) 
            cat("fastsimcoal exited normally\n")
    }
    else {
        stop("fastsimcoal exited with error ", err, "\n")
    }
    arp.file <- file.path(label, paste(label, "_1_1.arp", sep = ""))
    if (!file.exists(arp.file)) 
        stop("fastsimcoal did not generate output")
    if (!quiet) 
        cat("fastsimcoal: parsing output to gtypes\n")
    g <- myfscRead(arp.file, locus.params, label.haplotypes)
    if (delete.files) {
        if (!quiet) 
            cat("fastsimcoal: removing output files\n")
        unlink(label, recursive = TRUE, force = TRUE)
        file.remove(infile)
        file.remove("seed.txt")
    }
    return(g)
}



myfscRead = function (file, locus.params, label.haplotypes = FALSE) 
{
    .formatGenotypes <- function(x, ploidy) {
    	require(swfscMisc)
        nloci <- ncol(x) - 2
        loc.end <- seq(ploidy, nrow(x), by = ploidy)
        gen.data <- do.call(rbind, lapply(loc.end, function(i) {
            allele.i <- (i - ploidy + 1):i
            loci <- as.vector(x[allele.i, -(1:2)])
            id <- paste(x[allele.i, 2], collapse = ".")
            pop <- x[allele.i[1], 1]
            c(id, pop, loci)
        }))
        locus_names <- paste("Locus", zero.pad(1:nloci), sep = "_")
        locus_names <- paste(rep(locus_names, each = ploidy), 
            1:ploidy, sep = ".")
        colnames(gen.data) <- c("id", "pop", locus_names)
        gen.data
    }

    .formatDNA <- function(dna.seq, pop, locus.params, label) {
        #arp.file <- file.path(label, paste(label, "_1_1.arp", sep = ""))
        #num.chrom <- attr(locus.params, "num.chrom")
        num.chrom = system(paste("grep 'Total number of polymorphic sites:'", file, "| cut -f 2 -d : | cut -f 2 -d ' '"), intern = TRUE)
        chrom.pos <- if (is.null(num.chrom)) {
            tapply(locus.params$num.markers, locus.params$chromosome, 
                sum)
        }
        else {
            rep(sum(1), num.chrom)
        }
        chrom.pos <- cumsum(chrom.pos)
        chrom.pos <- cbind(start = c(1, chrom.pos[-length(chrom.pos)] + 
            1), end = chrom.pos)
        rownames(dna.seq) <- pop
        dna.seq <- tolower(dna.seq)
        new("multidna", lapply(1:nrow(chrom.pos), function(i) {
            as.matrix(dna.seq)[, chrom.pos[i, "start"]:chrom.pos[i, 
                "end"]]
        }))
    }
    #print("reading file")
    f <- readLines(file)
    start <- grep("SampleData=", f) + 1
    end <- which(f == "}") - 2
    pos <- cbind(start, end)
    .compileMatrix <- function(i, pos) {
        f.line <- f[pos[i, 1]:pos[i, 2]]
        f.line <- gsub("[[:space:]]+", "--", f.line)
        result <- do.call(rbind, strsplit(f.line, "--"))[, -2]
        cbind(rep(paste("Sample", i), nrow(result)), result)
    }
    #print("compiling matrix")
    data.mat <- do.call(rbind, lapply(1:nrow(pos), .compileMatrix, 
        pos = pos))
    ploidy <- attr(locus.params, "ploidy")
    data.type <- f[grep("DataType=", f)]
    data.type <- gsub("\tDataType=", "", data.type)
    switch(data.type, DNA = {
    	#print("splitting strings")
        dna.seq <- do.call(rbind, strsplit(data.mat[, 3], ""))
        if (attr(locus.params, "ploidy") == 2) {
            gen.data <- .formatGenotypes(cbind(data.mat[, 1:2], 
                dna.seq), ploidy)
            df2gtypes(gen.data, ploidy, description = file)
        } else {
        	#print("formatting DNA")
            dna.seq <- .formatDNA(dna.seq, data.mat[, 2], locus.params)
            #print("seq2gtype - SOMETHING IS WRONG HERE!!!  Has to do with creation of gen.data")
            g <- sequence2gtypes(dna.seq, strata = data.mat[,1], description = file)
            if (label.haplotypes) labelHaplotypes(g)$gtype else g
        }
    }, MICROSAT = {
        gen.data <- .formatGenotypes(data.mat, ploidy)
        df2gtypes(gen.data, ploidy, description = file)
    }, NULL)
}

