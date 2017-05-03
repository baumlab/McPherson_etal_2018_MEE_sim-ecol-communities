simul.comms.2<-function (s=c(5, 10, 15, 20, 25, 30, 35, 40), r=10, tr.pool=traits, 
		presence="random", abun.method="lnorm:1",
abundance="random", commin=NULL, commax=NULL, comnot=NULL, comnot.type=NULL, dups=0) {

# increase memory limit (4095Mb for 32-bit R; values up to 8Tb possible 
# for 64-bit R on 64-bit Windows)
memory.limit(10000)

# Load required packages and functions
library(scales)
library(moments)

# function to fill list of abundance vectors based on user specifications
fill.abun<-function(empty, comsize, presprob, abunprob, abundist, abunsd, ids){
 nn<-length(empty)/2
set<-sample(1:length(empty), size=(comsize-ids), replace=FALSE, prob=presprob)
		if(ids>0){set<-c(set, nn+sample(set, size=dups, replace=FALSE, prob=NULL))}
		sorted<-set[order(abunprob[set])]
if(abundist == "lnorm"){counts<-0.5+rlnorm(n=length(set), meanlog=0, sdlog=abunsd)}
if(abundist == "norm"){counts<-0.5+rnorm(n=length(set), mean=0, sd=abunsd)}
if(abundist == "unif"){counts<-0.5+runif(n=length(set), min=0, max=abunsd)}
if(abundist == "fixed"){counts<-rep(abunsd, times=length(set))}
empty[sorted]<-sort(counts)
return(empty)
} # end function fill.abun
# function to ensure positive abundance values (strictly > 0)
pos.abun<-function(vec){
		vmin<-min(vec)
		pres<-which(vec!=0)
		pmin<-min(abs(vec[pres]))
		nvec<-numeric(length=length(vec))
		nvec[pres]<-vec[pres]+abs(vmin)+pmin
		return(nvec)
} # end function pos.abun

# Error checks
t<-ncol(tr.pool)-1
if(nrow(tr.pool)<max(s)){
stop(paste("The number of species (rows) in 'tr.pool' must be greater than the", "maximum species richness level in 's'.", sep=" "))
}
if(ncol(tr.pool)<2){
		stop(paste("Matrix 'tr.pool' must have at least two columns: one for a trait and one for", 
" trait likelihood", sep=" "))
}
 if(!is.null(commin) && length(commin)<t){
stop(paste("Vector commin must contain", t, "values.", sep=" "))
}
if(!is.null(commax) && length(commax)<t){
stop(paste("Vector commax must contain", t, "values.", sep=" "))
}
if(!is.null(comnot) && length(comnot)<t){
stop(paste("Option 'comnot' must be a list that contains", t, "elements.", sep=" "))
}
if(!is.null(comnot.type) && length(comnot.type)<t){
 stop(paste("Option 'comnot.type' must be a list that contains", t, "elements.", sep=" "))
}
if(dups>(0.5*min(s))){
stop("The value in dups cannot exceed 50% of the minimum value specified in s!")
}

# assignments
 l.s <- length(s)
r.rep <- rep(r, l.s)
nb.sp <- rep(s, r.rep)
if(is.null(presence)){presence<-"random"}
presence<-match.arg(presence, c("random", "common", "rare"))
if(is.null(abundance)){abundance<-"random"}
abundance<-match.arg(abundance, c("random", "common", "rare"))
if(is.null(commin)){commin<-rep(NA,t)}
if(is.null(commax)){commax<-rep(NA,t)}
if(is.null(comnot)){comnot<-as.list(rep(NA,t))}
if(is.null(comnot.type)){comnot.type<-as.list(rep("single",t))}
if(is.null(dups)){dups<-0}
traits.final<-as.data.frame(tr.pool)

# interpret traits matrix
tr.names<-paste("TR", c(1:(ncol(traits.final)-1)), sep="")
dimnames(traits.final)[[2]]<-c(tr.names, "LIKELIHOOD")
if(is.null(dimnames(traits.final)[[1]])){
		sp.names<-paste("SP", 1:nrow(traits.final), sep="")
		dimnames(traits.final)[[1]]<-sp.names
}else{
		duptraits<-grepl("^DP[[:digit:]]", dimnames(traits.final)[[1]])
		if(sum(duptraits)==nrow(traits.final)/2){ 
			traits.final<-traits.final[duptraits==FALSE, ]
		}
}
p<-nrow(traits.final)

# check if the abundance distribution specified includes a standard deviation
if(is.null(abun.method)){abun.method<-"lnorm"}
stopper<-regexpr(":", abun.method, fixed=TRUE)
abun.dist<-ifelse(stopper==-1, abun.method, substr(abun.method, start=0, stop=stopper-1))
abun.dist<-match.arg(abun.dist, c("unif", "norm", "lnorm", "fixed"), several.ok=FALSE)
abun.sd<-as.numeric(ifelse(stopper==-1, 1, substr(abun.method, start=stopper+1, stop=nchar(abun.method))))

# Now sample species from species pool for each community, based on
# species richness specified in s, and repeated as specified in r

# First, create a list of abundance vectors, with enough elements to have a 
# separate vector of length p for each combination of species richness and
# replicate
abun <- list(rep(0, p))
if(dups>0){abun <- list(rep(0, 2*p))}
abun <- rep(abun, r * l.s)

# Next determine which species should be selected for non-zero abundances.
# Check how communities should reflect distributions specified for the overall
# species pool, and how abundances should be distributed among the species
# selected.
traits.prob<-traits.final$LIKELIHOOD
if(presence=="random"){pres.prob<-rep(1,p)}
if(presence=="common"){pres.prob<-traits.prob}
if(presence=="rare"){pres.prob<-rescale(1/traits.prob)}
if(abundance=="random"){abun.prob<-rep(1,p)}
if(abundance=="common"){abun.prob<-traits.prob}
if(abundance=="rare"){abun.prob<-rescale(1/traits.prob)}

# Also enforce specified minima and maxima and excluded values
if(sum(is.na(commin))<t || sum(is.na(commax))<t || sum(is.na(comnot))<t){
		for(i in 1:t){
if(!is.na(commin[i])){pres.prob<-ifelse(as.numeric(as.character(traits.final[,i]))<commin[i], 0, pres.prob)}
if(!is.na(commax[i])){pres.prob<-ifelse(as.numeric(as.character(traits.final[,i]))>commax[i], 0, pres.prob)}
			if(comnot.type[[i]]=="single"){
if(sum(!is.na(comnot[[i]]))>0){pres.prob<-ifelse(is.element(as.numeric(as.character(traits.final[,i])), as.numeric(as.character(comnot[[i]]))), 0, pres.prob)}
			}else{
				if(sum(!is.na(comnot[[i]]))>0){
					istart<-min(comnot[[i]], na.rm=TRUE)
					istop<-max(comnot[[i]], na.rm=TRUE)
pres.prob<-ifelse(findInterval(as.numeric(as.character(traits.final[,i])), c(istart, istop), rightmost.closed=TRUE)==1, 0, pres.prob)
				}
			}
		}
}
if(sum(pres.prob!=0)<=max(s))
{stop("Values in commin and/or commax and/or comnot are too restrictive!")
}

# adjust length of likelihood vectors if duplicate species are required
if(dups>0){pres.prob<-c(pres.prob,rep(0,p))}
if(dups>0){abun.prob<-rep(abun.prob,2)}

# assign abundance values in each element of abun to the relevant number of 
# species picked from the species pool, taking into account duplicates and  
# how species and their abundances should reflect the pool's trait distributions. 
abun<-mapply(fill.abun, empty=abun, comsize=nb.sp, MoreArgs=list( presprob=pres.prob, abunprob=abun.prob, abundist=abun.dist, abunsd=abun.sd, ids=dups))

# ensure that abundace values are positive (normal distribution can yield
# negative values, which may be problematic in further analysis)
if(min(abun)<0){abun<-apply(X=abun, MARGIN=2, FUN=pos.abun)}

# transpose the abundance matrix for a site-by-species matrix with rows
# equating to communities, and columns to species
abun <- t(abun)

# Now assign row and column names to the community matrix
names.sp<-dimnames(traits.final)[[1]]
if(dups>0){
		names.sp<-c(names.sp, paste("DP", c(1:p), sep=""))
		traits.final<-rbind(traits.final, traits.final)
		dimnames(traits.final)[[1]]<-names.sp
		traits.prob<-rep(traits.prob,2)
}
names.com<-paste("COM", c(1:(r * l.s)), sep = "")
dimnames(abun)<-list(names.com, names.sp)

# Finally, track min, max, range, skewness and kurtosis in trait values and 
# trait likelihood for the species pool and per community, plus user input values
comstats<-as.data.frame(matrix(NA, nrow=1+r*l.s, ncol=20+t*13))
dimnames(comstats)[[2]]<-c("Site", "SR_pool", "SR_site", "N_dups", "N_traits",
paste("T", rep(1:t, each=13), c("_dist", "_sd", "_type", "_levels", "_setmin", "_obsmin", "_setmax", "_obsmax", "_range", "_unique.values", "_var", "_skewness", "_kurtosis"), sep=""),
"likeli_min", "likeli_max", "likeli_range", "likeli_var", "likeli_skewness", "likeli_kurtosis", "presence", "abundance", "abun_dist", "abun_sd", "abun_min", "abun_max", "abun_var", "abun_skewness", "abun_kurtosis")
comstats$Site<-c("pool", names.com)
comstats$SR_pool<-p
comstats$SR_site<-c(p, nb.sp)
comstats$N_dups<-dups
comstats$N_traits<-t
comstats$presence<-presence
 comstats$abundance<-abundance
comstats$abun_dist<-abun.dist
comstats$abun_sd<-abun.sd
for(i in 0:(r*l.s)){
		#i<-2 	# for testing
		# deal with the special case of stats for the species pool
		if(i==0){
			indx<-c(1:nrow(traits.final))
			comstats$likeli_min[i+1]<-min(traits.prob[indx])
			comstats$likeli_max[i+1]<-max(traits.prob[indx])
			comstats$likeli_range[i+1]<-comstats$likeli_max[i+1]-comstats$likeli_min[i+1]
			comstats$likeli_var[i+1]<-var(traits.prob[indx])
			comstats$likeli_skewness[i+1]<-skewness(traits.prob[indx])
			comstats$likeli_kurtosis[i+1]<-kurtosis(traits.prob[indx])
			comstats$abun_min[i+1]<-NA
			comstats$abun_max[i+1]<-NA
			comstats$abun_var[i+1]<-NA
			comstats$abun_skewness[i+1]<-NA
			comstats$abun_kurtosis[i+1]<-NA
		}else{
			indx<-which(abun[i,]>0)
			comstats$likeli_min[i+1]<-min(traits.prob[indx])
			comstats$likeli_max[i+1]<-max(traits.prob[indx])
			comstats$likeli_range[i+1]<-comstats$likeli_max[i+1]-comstats$likeli_min[i+1]
			comstats$likeli_var[i+1]<-var(traits.prob[indx])
			comstats$likeli_skewness[i+1]<-skewness(traits.prob[indx])
			comstats$likeli_kurtosis[i+1]<-kurtosis(traits.prob[indx])
			comstats$abun_min[i+1]<-min(abun[i,indx])
			comstats$abun_max[i+1]<-max(abun[i,indx])
			comstats$abun_var[i+1]<-var(abun[i,indx])
			comstats$abun_skewness[i+1]<-skewness(abun[i,indx])
			comstats$abun_kurtosis[i+1]<-kurtosis(abun[i,indx])
			
		} # end if-else

		for(j in 1:t){
			#j<-1  # for testing
names.col<-paste("T", rep(j,each=13), c("_dist", "_sd", "_type", "_levels", "_setmin", "_obsmin", "_setmax", "_obsmax", "_range", "_unique.values", "_var", "_skewness", "_kurtosis"), sep="")
			comstats[i+1,names.col[1]]<-"??"
			comstats[i+1,names.col[2]]<-"??"
			comstats[i+1,names.col[3]]<-"??"
			comstats[i+1,names.col[4]]<-"??"
			comstats[i+1,names.col[5]]<-commin[j]
			comstats[i+1,names.col[7]]<-commax[j]
			
			# compute observed min, max, range, skewness and kurtosis per trait
			comstats[i+1, names.col[6]]<-min(as.numeric(as.character(traits.final[indx,j])))
			comstats[i+1, names.col[8]]<-max(as.numeric(as.character(traits.final[indx,j])))
comstats[i+1, names.col[9]]<-comstats[i+1, names.col[8]]-comstats[i+1, names.col[6]]
comstats[i+1, names.col[10]]<-length(unique(as.numeric(as.character(traits.final[indx,j]))))
			comstats[i+1, names.col[11]]<-var(as.numeric(as.character(traits.final[indx,j])))
comstats[i+1, names.col[12]]<-skewness(as.numeric(as.character(traits.final[indx,j])))
comstats[i+1, names.col[13]]<-kurtosis(as.numeric(as.character(traits.final[indx,j])))
		} # end for j
} # end for i
return(list(T=traits.final, A=abun, S=comstats))   
} # end function simul.comms.2