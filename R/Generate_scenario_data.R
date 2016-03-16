	                                                     
##################
# Main function to generate data (the OM)
# Code written by Dr. Kotaro Ono
##################
Generate_scenario_data <- function(Sim_Settings){
		
	# Libraries
  if (!require(pacman)) install.packages("pacman")
	pacman::p_load(MASS, RandomFields, fields, geoR, kml, fpc, gtools, tweedie)

  ## generate mvt matrix of an animal from a cell to the other 
	Mvt_matrix <- function(Pop_param, Depth_eff, Dist_eff, Lat_eff, Depth_func, Dist_func="Exp", Lat_func, x,y,data.bathym){
		Mvt_mat <- matrix(0,(length(x)*length(y)),(length(x)*length(y)))
		for (loc in 1:(length(x)*length(y))){ 
			Mvt_mat[loc,] <- Mvt_upd1(loc, Pop_param, Depth_eff=Depth_eff, Dist_eff=Dist_eff, Lat_eff=Lat_eff, Depth_func=Depth_func, Dist_func=Dist_func, Lat_func=Lat_func,x,y, data.bathym)
		}	
		return(Mvt_mat)
	}
		
  ## for updating the population distribution according to different different covariates (depth and/or dist), using joint prob
	Mvt_upd1 <- function(loc, Pop_param, Depth_eff, Dist_eff, Lat_eff, Depth_func, Dist_func="Exp", Lat_func, x,y, data.bathym){
		Di1 <- Pop_param[1]  	# 1st value for dist function
		Di2 <- Pop_param[2]  	# 2nd value for dist function
		De1 <- Pop_param[3]  	# 1st value for depth function
		De2 <- Pop_param[4]	 	# 2nd value for depth function
		Lat1 <- Pop_param[5]  	# 1st value for lat/long function
		Lat2 <- Pop_param[6]	# 2nd value for lat/long function
		
		a1 <- which(abs(data.bathym$x-data.bathym$x[loc])<=length(x)/2)
		a2 <- which(abs(data.bathym$x-data.bathym$x[loc])>length(x)/2)
		a3 <- which(abs(data.bathym$y-data.bathym$y[loc])<=length(y)/2)
		a4 <- which(abs(data.bathym$y-data.bathym$y[loc])>length(y)/2)
		Data.bathym <- list()
		Data.bathym$x[a1] <- data.bathym$x[a1]
		Data.bathym$x[a2] <- data.bathym$x[a2]#-length(x)
		Data.bathym$y[a3] <- data.bathym$y[a3]
		Data.bathym$y[a4] <- data.bathym$y[a4]#-length(y)
		Dist <- sqrt((Data.bathym$x-Data.bathym$x[loc])^2+(Data.bathym$y-Data.bathym$y[loc])^2)
		
		Prob1 <- 1
		if(Dist_eff=="T"){
			Prob1 <- dexp(Dist, rate=1/Di1)
		}
		
		Prob2 <- 1
		if(Depth_eff=="T"){
			if(Depth_func=="Norm") Prob2 <- dnorm(data.bathym$z, De1, De2*De1)
			if(Depth_func=="Exp") Prob2 <- dexp(data.bathym$z, rate=1/De1)
			if(Depth_func=="Lognorm") { SD <- sqrt(log(De2^2+1)); Prob2 <- dlnorm(data.bathym$z, (log(De1)-SD^2/2), SD) }
			if(Depth_func=="Unif") Prob2 <- dunif(data.bathym$z, De1, De2)
		}
		
		Prob3 <- 1
		if(Lat_eff=="T"){
			if(Lat_func=="Norm") Prob3 <- dnorm(data.bathym$x, Lat1, Lat2)
			if(Lat_func=="Exp") Prob3 <- dexp(data.bathym$x, rate=1/Lat1)
			if(Lat_func=="Lognorm") Prob3 <- dlnorm(data.bathym$x, (log(Lat1)-Lat2^2/2), Lat2)
			if(Lat_func=="Unif") Prob3 <- dunif(data.bathym$x, Lat1, Lat2)
		}
		
		Prob <- Prob1*Prob2*Prob3/sum(Prob1*Prob2*Prob3)
		return(Prob)
	}

  ## Calculate the initial population distribution
	Stable_pop_dist <- function(Mvt_mat_adult, B0){
		Pop1 <- rep(B0/(length(Range_X)*length(Range_Y)), (length(Range_X)*length(Range_Y)))
		Mvt_mat_adult1 <- t(Mvt_mat_adult)
		for (Time in 1:1000){
			Pop <- Mvt_mat_adult1%*%Pop1
			Pop1 <- Pop
		}
		Initpop_adult_mat <- matrix(Pop1, nrow=length(Range_X), ncol=length(Range_Y))
		return(Initpop_adult_mat)
	}

	attach(Sim_Settings)
    on.exit( detach(Sim_Settings) )

	#### Set the bathymetry field
	model_bathym <- RMgauss(var=SD_O^2, scale=SpatialScale)		
	Bathym <- RFsimulate(model_bathym, x=Range_X, y=Range_Y)
	data.bathym <- data.frame(ID = length(Range_X)*length(Range_Y), x = rep(Range_X, length(Range_Y)), y = rep(Range_Y, each=length(Range_X)), z=Bathym@data[,1])
	data.bathym$z = data.bathym$z+abs(min(data.bathym$z))	# to avoid getting negative depth here depth takes positive real number
	image.plot(Range_X, Range_Y, matrix(data.bathym$z,nrow=length(Range_X), ncol=length(Range_Y)))
		
	#### Set the pop mvt param
	Par_mvt_adult <- rbind(Fish_dist_par1, Fish_dist_par2, Fish_depth_par1, Fish_depth_par2, Fish_range_par1, Fish_range_par2);
	Mvt_mat_adult <- lapply(1:n_species, function(xxx) Mvt_matrix(Par_mvt_adult[,xxx], Depth_eff="T", Dist_eff="T", Lat_eff="T", Dist_func=func_mvt_dist[xxx], Depth_func=func_mvt_depth[xxx], Lat_func=func_mvt_lat[xxx], Range_X, Range_Y, data.bathym))
				
	#### Biomass dynamics: similarly to the work from Caruthers, we assume that there is a regional carrying capacity
	# Initial population distribution 
	Pop_adult <- sapply(1:n_species, function(x) Stable_pop_dist(Mvt_mat_adult[[x]], B0[x]))				
	
	if(plotting==TRUE){
		windows()
		nf <- layout(1:4)
		par(mar=c(1,1,1,1))
		fields::image.plot(matrix(Pop_adult[,1],length(Range_X),length(Range_Y)))
		fields::image.plot(matrix(Pop_adult[,2],length(Range_X),length(Range_Y)))
		fields::image.plot(matrix(Pop_adult[,3],length(Range_X),length(Range_Y)))
		fields::image.plot(matrix(Pop_adult[,4],length(Range_X),length(Range_Y)))
	}

	Biomass <- list()
	Biomass[[1]] <- Pop_adult
	Biomass[[1]] <- apply(Biomass[[1]], 2,function(x) replace(x, which(x<=0), 0))
		
	#### Vessel dynamics

	# generate total effort by year (total effort is randomly sampled with log normal error around a logistic type function)
	Sig_effort <- sqrt(log(CV_effort^2+1))
	Effort <- c(); Mean_Effort <- c()		
	Effort_area_year <- list()
	Sig_vessel <- sqrt(log(CV_vessel^2+1))
	qq_vessel <- rlnorm(Nvessels, log(qq_original)-Sig_vessel^2/2, Sig_vessel)		
		
	# use of tweedie distribution to generate catch data for individual boat
	I <- list()	
	xi <- xi	 # that controls for the value of power such that the variance is var[Y] = phi*mu^xi
	phi <- phi	  	
	Catch_area_year <- list(); Catch_area_year_ind <- c()

	#### Updating the population and generate catch data according to the above specifications
	for (iyear in 1:n_years){
		
		### do vessel change their fishing preference over the year (this ONLY controls vessel concentration factor) 
		if(Changing_preference==TRUE) Preference <- 2-exp(-0.1*iyear)
		
		### Assigning total effort in a year		
		Mean_Effort[iyear] <- (Tot_effort/(1+exp(-0.1*iyear)))		# mean effort by year
		Effort[iyear] <- trunc(rlnorm(1, log(Mean_Effort[iyear]-Sig_effort^2/2), Sig_effort))		# effort by year
			
		### Determining the predicted area specific CPUE
		map.size <- length(Range_X)*length(Range_Y)
			
		### changing the catchability "qq" every year depending on whether a specific area is closed or not
		if (iyear>=year.depth.restriction){
			catchability <- rep(1,map.size)
			if(length(depth.restriction)==2) catchability <- replace(catchability, which(bathym >= depth.restriction[1] & bathym <= depth.restriction[2]), 0)
			if(depth.restriction== "random") catchability <- replace(catchability, sample(1:map.size, size=20), 0)
			qq <- catchability
		} else { qq <- qq_original }
		
		### Continuous catch equation
		I[[iyear]] <- Biomass[[iyear]]*(1-exp(-qq)) * matrix(price_fish[iyear,], nrow=map.size, ncol=n_species, byrow=T)
		CPUE_tot <- apply(I[[iyear]],1,sum)
				
		### Distribute the total effort according to the above predicted CPUE
		p <- CPUE_tot^Preference/sum(CPUE_tot^Preference);min(p);max(p)		
		Effort_area_year[[iyear]] <- rmultinom(1, size=Effort[iyear], p=p)
			
		### Then calculate the total catch in each area  
		aa <- c(); AA <- c()			
		Catch_year_area_mat <- matrix(0, nrow=(length(Range_X)*length(Range_Y)), ncol=n_species)
		
		for (xyz in which(Effort_area_year[[iyear]]>0)){
			which.vessel.fished <- sample(1:Nvessels, Effort_area_year[[iyear]][xyz], replace=TRUE)
				
			# This is the code using the continuous function (as suggested by Andre: more realistic)
			if(do.tweedie==TRUE) rand.vessel <- matrix(sapply(1:n_species, function(x) rtweedie(Effort_area_year[[iyear]][xyz], qq_vessel[which.vessel.fished], xi=xi, phi=phi)),ncol=n_species)
			if(do.tweedie==FALSE) rand.vessel <- matrix(sapply(1:n_species, function(x) qq_vessel[which.vessel.fished]),ncol=n_species)
			catch_area_vessel <- sapply(1:n_species, function(Y) { hum <- Biomass[[iyear]][xyz,Y]*(1-exp(-sum(rand.vessel[,Y]))); ifelse(hum==0, catch_area_vessel <- rep(0,Effort_area_year[[iyear]][xyz]), catch_area_vessel <- hum*rand.vessel[,Y]/sum(rand.vessel[,Y])); return(catch_area_vessel)})
			catch_area_vessel <- data.frame(iyear, matrix(rep(c(data.bathym[xyz,-1]), each=Effort_area_year[[iyear]][xyz]),nrow=Effort_area_year[[iyear]][xyz]), which.vessel.fished, matrix(catch_area_vessel,ncol=n_species)); 
			catch_area_vessel <- matrix(sapply( catch_area_vessel, FUN=as.numeric ), nrow=nrow(catch_area_vessel))
			colnames(catch_area_vessel) <- c("year", "X", "Y", "depth", "vessel", paste0("Sp", 1:n_species))
			Catch_year_area_mat[xyz,] <- colSums(catch_area_vessel[,-c(1:5),drop=FALSE])
			
			AA <- rbind(AA, catch_area_vessel)
		}
		
		Catch_area_year_ind <- rbind(Catch_area_year_ind, AA)
		Catch_area_year[[iyear]] <- Catch_year_area_mat
		
   	### Then update the population
		Biomass[[iyear+1]] <- Biomass[[iyear]]+r*Biomass[[iyear]]*(1-Biomass[[iyear]]/Biomass[[1]])-Catch_area_year[[iyear]]
		Biomass[[iyear+1]] <- apply(Biomass[[iyear+1]],2,function(x) replace(x, which(is.na(x)==TRUE | x<=0),0))
		
		### population mvt at the end of year			
		if (no.mvt==FALSE)	Biomass[[iyear+1]] <- sapply(1:n_species, function(x) t(Mvt_mat_adult[[x]])%*%Biomass[[iyear+1]][,x])
		if (no.mvt==TRUE)	Biomass[[iyear+1]] <- Biomass[[iyear+1]]

		Biomass[[iyear+1]] <- apply(Biomass[[iyear+1]],2,function(x) replace(x, which(is.na(x)==TRUE | x<=0),0))

		if (Interactive == TRUE) print(iyear)
	}
		
	#### Plotting 
  if(plotting==TRUE) {
  	# plot total catch 1st year and last year 
		windows()
		nf <- layout(1:2)
		par(mar=c(1,1,2,1), oma=c(4,4,2,2))
		hist(Catch_area_year_ind[,'Sp1'], breaks=100)
		hist(Catch_area_year_ind[,'Sp2'], breaks=100)

	  # plot biomass change over time
		windows()
		for (iyear in 1:n_years){
			nf <- layout(1:5)
			par(mar=c(1,1,2,1), oma=c(4,4,2,2))
			fields::image.plot(matrix(Biomass[[iyear]][,1],length(Range_X),length(Range_X)), main=paste("Year", iyear))
			fields::image.plot(matrix(Biomass[[iyear]][,2],length(Range_X),length(Range_Y)))
			fields::image.plot(matrix(Biomass[[iyear]][,3],length(Range_X),length(Range_Y)))
			fields::image.plot(matrix(Biomass[[iyear]][,4],length(Range_X),length(Range_Y)))
			fields::image.plot(matrix(Effort_area_year[[iyear]] , length(Range_X),length(Range_X)))
			if(Interactive==TRUE) gtools::ask()
		}
			
	  # Time series of effort, raw CPUE, depletion
		windows()
		nf <- layout(1:4)
		par(mar=c(4,4,1,1))
		# Effort
		plot(Effort, ylim=c(min(Effort*0.8), max(Effort*1.2)), lwd=2, type="l", xlab="Year")
		lines(Mean_Effort, lwd=2)
		legend("topleft", legend=c("Sp1","Sp2","Sp3","Sp4"), col=c("black","red","green","blue"), lty=1, bty="n")
		# Total catch
		aa <- t(sapply(1:n_years, function(x) apply(Catch_area_year[[x]],2,sum)))   # total catch 
		matplot(aa, type="l", ylab="Total catch", lty=1, xlab="Year", lwd=2)
		# Raw CPUE (dash line) and biomass (solid line)
		CPUE_raw <- aa/Effort
		bb <- t(sapply(1:n_years, function(x) apply(Biomass[[x]],2,sum)/B0))
		matplot(bb, type="l",ylim=c(0,1), lty=1, lwd=2, ylab="Relative level") 	# the true biomass trajectory in solid line
		RAW_CPUE <- cbind(CPUE_raw[,1]/CPUE_raw[1,1], CPUE_raw[,2]/CPUE_raw[1,2], CPUE_raw[,3]/CPUE_raw[1,3], CPUE_raw[,4]/CPUE_raw[1,4])
		matlines(RAW_CPUE, type="l",ylim=c(0,1), lty=2, lwd=2) 	# the raw CPUE trajectory in dash line
		legend("topright", c("Biomass", "Raw CPUE"), lty=c(1,2), bty="n")
		# Relationship raw CPUE vs depletion
		plot(bb[,1], CPUE_raw[,1]/(CPUE_raw[1,1]), type="l", lwd=2, col=1, xlab="Biomass depletion level", ylab="Raw standardized CPUE", xlim=c(0,1), ylim=c(0,1), )
		lines(bb[,2], CPUE_raw[,2]/CPUE_raw[1,2], type="l", lwd=2, col="red")
		lines(bb[,3], CPUE_raw[,3]/CPUE_raw[1,3], type="l", lwd=2, col="green")
		lines(bb[,4], CPUE_raw[,4]/CPUE_raw[1,4], type="l", lwd=2, col="blue")
		abline(0,1, lty=2)
	}	
		
	#### Return data                                                                # v.names="CPUE", 
	#Raw_data <- reshape(Catch_area_year_ind, varying=paste0("Sp",1:n_species), timevar="Species", times=1:n_species, direction="long") 
	Raw_data <- NULL
	for(i in 1:n_species) Raw_data = rbind(Raw_data, cbind(Catch_area_year_ind[,c("year","X","Y","depth","vessel")], "Species"=i, "CPUE"=Catch_area_year_ind[,paste0("Sp",i)]))
	Bio <- t(sapply(1:n_years, function(x) colSums(Biomass[[x]])))
	Return = list("Data"= Raw_data, "Biomass"=Bio)
	return(Return)
}		
