
library(dplyr)
library(psych)
library(stringr)
##Set the working directory
#setwd('~/Dropbox (ASU)/Apr17')
setwd('C:/Users/fsa_0/Dropbox/Dropbox (ASU)/us-covid-19-data-master/April27')


###Import the original data
covid<-read.csv('us_covid19_confirmed.csv', check.names=F)
env<-read.csv('env_means_2020.02.24_to_2020.03.11.csv')
hum<-read.csv('Human_data.csv')
lock<-read.csv('Lockdown_Mar_18.csv')

##Data preparation
covid.merge1<-merge(covid,env,by=('FIPS'))
covid.merge2<-merge(covid.merge1,hum,by=('FIPS'))
covid.merge3<-merge(covid.merge2,lock,by=('Province_State'))

#Variables
variables<-covid.merge3[,c(2,115:127)]
#Remove all data but number of cases and FIPS
dat<-covid.merge3[,c(-1,-3:-11,-115:-127)]


##Coordinates
covid.coord<-covid.merge3[,c(2,10,9)]
#Lat = Y Long = X
names(covid.coord) <- c('FIPS','x','y')
str(covid.coord)




dates <- names((dat[,-1]))
dates1<-str_replace_all(dates, ('/'), '_')


end_dates<-c(57,67,77,87,97)

res3<-list()



#for (j in 1:length(end_dates)){

j=1
	#Select the target day
	sel_dat<-dat[,-1:end_dates[j]+1]
	#Calculate the maximum and count columns
	sel_dat2<-sel_dat[,-1]
	
	sel_dat$Count<-apply(sel_dat2,1,function(x) sum(x > 0))

	sel_dat[, "Max"] <- apply(sel_dat2, 1, max)
	#Calulate the number of columns
	b1<-ncol(sel_dat)
	#Select FIP, Max and county colummns
	max_count<-sel_dat[,c(1,b1,(b1-1))]	

	
	
	#Create a vector with FIPS code
	fips<-sel_dat$FIPS




	###Create a list to store the results 
	##OLS and geometric mean
	lc<-list()
	r_sq<-list()
	fip <- list()
	geo<-list()
	coun<-list()
	##Correlation between r and Lnr
	corr<-list()
	##Results from RF mmodels (Temperrature and Humidity)
	rf_rq_hum<-list()
	rf_rq_tem<-list()
	rf_mse_hum<-list()
	rf_mse_tem<-list()


	res3<-list()
	##Results from Autocorrelation analysis
	IncMSE.hum<-list()
	IncMSE.tem<-list()



		for (i in 1:nrow(sel_dat)){
		#Calulate the number of columns
		b<-ncol(sel_dat)
		
		#Select the number of cases for each county. Remove the Max and county colummns
		c1<-as.numeric(sel_dat[i,c(-b,-(b-1))])		
		#Select the cases greather or equal than 10
		y<-c1[-1][c1[-1]>=10]
	
		#Select counties with more than 5 cases
			if(length(y)>5){
			#Create a vector of days
			x<-seq(1:length(y))
			#Combine number of cases and days
			data.df<-data.frame(x,y)
			#Fit the OLS model
			model.0 <- lm(log(y) ~ x, data=data.df) 
			#Extract the Rsquared
			r_sq_lm_slope = summary(model.0)$r.squared
			#Extract the coefficients
			#Intercept
			alpha.0 <- (coef(model.0)[1])
			#Beta
			beta.0 <- coef(model.0)[2]
			##Calculate the geometric mean
			geo_mean<-geometric.mean(y/lag(y),na.rm=TRUE) 
		
			#Store the results in lists
			lc[[i]]<- rbind(beta.0) 
    		r_sq[[i]]<-rbind(r_sq_lm_slope)
    		geo[[i]]<- rbind(geo_mean)
			fip[[i]]<-rbind(fips[i])}
		}
	

	###Prepare the table with the results
	##Beta
	dc<- do.call(rbind, lc)
	dc <- data.frame(dc)

	####Rsquared
	r_sq1<- do.call(rbind, r_sq)
	r_sq1<- data.frame(r_sq1)

	###Geometric mean
	dgeo<- do.call(rbind, geo)
	dgeo <- data.frame(dgeo)


	##FIPS
	fips2<- do.call(rbind, fip)
	fips2 <- data.frame(fips2)


	res<-cbind(fips2,dc,r_sq1,dgeo)
	names(res)<-c('FIPS','Ln_r','R.sqrd','r')


	###Remove all cases with an rsquared lower that 0.8	
	optimized<-subset(res,R.sqrd>0.8)
	##Combine the datasets with the variables
	merge.r<-merge(variables,optimized,by=('FIPS'),all.x=T)
	merge.r2<-merge(merge.r,max_count,by=('FIPS'))

	#For counties with no cases or growth, add 0.
	#Remove all other counties woth NAs
	m0<-subset(merge.r2,r>0)
	##Select counties with zero and one case
	m1<-subset(merge.r2,Max==0 | Max==1)
	##Add 0 to NAs
	m1$r <-0
	#Combine counties with transmission and no transmission
	m2<-rbind(m0,m1)

	#Save the output
	#write.csv(m3,'Results_zero_counties.csv',row.names=F)
	write.csv(m2,paste('Matrix_',dates1[end_dates[j]],'.csv',sep=''),row.names=F)
	
	###Correlation between r and Ln_r
	corr[[j]]<-rbind(cor(m0$r,m0$Ln_r))

    #m3<-m2[order(m2$FIPS),]
	m3<-merge(covid.coord,m2,by=('FIPS'))
	
	#data to create the autocorrelated variable
	covid.coordinates <- m3[, c("x", "y")]

	#data to fit the model
	#covid.hum <- m3[, c("r", "Precipitation", "Specific.humidity", "Population", "Route.density", "House.density", "Lockdown")]
	#covid.tem <- m3[, c("r", "Temperature", "Population", "Route.density", "House.density", "Lockdown")]

	#Data with no Lockdown
	covid.hum.nolock <- m3[, c("r", "Precipitation", "Specific.humidity", "Population", "Route.density", "House.density")]
	covid.tem.nolock <- m3[, c("r", "Temperature.C", "Population", "Route.density", "House.density")]
	
	
	###################################################################
	###################################################################
	###################################################################
	################# Modeling
	###################################################################
	###################################################################
	###################################################################

	#Data frame for RF analysis
	

	#testing multicollinearity in model
	#HH::vif(m3[,c("accu_conv_prec_at_surf","spec_humi_at_2m","Population","Dens_housing","Route_den","Lockdown_days_as_of_April_15")])
	#HH::vif(covid.rf[,c("Precipitation","Specific.humidity","Population","House.density","Route.density","Lockdown")])
	#Precipitation Specific.humidity        Population     House.density     Route.density          Lockdown 
    #     1.485115          1.522731          1.112809          1.129548          1.064517          1.030594
 
	library(randomForest)

###############################
	###Random Forests all cases
###############################

	set.seed(20)
	rf_hum<-randomForest(r~ Population+
					Precipitation+
					Specific.humidity+
					 House.density +
					Route.density,data=covid.hum.nolock, importance=T)
	
	rf_hum	
	
	set.seed(20)
	rf_tem<-randomForest(r~ Population+
					Temperature.C+
					 House.density +
					Route.density,data=covid.tem.nolock, importance=T)
	rf_tem	

	rf_rq_hum[[j]]<-rbind(rf_hum$rsq[500])
	rf_mse_hum[[j]]<-rbind(rf_hum$mse[500])
	rf_rq_tem[[j]]<-rbind(rf_tem$rsq[500])
	rf_mse_tem[[j]]<-rbind(rf_tem$mse[500])


	##Save the importance graph
	pdf(paste('Importance_all_',dates1[end_dates[j]],'.pdf',sep=''), width=10, height=5)
	par(mfrow=c(1,2))
	varImpPlot(rf_hum,type=1)
	varImpPlot(rf_tem,type=1)
	dev.off()

	##Extractor function for variable importance measures as produced by randomForest
	#imp<-importance(rf)

	##Extractor function for variable importance measures as produced by randomForest
	imp_hum<-importance(rf_hum)
	imp_tem<-importance(rf_tem)

	###Save the partial dependet plots
	impvar_hum <- rownames(imp_hum)[order(imp_hum[, 1], decreasing=TRUE)]
	impvar_tem <- rownames(imp_tem)[order(imp_tem[, 1], decreasing=TRUE)]

	pdf(paste('Partial_dependent_plots_hummidity_all_',dates1[end_dates[j]],'.pdf',sep=''))
	op <- par(mfrow=c(2, 3))
		for (m in seq_along(impvar_hum)) {
		partialPlot(rf_hum, covid.hum.nolock, main=NULL, impvar_hum[m], ylab='Fitted function',xlab=impvar_hum[m],
		ylim=c(0,1.1))}
	dev.off()

	pdf(paste('Partial_dependent_plots_temperature_all_',dates1[end_dates[j]],'.pdf',sep=''))
	op <- par(mfrow=c(2, 3))
		for (m in seq_along(impvar_tem)) {
		partialPlot(rf_tem, covid.tem.nolock, main=NULL, impvar_tem[m], ylab='Fitted function',xlab=impvar_tem[m],
		ylim=c(0,1.1))}
	dev.off()

	###Prepare the table with the results
	#co
	dcor<- do.call(rbind, corr)
	dcor <- data.frame(dcor)

	#Rsquared
	r_rq_hum<- do.call(rbind, rf_rq_hum)
	r_rq_hum<- data.frame(r_rq_hum)
	r_rq_tem<- do.call(rbind, rf_rq_tem)
	r_rq_tem<- data.frame(r_rq_tem)

	
	#MSE
	r_mse_hum<- do.call(rbind, rf_mse_hum)
	r_mse_hum<- data.frame(r_mse_hum)
	r_mse_tem<- do.call(rbind, rf_mse_tem)
	r_mse_tem<- data.frame(r_mse_tem)

	
	res3[[j]]<-cbind(dates1[end_dates[j]],dcor,r_rq_hum,r_mse_hum,r_rq_tem,r_mse_tem)
	

				

	##############################################################
	##############################################################
	###RRANDOM FOREST MODEL WITH RANDOM VARIABLE ON EACH ITERATION
	##############################################################
	##############################################################

library(sp)
library(gstat)


	#HOW DO WE CREATE AUTOCORRELATED SURFACES?
############################################
#this is named "unconditional gaussian simulation". See more info here: http://santiago.begueria.es/2010/10/generating-spatially-correlated-random-fields-with-r/

#creating the raster area using coordinate range of covid data
raster.area <- expand.grid(
  x = seq(min(covid.coordinates$x) - 1, max(covid.coordinates$x) + 1, by = 1),
  y = seq(min(covid.coordinates$y) - 1, max(covid.coordinates$y) + 1, by = 1)
)

#building a dummy autocorrelation model
autocorrelation.model <- gstat(
  formula = z~1, 
  locations = ~x+y, 
  dummy = T, 
  beta = 1, 
  model = vgm(
    psill = 0.025, 
    range = sample(1:50, 1), #range (distance) of the spatial autocorrelation, in this case given in degrees. Indicates the distance at which data-points stop being autocorrelated
    model = 'Exp'
    ), 
  nmax=20
  )
#predicting the model to raster.area (this prediction can be done directly to covid.coordinates, we will use it that way to fit the models)
autocorrelated.map <- predict(
  autocorrelation.model, 
  newdata = raster.area, 
  nsim = 1
  )
gridded(autocorrelated.map) = ~x+y
spplot(obj=autocorrelated.map[1])


#function to generate autocorrelated variable
autocorrelated_variable <- function(covid.coordinates){
  
  #building a dummy autocorrelation model
  autocorrelation.model <- gstat(
    formula = z~1, 
    locations = ~x+y, 
    dummy = T, 
    beta = 1, 
    model = vgm(
      psill = 0.025, 
      range = sample(1:50, 1), #range of the spatial autocorrelation, max set to 50 because for this area and resolution seems to provide good outcomes
      model = 'Exp'
    ), 
    nmax=20
  )
  #predicting the model to raster.area (this prediction can be done directly to covid.coordinates, we will use it that way to fit the models)
  autocorrelated.variable <- predict(
    autocorrelation.model, 
    newdata = covid.coordinates, 
    nsim = 1
  )
  
  return(autocorrelated.variable$sim1)
  
}

library(viridis)
#RANDOM FOREST MODEL WITH RANDOM VARIABLE ON EACH ITERATION
iterations <- 1000
IncMSE.hum <- list()
IncMSE.tem <- list()

for(i in 1:iterations){

  #OLD METHOD
  #adds autocorrelated noise (the autocorrelation gets a random length from 100 to 1000 cases)
  # covid$autocor.noise <- as.vector(stats::filter(runif(nrow(covid)), 
  #        filter = rep(1, sample(100:1000, 1)), 
  #        circular = TRUE))
  
  #NEW METHOD
  covid.hum.nolock$autocor.noise <- autocorrelated_variable(covid.coordinates)
  covid.tem.nolock$autocor.noise <- autocorrelated_variable(covid.coordinates)

  #fits model (ranger is WAY faster, but you can put here the other randomForest)
  m.hum <- ranger::ranger(
    data = covid.hum.nolock, 
    dependent.variable.name = "r",
    num.trees = 500,
    min.node.size = 10,
    mtry = 3,
    importance = "permutation",
    scale.permutation.importance = TRUE,
    num.threads = 8 #set here your number of cores
    )
    
    m.tem <- ranger::ranger(
    data = covid.tem.nolock, 
    dependent.variable.name = "r",
    num.trees = 500,
    min.node.size = 10,
    mtry = 3,
    importance = "permutation",
    scale.permutation.importance = TRUE,
    num.threads = 8 #set here your number of cores
    )

    
  
  #getting importance components
  IncMSE.hum[[i]] <- m.hum$variable.importance
  IncMSE.tem[[i]] <- m.tem$variable.importance


}

#to dataframe
IncMSE.df.hum <- as.data.frame(do.call("rbind", IncMSE.hum))
IncMSE.df.tem <- as.data.frame(do.call("rbind", IncMSE.tem))


#boxplots
pdf(paste('RSME_all_',dates1[end_dates[j]],'.pdf',sep=''))
par(mfrow = c(2, 1), mar = c(3, 12, 4, 2) + 0.1)
boxplot(IncMSE.df.hum, 
        main = "Variable importance", 
        horizontal = TRUE, 
        las = 1, 
        cex.axis = 1, 
        cex.lab = 1.2,
        xlab = "% Increment in MSE", 
        notch = TRUE,
        col = viridis(ncol(covid.hum.nolock), alpha = 0.7)
        )
abline(v = quantile(IncMSE.df.hum$autocor.noise, probs = c(0.5, 1)), col = "gray50", lwd = c(3,2), lty = c(1,2))

boxplot(IncMSE.df.tem, 
        main = "Variable importance", 
        horizontal = TRUE, 
        las = 1, 
        cex.axis = 1, 
        cex.lab = 1.2,
        xlab = "% Increment in MSE", 
        notch = TRUE,
        col = viridis(ncol(covid.tem.nolock), alpha = 0.7)
        )
abline(v = quantile(IncMSE.df.tem$autocor.noise, probs = c(0.5, 1)), col = "gray50", lwd = c(3,2), lty = c(1,2))
dev.off()


results<- do.call(rbind, res3)

write.csv(results,paste('Results_RF_models',dates1[end_dates[j]],'.csv',sep=''),row.names=F)



