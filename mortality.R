##################################
#Script to estimate mortality for copepodite stage-pairs using a statistical regression approach 
#Created by Kristina Kvile, August 2015
#Written for R (https://www.r-project.org/)
#Requires the library mgcv (Mixed GAM Computation Vehicle with GCV/AIC/REML Smoothness Estimatio,n, S.Wood)
#https://cran.r-project.org/web/packages/mgcv/index.html

library(mgcv)
load("dat.rda") #Load data
#Data contain the following variables: 
#Sampling time and position:
# 1.Month
# 2.Year
# 3.Lat (Latitude, decimal degree)
# 4.Lon (Longitude, decimal degree)
# 5.Day (Julian day of sampling)
# 5.Season (Spring or summer)
# 6.Temp (Ambient temperature at station)
#Abundance of copepodite stages at station:
# 7.CI
# 8.CII
# 9.CIII
# 10.CIV
# 11.CV
# 12. CVI   
stages<-c("CI","CII","CIII","CIV","CV","CVI") #Copepodite stages
trans<-c("CI-CII","CII-CIII","CIII-CIV") #Copepodite stage pairs for which we estimate mortality
mort_true<-c(0.09,0.105,0.075,0.03,0.025,0.02) #Stage-specific mortalities per copepodite stage (CI-CVI) used in the simulation
mort_avg<-sapply(1:(length(mort_true)-1), function(i) mean(mort_true[i:(i+1)])) #Averaged per stage pair

#Coefficients to estimate development time per stage (Corkett, 1986)
aI<- 6419;  aII<-8014;  aIII<-9816;  aIV<-11601;  aV<-13526;  aVI<-17477;  alfa<-10.6;  b<-(-2.05)

#Estimate stage-specific ages from day-of-spawning based on temperature at sampling station ("Temp")
#Age of stage i is the midpoint between the age of entry until stage i and i+1
for (i in c(1:dim(dat)[1])){ 
  dat$Age.CI[i]<-median(seq(aI*(dat$Temp[i]+alfa)^b,aII*(dat$Temp[i]+alfa)^b,length.out = 10))
  dat$Age.CII[i]<-median(seq(aII*(dat$Temp[i]+alfa)^b,aIII*(dat$Temp[i]+alfa)^b,length.out = 10))
  dat$Age.CIII[i]<-median(seq(aIII*(dat$Temp[i]+alfa)^b,aIV*(dat$Temp[i]+alfa)^b,length.out = 10))
  dat$Age.CIV[i]<-median(seq(aIV*(dat$Temp[i]+alfa)^b,aV*(dat$Temp[i]+alfa)^b,length.out = 10))
  dat$Age.CV[i]<-median(seq(aV*(dat$Temp[i]+alfa)^b,aVI*(dat$Temp[i]+alfa)^b,length.out = 10))
}

#Estimate stage-specific day-of-spawning as the day of sampling minus the estimated age
dat$Spd.CI<-dat$Day-dat$Age.CI;  
dat$Spd.CII<-dat$Day-dat$Age.CII; 
dat$Spd.CIII<-dat$Day-dat$Age.CIII; 
dat$Spd.CIV<-dat$Day-dat$Age.CIV; 
dat$Spd.CV<-dat$Day-dat$Age.CV; 

#Estimate stage-durations from stage i to i+1
dat$Duration.CI<-aII*(dat$Temp+alfa)^b    -aI*(dat$Temp+alfa)^b; 
dat$Duration.CII<-aIII*(dat$Temp+alfa)^b  -aII*(dat$Temp+alfa)^b;  
dat$Duration.CIII<-aIV*(dat$Temp+alfa)^b  -aIII*(dat$Temp+alfa)^b;
dat$Duration.CIV<-aV*(dat$Temp+alfa)^b    -aIV*(dat$Temp+alfa)^b; 
dat$Duration.CV<-aVI*(dat$Temp+alfa)^b    -aV*(dat$Temp+alfa)^b;

###########################################
#Statistical Regression Approach (SRA)
#Covariates:
#1. Day of spawning = sampling day minus stage-specific age at middle of the stage -> Seasonal variation in spawning time
#2. Sampling location (Lon,Lat) -> Horizontal transport from spawning locatoin 
#3. Approx. stage-specific age, centered around zero -> Estimates mortality (-m)

#Matrix to store mortality estimates per stage-pair (CV-CVI not included)
mortalities<-matrix(NA,ncol=length(trans)) #Overall mortality
colnames(mortalities)<-trans

#Estimate mortality per stage pair
for (i in 1:length(trans)) {  
  #Create data-frame of to successive stages, i and i+1
  combined.abundance<-data.frame(cbind(dat[,stages[i]],dat[,stages[i+1]])) 
  combined.ages<-data.frame(cbind(dat[,paste0("Age.",stages[i])],dat[,paste0("Age.",stages[i+1])])) 
  combined.spds<-data.frame(cbind(dat[,paste0("Spd.",stages[i])],dat[,paste0("Spd.",stages[i+1])])) 
  combined.duration<-data.frame(cbind(dat[,paste0("Duration.",stages[i])],dat[,paste0("Duration.",stages[i+1])]))
  Data<-data.frame(cbind(rep(dat$Year,2), #Repeat common variables for the two stages observed in the same station
                         rep(dat$Day,2),
                         rep(dat$Lon,2),
                         rep(dat$Lat,2),
                  stack(combined.abundance)[1],#Add information for stage i and i+1 in the same column
                  stack(combined.ages)[1],
                  stack(combined.spds)[1]),
                  stack(combined.duration)[1]) 
  colnames(Data)<-c("Year","Day","Lon","Lat","Abundance","Age","Spd","Dur")  
  Data<-Data[!is.na(Data$Abundance) & Data$Abundance>0,] #Remove samples with zero abundance or missing values
  Data$Year<-as.factor(Data$Year) #Make sure year is a factor
  Data$Abundance<-Data$Abundance/Data$Dur #Scale abundance by stage duration
  Data$Age<-scale(Data$Age,center=TRUE,scale=FALSE)  #Center the age variable around zero
  Data$logAbundance<-log(Data$Abundance) #Log transform abundance
  m <- gam(Data$logAbundance~s(Spd)+te(Lon,Lat,k=5)+s(Year,bs="re")+Age+s(Year,bs="re",by=Age),data=Data) #SRA model
  mortalities[i]<- -coef(m)[2] #The negative of the age-coefficient gives mortality
}

###########################################
#Bootstrap procedure for the SRA estimates
B <- 1000 #The number of boostrap samples (note that 1000 iterations takes some time...)
bootmat<-matrix(NA,nrow=2,ncol=length(trans)) #Matrix to store bootstrap results
colnames(bootmat)<-trans
rownames(bootmat)<-c("Up","Low") #Upper and lower confidence limits for model coefficients

#Do the bootstrap per stage-pair
for (i in c(1:4)) {  
  #Data-frame of to successive stages, same as above
  combined.abundance<-data.frame(cbind(dat[,stages[i]],dat[,stages[i+1]]))
  combined.ages<-data.frame(cbind(dat[,paste0("Age.",stages[i])],dat[,paste0("Age.",stages[i+1])]))
  combined.spds<-data.frame(cbind(dat[,paste0("Spd.",stages[i])],dat[,paste0("Spd.",stages[i+1])]))
  combined.duration<-data.frame(cbind(dat[,paste0("Duration.",stages[i])],dat[,paste0("Duration.",stages[i+1])]))
  Data<-data.frame(cbind(rep(dat$Year,2),
                         rep(dat$Day,2),
                         rep(dat$Lon,2),
                         rep(dat$Lat,2),
                         stack(combined.abundance)[1],
                         stack(combined.ages)[1],
                         stack(combined.spds)[1]),
                     stack(combined.duration)[1])
  colnames(Data)<-c("Year","Abs.day","Lon.dec","Lat.dec","Abundance","Age","Spd","Dur")
  Data<-Data[!is.na(Data$Abundance) & Data$Abundance>0,] 
  Data$Abundance<-Data$Abundance/Data$Dur 
  Data$Year.num<-as.numeric(as.character(Data$Year)) 
  Data$Age<-scale(Data$Age,center=TRUE,scale=FALSE)  
 
  #Doing the bootstrap with year as the sampling unit
  yrs <- as.character(levels(Data$Year))
  #Splitting the data by year
  Data.list <- split(Data,Data$Year)
  names(Data.list) <- yrs
  #To store the resulting bootstrap vectors:
  boot.av.m <- NULL

  for(b in 1:B){
    #Sample years
    yrs.B <- sample(yrs,size=length(yrs),replace=T)
    #Constructing a bootstrap data set of the sampled years:
    Data.b <- NULL
    for(j in 1:length(yrs.B)){
      Data.b.j <- Data.list[[yrs.B[j]]] #Taking the data for year #j in the sample
      Data.b.j$year <- yrs[j]  #Renaming the year the #j in the original order
      Data.b <- rbind.data.frame(Data.b,Data.b.j) #Creating a dataframe of samples
    }
    #Redefining the year variable (same number of levels as original, some might have zero occurence):
    Data.b$Year <- factor(x=as.character(Data.b$year),levels=levels(Data$Year))
    #Estimate mortality from the resampled dataset 
    m.b <- gam(log(Data.b$Abundance)~s(Spd)+te(Lon.dec,Lat.dec)+s(Year,bs="re")+Age+s(Year,bs="re",by=Age),data=Data.b)
    av.m <- coef(m.b)["Age"] 
    #Store the estimates
    boot.av.m <- c(boot.av.m,av.m)
  }
  boot.av.m <-(-boot.av.m) #The negative of the estimates give mortality
    
  ### Calculation of upper and lower confidence limits for model coefficients (for standard errors, replace 'quantile(x,0.025)' with 'sd(x)')
  bootmat["Low",i] <- quantile(boot.av.m,0.025)
  bootmat["Up",i] <- quantile(boot.av.m,0.975)
}
