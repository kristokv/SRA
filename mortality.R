##################################
#Mortality estimation for copepodite stage-pairs using a statistical regression approach (SRA) 
#Created by Kristina Kvile, Aug. 2015
#Written for R (https://www.r-project.org/)
#Requires the mgcv library (Mixed GAM Computation Vehicle with GCV/AIC/REML Smoothness Estimation, S.Wood)
#https://cran.r-project.org/web/packages/mgcv/index.html

library(mgcv) 
load("dat.rda")

#dat contains the following variables: 
#Sampling time and position:
# 1.Month, 2.Year, 3.Latitude (decimal degree), 4. Longitude (decimal degree)
# 5.Julian day of sampling, 6.Season (Spring or summer), 7.Temperature at station
#Abundance of copepodite stages at station:
# 8.CI, 9.CII, 10.CIII, 11.CIV, 12.CV, 13. CVI   

stages<-c("CI","CII","CIII","CIV","CV") #Names of stages 
stage_pairs<-c("CI-CII","CII-CIII","CIII-CIV","CIV-CV") #Names of stage pairs 

#Coefficients to estimate development time per stage (Corkett, 1986)
aI<-6419; aII<-8014; aIII<-9816; aIV<-11601; aV<-13526; aVI<-17477; alfa<-10.6; b<-(-2.05)

#Stage-specific ages: Bêlehrádek temperature functions (D=a(T+alfa)^b)
#Age of stage i considered the midpoint between the age of entry into stage i and i+1
dat$Age.CI<-apply(cbind(aI*(dat$Temp+alfa)^b,aII*(dat$Temp+alfa)^b),1,FUN=mean)
dat$Age.CII<-apply(cbind(aII*(dat$Temp+alfa)^b,aIII*(dat$Temp+alfa)^b),1,FUN=mean)
dat$Age.CIII<-apply(cbind(aIII*(dat$Temp+alfa)^b,aIV*(dat$Temp+alfa)^b),1,FUN=mean)
dat$Age.CIV<-apply(cbind(aIV*(dat$Temp+alfa)^b,aV*(dat$Temp+alfa)^b),1,FUN=mean)
dat$Age.CV<-apply(cbind(aV*(dat$Temp+alfa)^b,aVI*(dat$Temp+alfa)^b),1,FUN=mean)

#Stage-specific day-of-spawning: difference between sampling day and age
dat$Spd.CI<-dat$Day-dat$Age.CI  
dat$Spd.CII<-dat$Day-dat$Age.CII 
dat$Spd.CIII<-dat$Day-dat$Age.CIII 
dat$Spd.CIV<-dat$Day-dat$Age.CIV 
dat$Spd.CV<-dat$Day-dat$Age.CV

#Stage-durations from stage i to i+1
dat$Duration.CI<-aII*(dat$Temp+alfa)^b-aI*(dat$Temp+alfa)^b 
dat$Duration.CII<-aIII*(dat$Temp+alfa)^b-aII*(dat$Temp+alfa)^b  
dat$Duration.CIII<-aIV*(dat$Temp+alfa)^b-aIII*(dat$Temp+alfa)^b
dat$Duration.CIV<-aV*(dat$Temp+alfa)^b-aIV*(dat$Temp+alfa)^b 
dat$Duration.CV<-aVI*(dat$Temp+alfa)^b-aV*(dat$Temp+alfa)^b

#Matrix to store mortality estimates per stage-pair 
mortalities<-matrix(NA,ncol=length(stage_pairs)) 
colnames(mortalities)<-stage_pairs

#Estimate mortality per stage pair
for (i in 1:length(stage_pairs)) {  
  #Create data frame for two successive stages, i and i+1
  Data<-data.frame(
   #1: Repeat variables common for two stages in the same station:
      rep(dat$Year,2),rep(dat$Day,2),rep(dat$Lon,2),rep(dat$Lat,2),
   #2: Stack data for stages i and i+1 in one column:
      c(dat[,stages[i]],dat[,stages[i+1]]), #Abundance
      c(dat[,paste0("Age.",stages[i])],dat[,paste0("Age.",stages[i+1])]), #Age
      c(dat[,paste0("Spd.",stages[i])],dat[,paste0("Spd.",stages[i+1])]), #Spawning day
      c(dat[,paste0("Duration.",stages[i])],dat[,paste0("Duration.",stages[i+1])])) #Stage duration
  colnames(Data)<-c("Year","Day","Lon","Lat","Abundance","Age","Spd","Dur")  
  Data<-Data[!is.na(Data$Abundance) & Data$Abundance>0,] #Remove zeroes or missing values
  Data$Year<-as.factor(Data$Year) #Make sure year is a factor
  Data$Abundance<-Data$Abundance/Data$Dur #Scale abundance by stage duration
  Data$Age<-scale(Data$Age,center=TRUE,scale=FALSE) #Center age around zero
  Data$logAbundance<-log(Data$Abundance) #Log transform abundance
  m <- gam(Data$logAbundance~s(Spd)+te(Lon,Lat,k=5)+s(Year,bs="re")+Age+s(Year,bs="re",by=Age),data=Data) #SRA model
  mortalities[i]<- -coef(m)[2] #The negative of the age-coefficient gives mortality
}
