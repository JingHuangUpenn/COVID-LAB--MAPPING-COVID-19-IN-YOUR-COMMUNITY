
library(ggplot2)
library(EpiEstim)
library(nlme)
library(dlnm)
library(dplyr)
library(tsModel)
library(sjstats)
library(corrplot)
library(MASS)
library(tidyr)

# set working directory
setwd('C:/Users/Jing Huang/Dropbox/000_Collaboration/00_CHOP_General/Covid/data_analysis/st_prediction')
setwd('/Users/jing14/Dropbox/000_Collaboration/00_CHOP_General/Covid/data_analysis/st_prediction')
# source R function
source('./functions/covid-function.R')


## read case data
mydat <- read.csv("./rawdata/factnew_2020-08-11.csv")
length(unique(mydat$fips))#745
mydat$date <- as.Date(mydat$date,"%Y-%m-%d")
names(mydat)[names(mydat)=='STATE_NAME']='state'
names(mydat)[names(mydat)=='COUNTY_NAME']='county'

mydat <- mydat[,c('fips','date','county','State','state','cases')]
mydat <- mydat[!(is.na(mydat$cases)),]
mydat <- mydat[mydat$cases>0,]


mydat <- arrange(mydat,fips,date)
summary(mydat$date)#3/1 to 8/6
mydat <- distinct(mydat)

## testing data
n.test <- read.csv('./rawdata/states-daily20200327.csv')

my.max <- mydat%>% group_by(fips)%>%
  summarize(county = unique(county), state= unique(state),
            state2 = unique(State), max_cases = max(cases,na.rm = TRUE),
            ob_days= sum(cases>5))

my.max <- my.max%>%arrange(desc(max_cases))
my.county <- my.max$county
length(my.county)#745

### specify assumptions
# set 1
# Assumption 1: generation time
GT.ncov  <- discr_si(seq(0, 30), 7.5, 3.4)
# Assumptions 2: sliding window, 3-5 days
nd <- 3
# Assumptions 3: prior mean of R = 2.5
meanprior <- 2.5


#estimate R 
gen.r(nd = nd, meanprior = meanprior, GT.ncov= GT.ncov)



# set the basic 569 counties to train the model
mainfips <- read.csv("./rawdata/county569train_0728.csv")
mainfips <- unique(mainfips$fips)
length(mainfips)#569

nbdat <- read.csv('./rawdata/Near_571_updated.csv')

myR <- read.csv('./derived/estimated_R_2020-08-12.csv')
myR$date<-as.Date(myR$date,'%Y-%m-%d')

reg.dat <- read.csv('./rawdata/temp_casetest_soc_08_10.csv')
reg.dat$date <- as.Date(reg.dat$date,"%m/%d/%Y")

reg.dat <- distinct(reg.dat)

### merge R with gis, temp
reg.dat <- left_join(reg.dat,myR)

czone <- read.csv('./rawdata/DOE_ClimateZones_updated.csv')
names(czone)[names(czone)=='FIPS']='fips'
czone <- czone[,-1]
head(czone)

reg.dat <- left_join(reg.dat,czone,by=c("fips"))

#drop Puerto Rico
reg.dat <- reg.dat[reg.dat$fips!=72127,]
table((distinct(reg.dat[,c("fips","BA_Climate_Zone")]))$BA_Climate_Zone)
reg.dat$zones <- as.factor(reg.dat$BA_Climate_Zone)
class(reg.dat$zones)#factor

levels(reg.dat$zones) <- c("Cold/Very Cold","Hot-Dry/Mixed-Dry",
                           "Hot-Humid/Mixed-Humid","Marine","Hot-Dry/Mixed-Dry",
                           "Hot-Humid/Mixed-Humid","Subarctic","Cold/Very Cold" )

levels(reg.dat$zones) <- c("Cold/Very Cold","Hot-Dry/Mixed-Dry",
                           "Hot-Humid/Mixed-Humid","Marine","Hot-Dry/Mixed-Dry",
                           "Hot-Humid/Mixed-Humid","Cold/Very Cold" )

# standardize county level covariates
tt.dat <- distinct(reg.dat[,c("fips","DIABETES","OBESE","PopDens",
                              "PovRatioLT2","Uninsured","Crowded","SMOKE")])
tt.dat <- tt.dat[complete.cases(tt.dat),]
tt.dat$popdens_s <- scale(log(tt.dat$PopDens))

reg.dat <- left_join(reg.dat,tt.dat)

reg.dat <- reg.dat %>% dplyr::group_by(fips) %>% tidyr::fill(R)

## generate variable y = log(R)
reg.dat$y <- log(reg.dat$R)


## 
reg.dat$diabetes_s <- reg.dat$DIABETES/100
reg.dat$obese_s <- reg.dat$OBESE/100
reg.dat$smoke_s <- reg.dat$SMOKE/100



## bound visitation 
reg.dat$daily_visitation_diff[reg.dat$daily_visitation_diff>0.1]=0.1

## generate lag socialing distancing measure

reg.dat$visit_mean<- runMean(reg.dat$daily_visitation_diff,lags=4:14,group=reg.dat$fips)
reg.dat$visit_mean3<- runMean(reg.dat$daily_visitation_diff,lags=0:2,group=reg.dat$fips)


## temp
reg.dat$dbt_mean<- runMean(reg.dat$dry_bulb_temp,lags=4:14,group=reg.dat$fips)
reg.dat$hum_mean<- runMean(reg.dat$abs_humid_est,lags=4:14,group=reg.dat$fips)

## center temp at 53 and dtemp at 56
#first, imput random missing temp
reg.dat<-reg.dat%>%
  group_by(fips) %>%
  mutate(lag.wet = dplyr::lag(wet_bulb_temp_F, n = 1, default = NA))
#View(reg.dat[reg.dat$fips==42101,c("fips","date","wet_bulb_temp_F","lag.wet")])
#reg.dat$wet_bulb_temp_F<- ifelse(reg.dat$date==as.Date("2020-07-31",'%Y-%m-%d')&
#                                   is.na(reg.dat$wet_bulb_temp_F), reg.dat$lag.wet, reg.dat$wet_bulb_temp_F)


reg.dat$temp <- reg.dat$wet_bulb_temp_F-mean(reg.dat$wet_bulb_temp_F,na.rm = T)
reg.dat$dtemp <- reg.dat$dry_bulb_temp-mean(reg.dat$dry_bulb_temp,na.rm = T)


## create lag y for AR1
reg.dat<-reg.dat%>%
  group_by(fips) %>%
  mutate(lag.y = dplyr::lag(y, n = 1, default = NA),
         lag.R = dplyr::lag(R, n = 1, default = NA))



# remove rows will all variables missing
reg.dat <- reg.dat[rowSums(is.na(reg.dat)) != ncol(reg.dat), ]

#View(reg.dat)
#reg.dat <- reg.dat[reg.dat$y>=]
reg.dat$Negative[is.na(reg.dat$Negative)]=0
reg.dat$Positive[is.na(reg.dat$Positive)]=0

####  for county level testing all in 08-10
#Total_test
#Positive_county
#Negative_county

reg.dat$Negative_county[is.na(reg.dat$Negative_county)]=0
reg.dat$Positive_county[is.na(reg.dat$Positive_county)]=0
reg.dat$Total_test[is.na(reg.dat$Total_test)]=0

reg.dat$Total_test[reg.dat$fips==44000] <- reg.dat$Negative[reg.dat$fips==44000] + reg.dat$Positive[reg.dat$fips==44000] 

### testing
#reg.dat$tottest <- reg.dat$Negative+reg.dat$Positive
reg.dat<-reg.dat%>%
  group_by(fips) %>%
  mutate(dailytest = c(first(Total_test),diff(Total_test)),
         dailypos = c(first(Positive_county),diff(Positive_county)))

reg.dat$dailytest[reg.dat$dailytest<0]=0
reg.dat$dailypos[reg.dat$dailypos<0]=0

reg.dat$positive7 <- runMean(reg.dat$dailypos,lags=0:6,group=reg.dat$fips)*7
reg.dat$daily7 <- runMean(reg.dat$dailytest,lags=0:6,group=reg.dat$fips)*7
reg.dat$posrate7 <- reg.dat$positive7/reg.dat$daily7
reg.dat$posrate <- reg.dat$dailypos/reg.dat$dailytest

reg.dat$positive3 <- runMean(reg.dat$dailypos,lags=0:2,group=reg.dat$fips)*3
reg.dat$daily3 <- runMean(reg.dat$dailytest,lags=0:2,group=reg.dat$fips)*3
reg.dat$posrate3 <- reg.dat$positive3/reg.dat$daily3

reg.dat$positive5 <- runMean(reg.dat$dailypos,lags=0:4,group=reg.dat$fips)*5
reg.dat$daily5 <- runMean(reg.dat$dailytest,lags=0:4,group=reg.dat$fips)*5
reg.dat$posrate5 <- reg.dat$positive5/reg.dat$daily5

reg.dat<-reg.dat%>%
  group_by(fips) %>%
  mutate(posrate_diff7 = c(NA,diff(posrate7)),
         lag.posrate7 = dplyr::lag(posrate7, n = 1, default = NA),
         posdiff_rel7 = (posrate_diff7)/(lag.posrate7),
         posrate_diff5 = c(NA,diff(posrate5)),
         lag.posrate5 = dplyr::lag(posrate5, n = 1, default = NA),
         posdiff_rel5 = (posrate_diff5)/(lag.posrate5),
         posrate_diff3 = c(NA,diff(posrate3)),
         lag.posrate3 = dplyr::lag(posrate3, n = 1, default = NA),
         posdiff_rel3 = (posrate_diff3)/(lag.posrate3),
         posrate_diff = c(NA,diff(posrate)),
         lag.posrate = dplyr::lag(posrate, n = 1, default = NA),
         posdiff_rel = (posrate_diff)/(lag.posrate)
  )



# adjust cap, change in 06/29, 
reg.dat$posdiff_rel7[reg.dat$posdiff_rel7>=2]= max(reg.dat$posdiff_rel7[reg.dat$posdiff_rel7<=2],na.rm = T)
reg.dat$posdiff_rel5[reg.dat$posdiff_rel5>=2]= max(reg.dat$posdiff_rel5[reg.dat$posdiff_rel5<=2],na.rm = T)
reg.dat$posdiff_rel3[reg.dat$posdiff_rel3>=2]= max(reg.dat$posdiff_rel3[reg.dat$posdiff_rel3<=2],na.rm = T)
reg.dat$posdiff_rel[reg.dat$posdiff_rel>=20]= max(reg.dat$posdiff_rel[reg.dat$posdiff_rel<=20],na.rm = T)



#rescale humidity data
reg.dat$hum_mean_logit <- logit(reg.dat$hum_mean)


# case prevalence
reg.dat$casepct <- reg.dat$cases/reg.dat$TotPop


# case percentage
reg.dat<-reg.dat%>%
  group_by(fips) %>%
  mutate(lag.casepct7 = dplyr::lag(casepct, n = 7, default = NA))

reg.dat<-reg.dat%>%
  group_by(fips) %>%
  mutate(statemasking = dplyr::lag(statemask, n = 7, default = NA))


## fill visit_mean and testing positivity
reg.dat <- reg.dat %>% dplyr::group_by(fips) %>% tidyr::fill(visit_mean)
reg.dat <- reg.dat %>% dplyr::group_by(fips) %>% tidyr::fill(visit_mean3)
#summary(reg.dat$posdiff_rel7)
reg.dat$posdiff_rel7[is.nan(reg.dat$posdiff_rel7)]=NA
reg.dat$posrate7[is.nan(reg.dat$posrate7)]=NA

#View(reg.dat[,c("fips","date","posdiff_rel7")])

reg.dat <- reg.dat %>% dplyr::group_by(fips) %>% tidyr::fill(posdiff_rel7)
reg.dat <- reg.dat %>% dplyr::group_by(fips) %>% tidyr::fill(lag.casepct7)
reg.dat <- reg.dat %>% dplyr::group_by(fips) %>% tidyr::fill(posrate7)


nkd=2
vk <- c(33,55,67)-53
vkd <- c(50,60,70,80)-56
vkh <- c(0.0057,0.0087,0.013)

lk <- equalknots(4:14, nk=nkd)
cb1 <- crossbasis(reg.dat$temp, lag=c(4,14), group=reg.dat$fips,
                  argvar=list(fun="ns", knots=vk), 
                  arglag=list(fun="ns", knots=lk))

data <- lm(y~ fips+date+ R+lag.R+lag.y+cb1+cb2 +visit_mean+ visit_mean3+  popdens_s+hum_mean
           +posdiff_rel7+zones+posrate7+ PopDens+COUNTY_NAME+State+ lag.casepct7+
             diabetes_s +PopO64+Crowded+TotPop, data=reg.dat,
           method="model.frame")

newfips <- unique(data$fips[!data$fips%in%mainfips])
length(mainfips)#571  --- 569
length(newfips)#175
imp.fips1 <- distinct(data[data$fips%in%mainfips,c("fips","fips")])
names(imp.fips1) <- c('fips','imp_fips')
dim(imp.fips1)
head(imp.fips1)

imp.fips2 <- nbdat[nbdat$FIPS%in%newfips,c("FIPS","NEAR_FIPS")]
dim(imp.fips2)
names(imp.fips2) <- names(imp.fips1) 
class(imp.fips2$imp_fips) <- class(imp.fips1$imp_fips)
head(imp.fips2)

imp.fips <- bind_rows(imp.fips1,imp.fips2)

imp.zones <- distinct(reg.dat[reg.dat$fips%in%mainfips,c("fips","zones")])
names(imp.zones) <- c('imp_fips','imp_zones')

imp.fips <- left_join(imp.fips,imp.zones)

data <- left_join(data,imp.fips)

data$number = 1
data <- data %>%
  group_by(fips) %>%
  mutate(ticker = cumsum(number))%>%
  mutate(ticker_s = ticker)

reg.dat$date <- as.Date(reg.dat$date)
casedt <- max(reg.dat$date[!is.na(reg.dat$cases)])

temp.pred <- data[data$date>casedt & data$date<=(casedt+28),]
traindata <- data[data$date<=casedt,]


## set y <-3 to NA in training 
traindata$y[(!is.na(traindata$y) & traindata$y< -3)]=NA
traindata$lag.y[(!is.na(traindata$lag.y) & traindata$lag.y< -3)]=NA

traindata <- traindata[traindata$ticker>7,]
case.dat <- read.csv("./rawdata/factnew_2020-08-11.csv")
case.dat$date <- as.Date(as.Date(case.dat$date,"%Y-%m-%d"),"%Y-%m-%d")
head(case.dat$date)

case.dat$cases[!is.na(case.dat$cases_nyt)] =case.dat$cases_nyt[!is.na(case.dat$cases_nyt)]

case.dat <- case.dat[c('fips','date','cases')]

case.dat<- distinct(case.dat)
case.dat <- bind_rows(reg.dat[reg.dat$date>max(case.dat$date),c("fips","date","cases")],case.dat)
case.dat <- arrange(case.dat,fips,date)
### change future visit_mean to last 7 day average
### change future positivity to median
medianposrate7 <- traindata%>%group_by(fips)%>%filter(date >= casedt-14)%>%summarise(median_posrate7 = median(posrate7,na.rm = T),
                                                                                     median_posdiff_rel7 = median(posdiff_rel7,na.rm = T))

testdata <- temp.pred
testdata$number = 1
testdata <- testdata %>%
  group_by(fips) %>%
  mutate(ticker = cumsum(number))%>%
  mutate(ticker_s = ticker)

testdata <- left_join(testdata,medianposrate7)

testdata$posdiff_rel7[testdata$ticker>=3]= testdata$median_posdiff_rel7[testdata$ticker>=3]
testdata$posrate7[testdata$ticker>=3]= testdata$median_posrate7[testdata$ticker>=3]

testdata$posrate7[testdata$posrate7==Inf]=testdata$median_posrate7[testdata$posrate7==Inf]
traindata2 <- traindata[traindata$fips%in%mainfips,]

min(traindata2$date);max(traindata2$date)
length(unique(traindata2$fips))

myfips <- unique(testdata$fips)
length(myfips)
length(mainfips)

f2 <- lme(y~  lag.y+cb1 +visit_mean +  popdens_s+
            +posdiff_rel7+posrate7+ 
            diabetes_s +PopO64+Crowded,method="ML",na.action=na.omit,
          random=~1+lag.y|imp_zones/imp_fips, control = lmeControl(opt = 'optim'),
          data = traindata2)

myst <- max(traindata$date)
testdata$date <- as.Date(testdata$date)
nd = 3; pdays =c(max(testdata$date)-min(testdata$date))+1;
GT.ncov  <- discr_si(seq(0, 30), 7.5, 3.4)

t0<-Sys.time()
pred.r.case.st(nd = nd, pdays = pdays,fmodel = f2,myst = myst,
               temp.pred = testdata,case.dat = case.dat, myfips = myfips,
               name='f2')
Sys.time()-t0

