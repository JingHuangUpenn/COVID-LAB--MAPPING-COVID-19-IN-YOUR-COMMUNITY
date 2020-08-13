expit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}
gen.r <- function(nd = 5, meanprior = 2.5, GT.ncov = discr_si(seq(0, 30), 6, 3.8)){
  # vector to save counties without R results yet
  sk.ct <- data.frame(id = double(),county = character(),state = character())
  
  for(i in c(1:length(my.county)) ){
    fips = my.max$fips[i]
    jj = as.character(my.max$county[i])
    if (jj =="Unknown"){next}
    zz = as.character(my.max$state2[i] )
    if(is.na(zz)){next}
    zz.full = as.character(my.max$state[i])
    jj.fips = my.max$fips[i]
    new.n.test <- n.test[n.test$state==zz,c('date','total','positive','negative','hospitalized')]
    new.n.test$date <- as.Date(as.character(new.n.test$date),"%Y%m%d")
    
    new.n.test <- new.n.test[order(new.n.test$date),]
    
    pos.rate = diff(new.n.test$positive)/diff(new.n.test$total)
    
    
    # find the last date that need to smooth
    # criteria: the most recent day that is earlier than March 27, and with positive rate <= 0.1 or >=0.9
    ind.d <- (new.n.test$date[-1])[max( which((pos.rate>=0.9| pos.rate<=0.05) & (new.n.test$date[-1]) <= as.Date("2020-03-27","%Y-%m-%d") ) )]+1
    
    ## smooth the incidence and impute the earlier days
    tmp.dat <- mydat[mydat$fips==jj.fips,c('date','cases')]

    if(dim(tmp.dat)[1]<1){
      sk.ct.tt <- data.frame(id = i,county = jj,state = zz);
      sk.ct <- rbind(sk.ct, sk.ct.tt)
      next}
    tmp.dat$date <- as.Date(tmp.dat$date,'%m/%d/%Y')
    tmp.dat$I <- c(tmp.dat$cases[1],diff(tmp.dat$cases))
    
    
    #### remove the previous criteria if(mean(tmp.dat$I)<5 | sum(tmp.dat$I >=5)<3 ) 
    
    
    ## kernal smoothing
    if (!is.na(ind.d) & ind.d >=  min(tmp.dat$date)){
      ks <- ksmooth(c(1:length(tmp.dat$I)),
                    tmp.dat$I, "box", bandwidth = 7)
      
      tt <- unlist(lapply(c(1:length(tmp.dat$I)), function(x){which.min(abs(x-ks$x)) }))
      
      new.i = c((round(ks$y[tt]))[which(tmp.dat$date <= ind.d)],tmp.dat$I[c((max(which(tmp.dat$date <= ind.d))+1):length(tmp.dat$I))])
      new.time = tmp.dat$date
      new.n = cumsum(new.i)
      
    }else{
      new.i = tmp.dat$I
      new.time = tmp.dat$date
      new.n = cumsum(new.i)
    }
    
    
    t_start <- seq(max(2,min(which(new.n>=5))), length(new.i)-nd+1)
    t_end <- t_start+nd-1
    i.dat <- data.frame(I = new.i, dates = new.time)
    i.dat$I[i.dat$I<0]=0
    i.dat.tt <- data.frame(dates = seq.Date(min(i.dat$dates),max(i.dat$dates),1), 
                           I.add = rep(0, length(seq.Date(min(i.dat$dates),max(i.dat$dates),1)) ))
    i.dat <- left_join(i.dat.tt,i.dat)
    i.dat[is.na(i.dat)]=0
    i.dat$I <- i.dat$I.add + i.dat$I
    i.dat <- i.dat[,-which(colnames(i.dat)=='I.add')]
    res <- estimate_R(incid = i.dat,
                      method = "non_parametric_si",
                      config = make_config(list(si_distr = GT.ncov,
                                                t_start = t_start,
                                                t_end = t_end,
                                                mean_prior=meanprior)))
    tmp.r.est<- data.frame(
      fips = rep(fips,length(res$R$`Median(R)`)),
      date= res$dates[res$R$t_start],
      R = res$R$`Median(R)`,
      R_lower = res$R$`Quantile.0.025(R)`,
      R_upper = res$R$`Quantile.0.975(R)`,
      R_std = res$R$`Std(R)`)
    
    if(i ==1){
      r.est <- tmp.r.est
    }else{
      r.est <- rbind(r.est,
                     tmp.r.est)
    }
    
  }
  write.csv(r.est,paste('./estimated_R_',Sys.Date(),'.csv',sep=""),row.names=F)
  write.csv(sk.ct,paste('./skipped_counties_',Sys.Date(),'.csv',sep=""),row.names=F)
}


pred.r.case.st <- function(nd = nd, pdays = 7,
                           fmodel = fmodel,
                           temp.pred = temp.pred, 
                           case.dat = case.dat,
                           myfips= myfips, myst=myst,
                           name = 'dry'){
  
  zzpred <-temp.pred
  
  zzpred$y_est<- zzpred$y
  zzpred$lag.y_est <- zzpred$lag.y
  
  zzpred <- zzpred[zzpred$date>myst,]
  
  for(i in c(1:length(myfips))){
    fips <- myfips[i]
    
    
    newdat <-zzpred[zzpred$fips==fips,]
    
    if(dim(newdat)[1]==0){#print(mystation[i,]);
      next}
    
    
    
    mycase <- case.dat[case.dat$fips==fips,]
    
    mycase$inc <- c(mycase$cases[1],diff(mycase$cases))
    
    tt.i <- mycase[mycase$date<=myst,c('date','inc')]
    tt.i$imp <- rep(0,dim(tt.i)[1])
    tt.i$cases = cumsum(tt.i$inc)
    
    
    for(j in 1:(dim(newdat)[1])){   
      
      newdat$y[j] <-c(predict(fmodel,newdat[j,]))
      
      
      if(j != dim(newdat)[1]){
        newdat$lag.y[j+1] <- newdat$y[j]  
      }
      
      newdat$R[j]<- exp(newdat$y[j])
      
      
      lam <- sum((GT.ncov)[1:(min(length(tt.i$inc),31))]*rev(tt.i$inc)[1:(min(length(tt.i$inc),31))],na.rm = T)
      
      xx <- lam*newdat$R[j]
      
      tt.i<- rbind(tt.i,
                   data.frame(date= tail(tt.i$date,1)+1,
                              inc = xx,imp=1,
                              cases = tail(tt.i$cases,1)+xx
                   ))
      if(j != dim(newdat)[1]){
        newdat$lag.casepct7[j+1] <- tail(tt.i$cases,7)[1]/newdat$TotPop[j+1]
        newdat$lag.casepct14[j+1] <- tail(tt.i$cases,14)[1]/newdat$TotPop[j+1]
        newdat$casepct[j+1] <- tail(tt.i$cases,1)[1]/newdat$TotPop[j+1]
      }
    }
    
    
    
    if(i ==1){
      all.pred <- data.frame(fips = rep(fips,dim(newdat)[1]), date = newdat$date,
                             obs_case =mycase$cases[mycase$date%in%newdat$date],
                             obs_inc = mycase$inc[mycase$date%in%newdat$date],
                             pred_R = newdat$R,est_R = exp(newdat$y_est),
                             pred_inc = tt.i$inc[tt.i$imp==1],
                             pred_case = tt.i$cases[tt.i$imp==1]#,
                            
      )
      
    }else{
      all.pred <- rbind(all.pred,
                        data.frame(fips = rep(fips,dim(newdat)[1]), date = newdat$date,
                                   obs_case =mycase$cases[mycase$date%in%newdat$date],
                                   obs_inc = mycase$inc[mycase$date%in%newdat$date],
                                   pred_R = newdat$R,est_R = exp(newdat$y_est),
                                   pred_inc = tt.i$inc[tt.i$imp==1],
                                   pred_case = tt.i$cases[tt.i$imp==1]#,
                                  
                        )
      )
    }
    
  }
  
  write.csv(all.pred,paste('./pred_st_',name,'_',Sys.Date(),'.csv',sep=''),row.names = F)
  
}



