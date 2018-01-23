dynamicIntegrateAUC <- function(survival.time, survival.status, start=NULL, marker, cutoffTime) {
   if(is.null(start))
      kmfit <- survfit(Surv(time=survival.time, event=survival.status) ~ 1)
   else
      kmfit <- survfit(Surv(time=start, time2=survival.time, event=survival.status) ~ 1)
      
   mmm <- MeanRank( survival.time= survival.time, survival.status= survival.status, marker= marker, start=start )

   #Get overlap between survival function and mmm
   meanRanks <-  mmm$mean.rank[which(mmm$time <= cutoffTime)]
   survTimes <- mmm$time[mmm$time <= cutoffTime]
   timeMatch <- match(survTimes, kmfit$time)

   S_t <- kmfit$surv[timeMatch]

   #Calculate weights for c-index
   f_t <- c( 0, (S_t[-length(S_t)] - S_t[-1]) )
   S_tao <- S_t[length(S_t)]
   weights <- (2*f_t*S_t)/(1-S_tao^2)
   
   return( sum(meanRanks * weights) ) #C-index
}

