library(doSNOW)
library(foreach)
library(tidyverse)

rm(list=ls())
setwd("/home/alison/Dropbox/03_IPAC/R_code/VisitsBasedModel")

SEED <-1
nruns <- 200 #was 200
sim_length <- 10 #was 730
warmup <-100 #was 1000

#####input variables#####
wait_target <- 0 # was 1 - days waiting target before service starts
PerTar <- 1.0 #percentage to be serviced within target days
lambda <- 12.6 #P1 arrival rate non-covid uplift *adjusted*

#set the minimum and maximum LOS and ISR. It will simulate required capacity for all these combinations.
ISR <- 4
endSR <- 2
avg_LOS <- 12.2
weekly_patients <- lambda *7 #number of arrivals per week

cost_acute <- 346
cost_community <- 125

add<-0 #initialising number of slots to add to capacity

###initialising Final df###

rows<-(avg_LOS+1)*4 #defining how many rows are needed for Final data frame
u<- 0 # initialising counter for final data frame
Final <- data.frame(LOS = integer(rows),
                    ISR = integer(rows),
                    nruns = integer(rows),
                    sim_length = integer(rows),
                    warm_up=integer(rows),
                    capacity = integer(rows),
                    added_cap = integer(rows),
                    target = integer(rows),
                    percentageTar=numeric(rows),
                    sd_Tar = numeric(rows),
                    min_Tar = numeric(rows),
                    max_Tar =numeric(rows),
                    minCI_Tar = numeric(rows),
                    plusCI_Tar = numeric(rows),
                    avg_wait = numeric(rows),
                    avg_delay_cost = numeric(rows),
                    minCI_cost = numeric(rows),
                    maxCI_cost = numeric(rows),
                    sd_wait = numeric(rows),
                    max_wait = numeric(rows),
                    minCI_wait = numeric(rows),
                    plusCI_wait = numeric(rows),
                    avg_Q = numeric(rows),
                    sd_Q = numeric(rows),
                    min_Q = numeric(rows),
                    max_Q = numeric(rows),
                    minCI_Q = numeric(rows),
                    plusCI_Q = numeric(rows),
                    res_used= numeric(rows),
                    mean_idle= numeric(rows),
                    sd_idle = numeric(rows),
                    min_idle = numeric(rows),
                    max_idle = numeric(rows),
                    minCI_idle = numeric(rows),
                    plusCI_idle = numeric(rows),
                    in_sys_cap = numeric(rows),
                    stringsAsFactors=FALSE)

#loop through range of capacities
  while (TRUE){
    
    if (add==0){ #only calculate capacity if its the first iteration for the parameter combo
      
      if(u>0){ 
        same <- which(Final$LOS==avg_LOS) #if the ISR is different, then take last capacity from when same LOS was last run
        if (Final[u,2]==ISR) { #if just starting or if ISR is same as ISR in previous iteration
          n_slots <- Final$capacity[u] #start capacity from previous parameter combination}
        } else if (length(same)>0){
          n_slots <- max(Final$capacity[same]) #start from maximum simulated capacity from when same LOS was last run
        } else {n_slots <- 480} #ceiling(lambda * avg_LOS* (ISR/2 + endSR/2)+ISR)} #start from calculate minimum capacity if else
      } else {  
        n_slots <-480 # ceiling(lambda * avg_LOS* (ISR/2 + endSR/2)+ISR) #start from calculate minimum capacity if else
      }
    }

    ####runs####
    cl<-makeCluster(4)
    registerDoSNOW(cl)
    RESULTS<-foreach(run=1:nruns,.combine="rbind") %dopar% {
      set.seed(nruns*(SEED-1)+run)
      
      #distributions
      #arrival distribution
      dis_arrival <- function(){
        x<- round(rpois(1,lambda=lambda))
        if (x<=0){ #if x=0 or negative, then x = 0, no one arrives that day
          x <-0
          return(x)
        } else {return(x)}
      }
      
      #LOS distribution
      dis_los <- function(){
        x<- round(rnorm(1,mean=avg_LOS,sd=2))
        if (x<=0){ #if x=0 or negative, then x = 1, los needs to be at least one day
          x <-1
          return(x)
        } else {return(x)}
      }
      
      #ISR distribution
      dis_init_slots <- function(){
        x<- round(rnorm(1,mean=ISR,sd=1))
        if (x<=0){ #if x=0 or negative, then x = 1, there needs to be at least one visit per day to start with
          x <-1 
          return(x)
        } else if (x>6){
          x<-6
          return(x)
        } else if (x > n_slots){ #ISR can not be greater then the total number of slots available in one day
          x <-n_slots
          return(x)
        } else {return(x)}
      }
      
      #end SR distribution 
      dis_end_slots <- function(){
        x<- round(rnorm(1,mean=endSR,sd=1))
        if (x<=0){ #if x=0 or negative, then x = 1, there needs to be at least one visit per day to start with
          x <-1 
          return(x)
        } else if (x>ISR){
          x<-ISR
          return(x)
        } else if (x > n_slots){ #ISR can not be greater then the total number of slots available in one day
          x <-n_slots  #unsure if this is needed here
          return(x)
        } else {return(x)}
      }
      
      #####output variables#####
      ent_sys <- 0 # number of entities that entered the system
      left_sys <-0 # number of entities that left the system
      
      #output after warm up period
      output<-data.frame(RUNX=integer(sim_length), #run number x
                         day= integer(sim_length), #output per day
                         q_length = integer(sim_length), #number of patients in the queue
                         res_used=numeric(sim_length), #used slots
                         res_idle=numeric(sim_length), #idle slots
                         in_sys=numeric(sim_length) #number of patinets in the system
      ) 
      
      
      #####creating necessary data structures#####
      #patient list
      patients<-data.frame(id=integer((sim_length+warmup)*2),            #patient id
                           los=integer((sim_length+warmup)*2),           #length of stay
                           arrival_time =integer((sim_length+warmup)*2), # day in the simulation the entity arrived
                           start_service=integer((sim_length+warmup)*2), # day actual service started
                           end_service=integer((sim_length+warmup)*2),   # day service ended
                           wait_time=integer((sim_length+warmup)*2),     # number of days spent in the queue
                           exit=logical((sim_length+warmup)*2),          # boolean variable, TRUE if the entity has left the system
                           stringsAsFactors=FALSE)
      npat<-0 #initialising counter for patients dataframe
      
      #list with required visit vectors for each patient
      req_visits <- list()
      
      #resources
      resources <- matrix(nrow=(sim_length+warmup)*2, ncol = 1) #times 2 to make the calculations for resources work at the end of the simulation
      
      resources[,] <- n_slots
      
      #vector for storing waiting time, kept for each patient who left the system
      waittime_vec <- data.frame(RUNX=integer(),
                                 start_service= integer(),
                                 waittime = integer(),
                                 stringsAsFactors=FALSE)
      
      #####simulation#####
      id<-0  

      for (t in 1:(sim_length+warmup)) {
        #arrivals to service
        narr<-dis_arrival()
        if(narr>0){
          ent_sys <- ent_sys + narr
          
          #for each arrived patient
          for (j in 1:narr) {
            id<-id+1
            npat<-npat+1
            los<- dis_los()
            arrival_time <- t
            exit <-FALSE
            patients[npat, ] <- c(id,los,arrival_time,NA, NA, 0, exit)
            
            #initial slots and creating required visits vector
            init_slots <- dis_init_slots()
            end_slots <- dis_end_slots()
            visit_vector <- round(seq(init_slots,end_slots,length.out = los))
            
            req_visits[[id]] <- visit_vector
            
            #planning service, checking resources
            tt<-t #temporary t for incrementing when no resources available
  
            while (is.na(patients$start_service[npat])==TRUE){
              if (all((resources[((tt):((tt)+patients$los[npat]-1)),]>= req_visits[[id]])==TRUE)){
                patients$start_service[npat] <- tt
                patients$end_service[npat] <- patients$start_service[npat]+(patients$los[npat]-1)
                
                #decrease capacity
                resources[((tt):((tt)+patients$los[npat]-1)),] <- resources[((tt):((tt)+patients$los[npat]-1)),] - req_visits[[id]]
              } else {
                tt<-tt+1 #if no sufficient resources, check for starting on the next day
              } 
            }          
          }
        }
        
        #increase wait time for patients in the queue
        in_q<-which((patients$start_service>t)&(patients$id>0))
        if (length(in_q)>0){
          patients[in_q,6]<- patients[in_q,6]+1
        }
        
        
        #recording output from the day warm up period has finished
        if (t>warmup){ #only start recording after the warm up period
          if (npat>0 & nrow(waittime_vec)>0) {
            output[t-warmup, ]<- c(RUNX=run,
                                   day= t,
                                   q_length = length(in_q),
                                   res_used= 1- (resources[t,]/n_slots),
                                   res_idle= resources[t,]/n_slots,
                                   in_sys = (ent_sys - left_sys))
          } else if (npat>0 & nrow(waittime_vec)==0) {
            output[t-warmup, ]<- c(RUNX=run,
                                   day= t,
                                   q_length = length(in_q),
                                   res_used= 1- (resources[t,]/n_slots),
                                   res_idle= resources[t,]/n_slots,
                                   in_sys = (ent_sys - left_sys))
          } else {
            output[t-warmup, ]<- c(RUNX=run,
                                   day= t,
                                   q_length = length(in_q),
                                   res_used= 1- (resources[t,]/n_slots),
                                   res_idle= resources[t,]/n_slots,
                                   in_sys = (ent_sys - left_sys))
          }
        }
        
        #remove patients whose service has ended from the patients table
        remove <- which(patients$end_service==t)
        if(length(remove)>0){
          if(t>=warmup){
            df<-data.frame(RUNX = run, start_service= patients$start_service[remove], waittime= patients[remove,6])
            waittime_vec <- rbind(waittime_vec,df) #keeping waiting time
          }  
          patients <- patients[-remove,] #remove from patient list
          npat<- npat - length(remove)
          left_sys <- left_sys + length(remove)
        }
        
      }
      
      list<-list(output, resources, waittime_vec, req_visits)
      
      return(list)
      
    }
    
    stopCluster(cl)
    ################################################

    #creating dataframe for summary info
    summary <- data.frame(LOS = integer(nruns),
                          ISR = integer(nruns),
                          nruns = integer(nruns),
                          sim_length = integer(nruns),
                          warm_up=integer(nruns),
                          capacity = integer(nruns),
                          added_cap = integer(nruns),
                          target = integer(nruns),
                          percentage1=numeric(nruns),
                          mean_wait= numeric(nruns),
                          mean_cost= numeric(nruns),
                          min_CI_cost = numeric(nruns),
                          max_CI_cost = numeric(nruns),
                          q_length = numeric(nruns),
                          res_used= numeric(nruns),
                          res_idle= numeric(nruns),
                          in_sys = numeric(nruns))
    
    #splitting up RESULTS list in 2
    output<-RESULTS[,1]
    out<-do.call(rbind, output) #combining in one dataframe
    
    resources<-RESULTS[,2]
    res<-do.call(cbind, resources) 
    colnames(res)<- c(1:nruns)
    
    waittimes <- RESULTS[,3]
    wait<-do.call(rbind, waittimes)

    req_visits <- RESULTS[,4]
    visits <- do.call(rbind, req_visits)
    #summary of all runs
    for (k in 1:nruns){
      r.out <- which(out[,1]==k)
      k.wait <- which(wait[,1]==k)
      mean_acute_cost <- (cost_acute*(mean(wait$waittime[k.wait]))*weekly_patients)
      minCI_acute <- (cost_acute*(quantile(wait$waittime[k.wait], 0.05))*weekly_patients)
      maxCI_acute <- (cost_acute*(quantile(wait$waittime[k.wait], 0.95))*weekly_patients)
      mean_comm_cost <- (cost_community*n_slots)
      summary[k,]<- c(LOS = avg_LOS,
                      ISR = ISR,
                      nruns = nruns,
                      sim_length = sim_length,
                      warm_up=warmup,
                      capacity = n_slots,
                      added_cap = add,
                      target = wait_target,
                      percentage1=mean(wait$waittime[k.wait] <= wait_target),
                      mean_wait= round(mean(wait$waittime[k.wait]),2),
                      mean_cost= round(mean_acute_cost+mean_comm_cost,2),
                      min_CI_cost = round(minCI_acute+mean_comm_cost, 2),
                      max_CI_cost = round(maxCI_acute+mean_comm_cost, 2),
                      q_length = round(mean(out$q_length),2),
                      res_used= round(mean(out$res_used[r.out]),2),
                      res_idle= round(mean(out$res_idle[r.out]),2),
                      in_sys= round(mean(out$in_sys[r.out]),2) ) 
    }
    
    head(summary)

u<-u+1
    Final[u,] <- c(LOS = avg_LOS,
                   ISR = ISR,
                   nruns = nruns,
                   sim_length = sim_length,
                   warm_up=warmup,
                   capacity = n_slots,
                   added_cap = add,
                   target = wait_target,
                   percentageTar=mean(summary$percentage1),
                   sd_Tar = sd(summary$percentage1),
                   min_Tar = min(summary$percentage1),
                   max_Tar = max(summary$percentage1),
                   minCI_Tar = quantile(summary$percentage1,0.025),
                   plusCI_Tar = quantile(summary$percentage1,0.975),
                   avg_wait = mean(summary$mean_wait),
                   avg_cost = mean(summary$mean_cost),
                   minCI_cost = mean(summary$min_CI_cost),
                   maxCI_cost = mean(summary$max_CI_cost),
                   sd_wait = sd(summary$mean_wait),
                   max_wait = max(summary$mean_wait),
                   minCI_wait = quantile(summary$mean_wait,0.025),
                   plusCI_wait = quantile(summary$mean_wait,0.975),
                   avg_Q = mean(summary$q_length),
                   sd_Q = sd(summary$q_length),
                   min_Q = min(summary$q_length),
                   max_Q = max(summary$q_length),
                   minCI_Q = quantile(summary$q_length,0.025),
                   plusCI_Q = quantile(summary$q_length,0.975),
                   res_used= mean(summary$res_used),
                   mean_idle= mean(summary$res_idle),
                   sd_idle = sd(summary$res_idle),
                   min_idle = min(summary$res_idle),
                   max_idle = max(summary$res_idle),
                   minCI_idle = quantile(summary$res_idle,0.025),
                   plusCI_idle = quantile(summary$res_idle,0.975),
                   in_sys = mean(summary$in_sys))
    head(Final)

    #if the percentage of patients that has to wait longer than one day is more than target percentage then add capacity
   if (Final$percentageTar[u]<PerTar) { 
      add<-add+1
      n_slots <- n_slots+1
    } else {  #else increase LOS
      #avg_LOS<- avg_LOS+1
      #add<-0
      break}
  }

################# data reshape ######
visits <- data.frame(visits)

################# plots #############

png(filename = "9Costs per capacity.png", width = 841, height = 493)
plot_costSD <- ggplot(Final, aes(x = capacity, y = avg_delay_cost))
plot_costSD + geom_point() +
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "grey50", size = 0.5, linetype="solid"))+
  theme(text=element_text(size=20))+
  xlab("P1 capacity (maximum visits possible per day)")+
  ylab("Total acute delay cost + total P1 service cost")+
  ggtitle("Total cost of stable P1 capacity range")
dev.off()


png(filename = "fAverage delayed discharge_error bars.png", width = 841, height = 493)
plot_dtoc <- ggplot(Final, aes(x = capacity, y = avg_Q))
plot_dtoc + geom_point() +
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "grey50", size = 0.5, linetype="solid"))+
  theme(text=element_text(size=20))+
  geom_errorbar(aes(x= capacity, y= avg_Q, ymin = minCI_Q, ymax = plusCI_Q), width =.4, position = position_dodge(.9)) +
  xlab("P1 capacity [maximum visits possible per day]")+
  ylab("Average acute delay (number of patients)")+
  ggtitle("Average delay per stable P1 capacity range")
dev.off()

# png(filename = "fAverage delayed discharge.png", width = 841, height = 493)
# plot_cost <- ggplot(Final, aes(x = capacity, y = avg_Q))
# plot_cost + geom_point() +
#   theme_bw()+
#   theme(panel.border = element_rect(fill = NA, color = "grey50", size = 0.5, linetype="solid"))+
#   theme(text=element_text(size=20))+
#   xlab("P1 capacity [maximum visits possible per day]")+
#   ylab("Average acute delay (number of patients)")+
#   ggtitle("Average delay per stable P1 capacity range")
# dev.off()

png(filename = "fAverage DTOC wait.png", width = 841, height = 493)
plot_delay <- ggplot(Final, aes(x = capacity, y = avg_wait))
plot_delay + geom_point() +
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "grey50", size = 0.5, linetype="solid"))+
  theme(text=element_text(size=20))+
  xlab("P1 capacity [maximum visits possible per day]")+
  ylab("Average acute delay (days)")+
  ggtitle("Average delay per stable P1 capacity range")
dev.off()

#writing results into csv file
write.csv(Final, file = paste("7Base model_Target days_P1_arrivals", wait_target, PerTar, "Percent Simlength_",sim_length, "Seed_",SEED ," Runs_", nruns, " Lambda_", lambda, " LOS_",avg_LOS, " ISR_",ISR,"_",endSR, ".csv"))


