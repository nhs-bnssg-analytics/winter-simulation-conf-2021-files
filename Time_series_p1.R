library(doSNOW)
library(foreach)
library(tidyverse)

rm(list=ls())
setwd("/home/alison/Dropbox/03_IPAC/R_code/VisitsBasedModel")

#p1 L0S calc from Sirona Access and Flow weekly report 19/02/21 per LA
Bristol_los <- mean(c(11.3, 7.6, 8.7, 7.4, 10.9, 13.4, 10.7, 13.6, 15.9, 14.9))
NorthSomerset_los <- mean(c(19.4, 18.8, 14.9, 18, 16.6, 23, 13, 13, 15.4, 18.4))
SouthGlos_los <- mean(c(9.4, 9.3, 10.8, 6.5, 10.8, 10.7, 9.7, 9.7, 11.8))
pop_Bristol<-463400
pop_NorthSomerset<-213919
pop_SouthGlos<-282644
total_pop <- pop_Bristol + pop_NorthSomerset + pop_SouthGlos
avg_LOS <- (pop_Bristol*Bristol_los+pop_NorthSomerset*NorthSomerset_los+pop_SouthGlos*SouthGlos_los)/total_pop

SEED <-1
nruns <- 200 
sim_length <- 40 #days
warmup <-1000 

#####input variables#####

lambda <- 12.6 #P1 arrival rate non-covid uplift

#set length of stay, initial visits (ISR), and final visits (endSR)
ISR <- 4
endSR <- 2
avg_LOS <- avg_LOS

n_patients = 170 #number of patients system can take
n_slots  <- n_patients * mean(c(ISR, endSR)) #average visits

date_ts = data.frame(seq(as.Date(format(Sys.time(), "%Y-%m-%d")), by = "day", length.out = sim_length))

    ####runs##################################
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
        x<- round(rnorm(1,mean=ISR,sd=0.5))
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
        x<- round(rnorm(1,mean=endSR,sd=0.5))
        if (x<=0){ #if x=0 or negative, then x = 1, there needs to be at least one visit per day to start with
          x <-1 
          return(x)
        } else if (x>ISR){
          x<-ISR
          return(x)
         # } else if (x > n_slots){ #ISR can not be greater then the total number of slots available in one day
         #   x <-n_slots  
         #   return(x)
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
                         in_sys=numeric(sim_length) #number of patients in the system
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
      
      ##### SIMULATION ####
      
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
      
      list<-list(output, resources, waittime_vec)
      
      return(list)
      
    }
    
    stopCluster(cl)
    ###############################################
    
    #creating dataframe for summary info
    summary <- data.frame(LOS = integer(nruns),
                          ISR = integer(nruns),
                          nruns = integer(nruns),
                          sim_length = integer(nruns),
                          warm_up=integer(nruns),
                          capacity = integer(nruns),
                          mean_wait= numeric(nruns),
                          q_length = numeric(nruns),
                          res_used= numeric(nruns),
                          res_idle= numeric(nruns),
                          in_sys = numeric(nruns))
    
    #splitting up RESULTS list in 2
    output<-RESULTS[,1]
    out<-do.call(rbind, output)
    #combining in one dataframe
    
    resources<-RESULTS[,2]
    res<-do.call(cbind, resources) 
    colnames(res)<- c(1:nruns)
    
    waittimes <- RESULTS[,3]
    wait<-do.call(rbind, waittimes)

    #summary of all runs
    for (k in 1:nruns){ 
      r.out <- which(out[,1]==k)
      k.wait <- which(wait[,1]==k)
      summary[k,]<- c(LOS = avg_LOS,
                      ISR = ISR,
                      nruns = nruns,
                      sim_length = sim_length,
                      warm_up=warmup,
                      capacity = n_slots,
                      mean_wait= round(mean(wait$waittime[k.wait]),2),
                      q_length = round(mean(out$q_length[r.out]),2),
                      res_used= round(mean(out$res_used[r.out]),2),
                      res_idle= round(mean(out$res_idle[r.out]),2),
                      in_sys= round(mean(out$in_sys[r.out]),2) ) 
    }
########### to csv and manipulating output ##############
    
write.csv(out, file = "output_by_run.csv")
ts_output = read.csv("output_by_run.csv")
head(ts_output)

#take mean and quantiles per day across all runs
ts_output_quants <- ts_output %>% group_by(day) %>% 
  summarise(q05_q_length=quantile(q_length, 0.05), q5_q_length=quantile(q_length, 0.5), 
            q95_q_length=quantile(q_length, 0.95), mean_q_length=mean(q_length), sd_q_length=sd(q_length),
            q10_q_length=quantile(q_length, 0.1), q25_q_length=quantile(q_length, 0.25), 
            q75_q_length=quantile(q_length, 0.75), q90_q_length=quantile(q_length, 0.9),
            q05_in_sys=quantile(in_sys, 0.05), q5_in_sys=quantile(in_sys, 0.5), 
            q95_in_sys=quantile(in_sys, 0.95), mean_in_sys=mean(in_sys), sd_in_sys=sd(in_sys),
            q10_in_sys=quantile(in_sys, 0.1), q25_in_sys=quantile(in_sys, 0.25), 
            q75_in_sys=quantile(in_sys, 0.75), q90_in_sys=quantile(in_sys, 0.9), 
            mean_res_idle=mean(res_idle), mean_res_used=mean(res_used), 
            sd_res_used = sd(res_used),sd_res_idle = sd(res_idle),
            q05_res_used=quantile(res_used, 0.05), q5_res_used=quantile(res_used, 0.5), 
            q95_res_used=quantile(res_used, 0.95), q10_res_used=quantile(res_used, 0.1), 
            q25_res_used=quantile(res_used, 0.25), q90_res_used = quantile(res_used, 0.9),
            q75_res_used = quantile(res_used, 0.75))
head(ts_output_quants)
ts_output_quants <- cbind(date_ts, ts_output_quants) %>% rename_at(1, ~"Date")
ts_output <-cbind(date_ts, ts_output)%>% rename_at(1, ~"Date")

#####################   PLOTS   ########################

png(filename = paste("f_Optimal: resources utilisation_quants.png"), width = 841, height = 493)
ggplot(data=ts_output_quants, aes(x=Date, y=mean_res_used, group = 1)) +
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "grey50", size = 0.5, linetype="solid"))+
  geom_ribbon(aes(ymin = (q05_res_used), ymax = (q95_res_used)), fill = "cadetblue3", alpha = 0.3) +
  geom_ribbon(aes(ymin = (q10_res_used), ymax = (q90_res_used)), fill = "cadetblue3", alpha = 0.3) +
  geom_ribbon(aes(ymin = (q25_res_used), ymax = (q75_res_used)), fill = "cadetblue3", alpha = 0.3)+
  geom_line(color="black", linetype = "dashed")+
  theme(text=element_text(size=20))+
  ylab("Utilisation of visits")+
  xlab("Date")+
  theme(text=element_text(size=20))+
  ggtitle("Optimal capacity scenario: Utilisation of visits (median with quantiles)")
dev.off() 

png(filename = paste("f_Optimal_capacity_q_length_quant.png"), width = 841, height = 493)
ggplot(data=ts_output_quants, aes(x=Date, y=mean_q_length, group = 1)) +
  theme_bw()+
  #geom_ribbon(aes(ymin = 0, ymax =(mean_q_length+sd_q_length)), fill ="violet", alpha =0.3)+
  theme(panel.border = element_rect(fill = NA, color = "grey50", size = 0.5, linetype="solid"))+
  geom_ribbon(aes(ymin = (q05_q_length), ymax = (q95_q_length)), fill = "darkslategray4", alpha = 0.3) +
  geom_ribbon(aes(ymin = (q10_q_length), ymax = (q90_q_length)), fill = "darkslategray4", alpha = 0.3) +
  geom_ribbon(aes(ymin = (q25_q_length), ymax = (q75_q_length)), fill = "darkslategray4", alpha = 0.3)+
  #geom_line(aes(x=Date, y=q5_q_length), color = "black")+
  geom_line(color="black", linetype = "dashed")+
  theme(text=element_text(size=20))+
  ylab("Number of patients delayed")+
  theme(text=element_text(size=20))+
  xlab("Date")+
  ggtitle("Optimal capacity scenario: Number of acute patients delayed per day")
dev.off() 

png(filename = paste("f_Optimal_capacity_number_in_system_quant.png"), width = 841, height = 493)
ggplot(data=ts_output_quants, aes(x=Date, y=mean_in_sys, group = 1)) +
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "grey50", size = 0.5, linetype="solid"))+
  #geom_ribbon(aes(ymin = (mean_in_sys - sd_in_sys), ymax =(mean_in_sys+sd_in_sys)), fill ="violet", alpha =0.3)+
  geom_ribbon(aes(ymin = (q05_in_sys), ymax = (q95_in_sys)), fill = "darkslategray4", alpha = 0.3) +
  geom_ribbon(aes(ymin = (q10_in_sys), ymax = (q90_in_sys)), fill = "darkslategray4", alpha = 0.3) +
  geom_ribbon(aes(ymin = (q25_in_sys), ymax = (q75_in_sys)), fill = "darkslategray4", alpha = 0.3)+
  geom_line(color="black", linetype = "dashed")+
  theme(text=element_text(size=20))+
  theme(text=element_text(size=20))+
  #geom_line(aes(x=Date, y=q5_in_sys), color = "black")+
  ylab("Number of patients in P1 service")+
  xlab("Date")+
  ggtitle("Optimal capacity scenario: Number of P1 patients in service per day")
dev.off()