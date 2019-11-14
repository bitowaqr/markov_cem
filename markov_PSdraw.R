source("./markov_support_PS.R")
library(MCMCpack)


### PSA Set up
  #time1 = Sys.time()
  cycle_length = 84
  start_age = 16 
  DR_QALY = DR_COSTS = 0.035 # discount rate
  PSA_length = 100
  PSA_res = matrix(data=NA,ncol=5,nrow=PSA_length,dimnames = list(1:PSA_length,c("QALY base", "Costs base","QALY int","Costs int","Icer")))
  parameter_trace = matrix(data=NA,ncol=2,nrow=PSA_length)

# draw parameters outside of loop.
  psa.input <- c()
  for(x in 1:200){
  psa.input[[x]] <- draw_params()
  }
  
time1 = Sys.time()
## PSA LOOP
  for(r in 1:PSA_length){
    #cat("\r",r)
    # transistion matrices
    m = psa.input[[r]]
    mat_b = trans_mat(m$tb1,m$tb2,m$tb3,m$tb4)
    mat_int = trans_mat(m$ti1,m$ti2,m$ti3,m$ti4)
    
    # create empty matrices
    markov_trace_base = matrix(data=NA,ncol=dim(mat_b)[1]+1,nrow=1+cycle_length)
    markov_trace_base[1,] = c(m$initial_distribution,0)
    markov_trace_int = markov_trace_base
    aging = 0
    
    ### MARKOV LOOP
      for(i in 1:cycle_length){
        step_matrix_b = death_age_fmat(mat=mat_b,x=start_age+aging)
        step_matrix_int = death_age_fmat(mat=mat_int,x=start_age+aging)
        markov_trace_base[i+1,]= markov_trace_base[i,] %*% step_matrix_b
        markov_trace_int[i+1,]= markov_trace_int[i,] %*% step_matrix_int
        aging = aging + 1
        }
    
    ### EXTRACT RESULTS FROM MARKOV LOOP
    res_base  = get_CE(trace=markov_trace_base,
                       utils = m$state_utils,
                       costs = m$base_state_costs,
                       iv_rates = m$base_iv_days,
                       iv_disutil = m$iv_excer_disutil)
    res_int  = get_CE(trace=markov_trace_int,
                       utils = m$state_utils,
                       costs = m$int_state_costs,
                       extra_costs = m$int_costs_once,
                       iv_rates = m$int_iv_days,
                       iv_disutil = m$iv_excer_disutil)
    Icer = (res_base$state_costs_disc-res_int$state_costs_disc)/ (res_base$QALY_disc - res_int$QALY_disc)
    
    ## save results in PSA trace
    PSA_res[r,] = c(res_base$QALY_disc,
                    res_base$state_costs_disc,
                    res_int$QALY_disc,
                    res_int$state_costs_disc,
                    Icer)
    
    ## parameter trace for EVPPI
    parameter_trace[r,] = c(m$cons,m$gamma)
  }

time2 = Sys.time()
time.diff = round(time2-time1)
time.diff

  
  
####################
### PSA RESULTS
###########
    
  # CE-PLANE PLOT
  incr_Q = PSA_res[,"QALY int"] - PSA_res[,"QALY base"] 
  incr_C = PSA_res[,"Costs int"] - PSA_res[,"Costs base"] 
    library(ggplot2)
    ce_plane_plot = ggplot() +
      geom_point(aes(x=incr_Q,y=incr_C),col="darkblue",alpha=0.8) +
      geom_point(aes(x=mean(incr_Q),y=mean(incr_C)),size=2,col="red") +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      geom_line(aes(x=c(-0.5,0.5),y=c(-10000,10000))) +
      theme_minimal() +
      xlim(c(-max(abs(incr_Q)),max(abs(incr_Q)))) +
      ylim(c(-max(abs(incr_C)),max(abs(incr_C))))  +
      xlab("Incremental QALY") +
      ylab("Incremental Costs") 
    ce_plane_plot
  
  ## CEAC PLOT
    ceac_data = get_ceac()
    ceac_plot = ggplot() +
      geom_line(aes(x=ceac_data$wtp,y=ceac_data$prob_ce)) +
      ylim(c(0,1)) +
      ylab("Prob. cost effective") +
      xlab("Willingness to pay threshold") +
      theme_minimal()+
      ggtitle("Cost acceptibility curve")
    
  ## PSA results summary
    res_summary = data.frame(QALY_base = mean(PSA_res[,"QALY base"] ),
                             QALY_int = mean(PSA_res[,"QALY int"]),
                             incr_Qaly = mean(PSA_res[,"QALY int"] - PSA_res[,"QALY base"] ),
                             cost_base = mean(PSA_res[,"Costs base"]),
                             cost_int = mean(PSA_res[,"Costs int"]),
                             incr_costs = mean(PSA_res[,"Costs int"] - PSA_res[,"Costs base"]),
                             netbenefit = mean(incr_Q*20000 - incr_C),
                             row.names = "Point estimate")
    t(res_summary)
## finish
    