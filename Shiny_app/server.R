
library(shiny)
library(shinydashboard)
library(ggplot2)
library(ggvis)

source('FA_model.R')

server = function(input, output, session) {
  
  values <- reactive({
    
    v=list(gamma=input$gamma,
           delta=input$delta,
           phi=input$phi,
           h=input$h,
           beta=input$beta,
           k=input$k,
           p=input$p,
           q=input$q,
           n=input$n,
           m=input$m,
           lO=input$lO,
           P50=input$P50,
           shape= input$shape,
           lm = input$lm,
           n_int = input$n_int,
           omega=input$omega,
           Ea=input$Ea,
           temp=seq(input$temp[1],input$temp[2],l=input$lm),
           temp_ref=input$Tref,
           Topt=input$Topt,
           O2crit=input$O2crit,
           M=input$M,
           v=input$v)
    
    #browser()
    
    v
    
  })
  
  
  #######################################
  ############  data ####################
  #######################################
  
  plot_data_temp <- reactive({
    
    v = values()
    lm = 10^seq(0,6,l=v$lm)
    #browser()
    out = eval_tau_eq(gamma=v$gamma,
                      delta=v$delta,
                      phi=v$phi,
                      h=v$h,
                      beta=v$beta,
                      k=v$k,
                      p=v$p,
                      q=v$q,
                      n=v$n,
                      m=lm,
                      M=v$M,
                      v=v$v)
    
    pd <- data.frame(m=lm)
    pd$tau_uc <- out
    
    #browser()
    O2 = O2_supply(O2=seq(1,10,l=length(v$temp)),Topt=v$Topt,O2crit=v$O2crit,Tmax=v$temp[length(v$temp)],T=v$Topt,delta=v$lO,omega=v$shape,P50=v$P50,n=v$n)
    tau_o2 = sapply(eval_tau_max_o2(f=O2,omega = v$omega,
                                    gamma=v$gamma,
                                    delta=v$delta,
                                    phi=v$phi,
                                    h=v$h,
                                    beta=v$beta,
                                    k=v$k,
                                    p=v$p,
                                    q=v$q,
                                    n=v$n,
                                    m=v$m
    ),min,out)
    
    tau_o2[tau_o2<0] <- 0
    tau_o2[tau_o2>1] <- 1
    
    o2frame <- data.frame(O2=seq(1,10,l=v$lm),tau_o2 = tau_o2)
    
    O2_tcor <- O2_fact(v$temp,5)
    f=O2_supply(10*O2_tcor,Topt=v$Topt,O2crit=v$O2crit,Tmax=v$temp[length(v$temp)],T=v$temp,delta=v$lO,omega=v$shape,P50=v$P50,n=v$n)
    #browser()
    max_tau <- eval_tau_max_temp(f=f,
                                 temp=v$temp,
                                 Ea=v$Ea,
                                 omega = v$omega,
                                 gamma=v$gamma,
                                 delta=v$delta,
                                 phi=v$phi,
                                 h=v$h,
                                 beta=v$beta,
                                 k=v$k,
                                 p=v$p,
                                 q=v$q,
                                 n=v$n,
                                 m=v$m)
    tau <- eval_tau_eq_temp(temp=v$temp,
                            Ea=v$Ea,
                            gamma=v$gamma,
                            delta=v$delta,
                            phi=v$phi,
                            h=v$h,
                            beta=v$beta,
                            k=v$k,
                            p=v$p,
                            q=v$q,
                            n=v$n,
                            m=v$m,
                            M=v$M,
                            v=v$v)
    
    tau_max = apply(cbind(tau,max_tau),1,min)
    tau_max[tau_max<0] <- 0
    tau_max[tau_max>1] <- 1
    
    model_frame <- model_out(tau_max = tau_max,
                             temp=v$temp,
                             Ea=v$Ea,
                             gamma=v$gamma,
                             delta=v$delta,
                             phi=v$phi,
                             h=v$h,
                             beta=v$beta,
                             k=v$k,
                             p=v$p,
                             q=v$q,
                             n=v$n,
                             m=v$m)
    
    
    #browser()
    scope <- f*v$m^v$n-model_frame$Metabolism*v$omega
    scope[scope<0.0001] <- 0
    #browser()
    
    bind_cols(pd,
              o2frame,
              model_frame,
              data_frame(
                Temperature=v$temp,
                Realised = tau_max,
                Limit = sapply(sapply(max_tau,max,0),min,1),
                M = (tau_max*v$v+v$M)*v$m^(v$q-1),
                Optimum = sapply(sapply(tau,max,0),min,1),
                Scope=ifelse(scope<0.0001,0,scope),
                `MMR` = f*v$m^v$n,
                `Active M.` = model_frame$Metabolism*v$omega,
                `Std M.` = model_frame$Std*v$omega,
                Viable = as.numeric(tau_max>0.0001 & model_frame[['C for growth']]>0.0001))
    )
    
  })
  
  
  plot_data_growth <- reactive({
    
    v = values()
    lm = 10^seq(0,6,l=v$lm)
    tau_uc = eval_tau_eq(gamma=v$gamma,
                         delta=v$delta,
                         phi=v$phi,
                         h=v$h,
                         beta=v$beta,
                         k=v$k,
                         p=v$p,
                         q=v$q,
                         n=v$n,
                         m=lm,
                         M=v$M,
                         v=v$v)
    
    tau_uc[tau_uc<0] <- 0
    tau_uc[tau_uc>1] <- 1
    
    taus <- get_taus(v,tau_uc,10,v$temp_ref)
    
    model_out_par(tau_max = taus$tau_max,
                  tau = taus$tau,
                  tau_uc = tau_uc,
                  tau_o2 = taus$tau_o2,
                  M=v$M,
                  temp=v$temp_ref,
                  Ea=v$Ea,
                  gamma=v$gamma,
                  delta=v$delta,
                  phi=v$phi,
                  h=v$h,
                  beta=v$beta,
                  k=v$k,
                  p=v$p,
                  q=v$q,
                  n=v$n,
                  m = lm
    )
  })
  
  
  plot_data_growth_tO2 <- reactive({
    
    v = values()
    lm = 10^seq(0,6,l=v$lm)
    #browser()
    tau_uc = eval_tau_eq(gamma=v$gamma,
                         delta=v$delta,
                         phi=v$phi,
                         h=v$h,
                         beta=v$beta,
                         k=v$k,
                         p=v$p,
                         q=v$q,
                         n=v$n,
                         m=lm,
                         M=v$M,
                         v=v$v)
    
    tau_uc[tau_uc<0] <- 0
    tau_uc[tau_uc>1] <- 1
    
    #browser()
    n_int <- v$n_int
    O2_range <- seq(1,10,l=n_int)
    winfs <- data.frame(matrix(NA,n_int,5))
    winfs[,1] <- O2_range
    winfs[,2] <- seq(min(v$temp),max(v$temp),l=n_int)
    
    withProgress(message = 'Calculating Winf', value = 0, {
      winfs[,3:5] <- bind_rows(parallel::mclapply(1:n_int, function(i){#parallel::mclapply(1:n_int, function(i){
        taus <- get_taus(v,tau_uc,winfs[i,1],winfs[i,2],m=lm)
        
        #if(all(taus$tau_max<0.05)) return(data_frame(winf=0,winf_o2=0,winf_uc=0))
        mout <- model_out_par(tau_max = taus$tau_max,
                              tau = taus$tau,
                              tau_uc = tau_uc,
                              tau_o2 = taus$tau_o2,
                              M=v$M,
                              v=v$v,
                              temp=winfs[i,2],
                              Ea=v$Ea,
                              gamma=v$gamma,
                              delta=v$delta,
                              phi=v$phi,
                              h=v$h,
                              beta=v$beta,
                              k=v$k,
                              p=v$p,
                              q=v$q,
                              n=v$n,
                              m = lm,
                              ret_winf = T)
        # Increment the progress bar, and update the detail text.
        incProgress(1/n_int, detail = paste("Temp", round(winfs[i,2],2)))
        
        mout
        
      },mc.cores = 4))})
    
    colnames(winfs) <- c('O2','Temperature','Temp','Oxygen','Base')
    winfs
    
  })
  
  #######################################
  ############ plots ####################
  #######################################
  # 
  
  
  output$Tauplot <- renderPlot({
    
    plot_data_temp() %>%
      gather(Activity,value,Optimum, Limit,Realised) %>%
      mutate(Activity = relevel(as.factor(Activity), "Realised")) %>%
      ggplot() + 
      geom_line(aes(x=Temperature, y=value, linetype=Activity),alpha=0.7) +
      theme_cowplot()+
      theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.justification=c(0,1), 
        legend.position=c(0.05,0.85),
        legend.box = "horizontal",
        legend.text  = element_text(size=10),
        legend.title = element_text(size=12),
        legend.spacing = unit(0.005,'npc')) + 
      labs(x = expression("Mean temperature " ( degree*C)),
           y=expression('Activity ' (tau)))
    
  }, bg="transparent",type = "cairo-png")
  
  output$O2plot <- renderPlot({
    plot_data_temp() %>%
      tidyr::gather(Rate,value,MMR, `Active M.`, `Std M.`) %>%
      ggplot() + 
      geom_line(aes(x=Temperature, y=value, linetype=Rate),alpha=0.7) +
      theme_cowplot()+
      theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.justification=c(0,1), 
        legend.position=c(0.05,0.85),
        legend.text  = element_text(size=10),
        legend.title = element_text(size=12),
        legend.spacing = unit(0.005,'npc')#legend.background = element_rect(fill='white')
      )+
      scale_color_discrete(guide='none')+
      scale_shape_discrete(guide='none')+
      labs(x = expression("Mean temperature " ( degree*C)),
           y=expression(Oxygen~(g.yr^-1))) 
    
  }, bg="transparent",type = "cairo-png")
  
  
  output$Eplot <- renderPlot({
    
    plot_data_temp() %>%
      mutate(Consumption=sc(Consumption),
             Efficiency=sc(Efficiency),
             `C for growth` = ifelse(`C for growth`<0,0,`C for growth`),
             `Available Energy` = sc(`C for growth`)) %>% 
      tidyr::gather(Rate,value,`Feeding level`,`Available Energy`,M) %>%
      mutate(Rate = relevel(as.factor(Rate),'Feeding level')) %>%
      ggplot() + 
      geom_line(aes(x=Temperature, y=value, linetype=Rate),alpha=0.7) +
      theme_cowplot()+
      #lims(y=c(0,1))+
      theme(
        panel.background = element_blank(),
        plot.background = element_blank(),legend.justification=c(0,1), 
        legend.position=c(0.05,0.85),
        legend.text  = element_text(size=10),
        legend.title = element_text(size=12),
        legend.spacing = unit(0.005,'npc'))+
      scale_color_discrete(guide='none')+
      scale_shape_discrete(guide='none')+
      labs(x = expression("Mean temperature " ( degree*C)),
           y='Relative rate') 
    
  }, bg="transparent",type = "cairo-png")
  
  output$TGvis <- renderPlot({
    plot_data_growth_tO2() %>%
      mutate(length=lw(Temp)) %>% 
      ggplot() +
      geom_line(aes(x=Temperature, y=length),alpha=0.7) +
      #geom_point(data=growth_scenarios[seq.int(1,nrow(growth_scenarios),by = 5),],aes(x=Temperature, y=Temp, col=Scenario,shape=Scenario),alpha=0.5,size=2) +
      #theme_cowplot()+
      theme(
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.justification=c(0,1), 
        legend.position=c(0.05,0.85),
        legend.text  = element_text(size=10),
        legend.title = element_text(size=12),
        legend.spacing = unit(0.005,'npc'),
        legend.background = element_rect(fill='white'))+
      #scale_color_discrete(guide='none')+
      #scale_shape_discrete(guide='none')+
      #scale_y_log10()+
      labs(x = expression("Mean temperature " ( degree*C)),
           y='Adult length (cm)')
  }, bg="transparent",type = "cairo-png")
}