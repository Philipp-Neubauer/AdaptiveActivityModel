
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
           m=max(10^input$m)/4,
           max_m=max(input$m),
           min_m=min(input$m),
           slope=input$slope,
           tmax=input$tmax,
           tr=input$tr,
           c=input$c,
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
           v=input$v,
           nu=input$nu, 
           dt = input$dt)
    
    #browser()
    
    v
    
  })
  
  
  #######################################
  ############  data ####################
  #######################################
  
  plot_data_temp <- reactive({
    
    v = values()
    lm = 10^seq(v$min_m,v$max_m,l=v$lm)
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
                                 temp_ref=v$temp_ref,
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
                            temp_ref=v$temp_ref,
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
                             temp_ref=v$temp_ref,
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
                M = (tau_max*v$v+v$M)*v$m^v$nu,
                Optimum = sapply(sapply(tau,max,0),min,1),
                Scope=ifelse(scope<0.0001,0,scope),
                `MMR` = f*v$m^v$n,
                `Active M.` = model_frame$Metabolism*v$omega,
                `Std M.` = model_frame$Std*v$omega,
                Viable = as.numeric(tau_max>0.0001 & model_frame[['C for growth']]>0.0001))
    )
    
  })
  
  plot_data_growth_tO2 <- eventReactive(input$go,{
    
    v = values()
    lm = seq(10^v$min_m,10^v$max_m,l=v$lm)
    #browser()
    
    #browser()
    n_int <- v$n_int
    O2_range <- seq(1,10,l=n_int)
    
    #### Step 1 get mstar=mo+slope*age for all Temp
    #browser()
    
    mout <- model_out_growth(temp=seq(min(v$temp),max(v$temp),l=n_int),
                             temp_ref=v$temp_ref,
                             l=v$lm,
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
                             tmax = v$tmax,
                             slope=v$slope,
                             tr=v$tr,
                             v=v,
                             lm=lw(lm),
                             dt=v$dt)
    
    winfs <- mout$winfs
    
    ### Step 2 calc mat for selected response 
    
    #browser()
    
    mouts <- model_out_growth_check(temp=winfs$Temperature,
                                    temp_ref=v$temp_ref,
                                    Ea=v$Ea,
                                    gamma= seq(v$gamma*(1-v$c),v$gamma*(1+v$c),l=n_int),
                                    delta=v$delta,
                                    phi=v$phi,
                                    h=v$h,
                                    beta=v$beta,
                                    k=v$k,
                                    p=v$p,
                                    q=v$q,
                                    n=v$n,
                                    mstar = lw(lm[winfs$opt[which.min(abs(winfs$Temperature-v$temp_ref))]]),
                                    tmax = v$tmax,
                                    slope=v$slope,
                                    tr=v$tr,
                                    v=v,
                                    v$dt,
                                    lm=lw(lm))
    
    list(winfs=winfs,
         growth=mout$growth,
         g_growth=mouts$g_growth,
         t_growth=mouts$t_growth)
    
  })
  
  # 
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
    plot_data_growth_tO2()[['winfs']] %>%
      mutate(length=lw(l)) %>% 
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
  
  output$Gvis <- renderPlot({
    plot_data_growth_tO2()[['winfs']] %>%
      ggplot() +
      geom_line(aes(x=Temperature, y=L),alpha=0.7) +
      #geom_point(data=growth_scenarios[seq.int(1,nrow(growth_scenarios),by = 5),],aes(x=Temperature, y=Temp, col=Scenario,shape=Scenario),alpha=0.5,size=2) +
      #theme_cowplot()+
      theme(
        panel.background = element_blank(),
        plot.background = element_blank())+
      #scale_color_discrete(guide='none')+
      #scale_shape_discrete(guide='none')+
      #scale_y_log10()+
      labs(x = expression("Mean temperature " ( degree*C)),
           y=expression(G~(mm.y^-1))) 
  }, bg="transparent",type = "cairo-png")
  
  output$dPm <- renderPlot({
    
    dPm <- plot_data_growth_tO2()[['dPM']] %>%
      mutate(length=lw(Mass))
    
    ggplot(dPm) +
      geom_line(aes(x=length, y=dPm,col=as.factor(Temperature))) +
      geom_hline(aes(yintercept=1))+
      #geom_point(data=growth_scenarios[seq.int(1,nrow(growth_scenarios),by = 5),],aes(x=Temperature, y=Temp, col=Scenario,shape=Scenario),alpha=0.5,size=2) +
      #theme_cowplot()+
      viridis::scale_colour_viridis(discrete = T,guide='none')+
      theme(
        panel.background = element_blank(),
        plot.background = element_blank())+
      #scale_color_discrete(guide='none')+
      #scale_shape_discrete(guide='none')+
      coord_cartesian(ylim=c(0.5,max(dPm$dPm)),xlim=c(0,max(dPm$length[dPm$dPm>0.4])))+
      labs(y = 'd(P/m)/dm',
           x='Length (cm)')
  }, bg="transparent",type = "cairo-png")
  
  output$Pm <- renderPlot({
    
    Pm <- plot_data_growth_tO2()[['dPM']] %>%
      mutate(length=lw(Mass))
    
    ggplot(Pm) +
      geom_line(aes(x=length, y=Pm,col=as.factor(Temperature))) +
      geom_hline(aes(yintercept=1))+
      #geom_point(data=growth_scenarios[seq.int(1,nrow(growth_scenarios),by = 5),],aes(x=Temperature, y=Temp, col=Scenario,shape=Scenario),alpha=0.5,size=2) +
      #theme_cowplot()+
      viridis::scale_colour_viridis(discrete = T,guide='none')+
      theme(
        panel.background = element_blank(),
        plot.background = element_blank())+
      #scale_color_discrete(guide='none')+
      #scale_shape_discrete(guide='none')+
      coord_cartesian(ylim=c(0,max(Pm$Pm)),xlim=c(0,max(Pm$length[Pm$Pm>-1])))+
      labs(y = 'P/m',
           x='Length (cm)')
  }, bg="transparent",type = "cairo-png")
  
  # 
  output$am <- renderPlot({
    
    alloc <- plot_data_growth_tO2()$g_growth
    #browser()
    norm <- alloc %>% 
      group_by(Gamma) %>% 
      summarise(ts=t[ifelse(any(abs(allocs-0.5)<0.01),which.min(abs(allocs-0.5)),NA)-1],
                m=g_length[ifelse(any(abs(allocs-0.5)<0.01),which.min(abs(allocs-0.5)),NA)-1]) %>%
      filter(!is.na(ts))
    
    maxx <- max(norm$ts,na.rm = T)+max(norm$ts,na.rm = T)/2
    #
    ggplot() +
      geom_line(aes(x=t, y=g_length,col=as.factor(Gamma)),data=alloc,linetype=2) +
      geom_point(aes(x=ts, y=m),size=1.2,data=norm) +
      #geom_point(data=growth_scenarios[seq.int(1,nrow(growth_scenarios),by = 5),],aes(x=Temperature, y=Temp, col=Scenario,shape=Scenario),alpha=0.5,size=2) +
      #theme_cowplot()+
      viridis::scale_colour_viridis(discrete = T,guide='none')+
      theme(
        panel.background = element_blank(),
        plot.background = element_blank())+
      #scale_color_discrete(guide='none')+
      #scale_shape_discrete(guide='none')+
      coord_cartesian(xlim=c(0,maxx),
                      ylim = c(0,norm$m[norm$ts<maxx]))+
      labs(y = 'Length (cm)',
           x='Age (years)')
  }, bg="transparent",type = "cairo-png")
  
  
  output$alloc <- renderPlot({
    #browser()
    plot_data_growth_tO2()[['t_growth']] %>%
      ggplot() +
      geom_line(aes(x=t, y=allocs,col=as.factor(Temperature))) +
      geom_hline(aes(yintercept=1))+
      #geom_point(data=growth_scenarios[seq.int(1,nrow(growth_scenarios),by = 5),],aes(x=Temperature, y=Temp, col=Scenario,shape=Scenario),alpha=0.5,size=2) +
      #theme_cowplot()+
      viridis::scale_colour_viridis(discrete = T,guide='none')+
      theme(
        panel.background = element_blank(),
        plot.background = element_blank())+
      #scale_color_discrete(guide='none')+
      #scale_shape_discrete(guide='none')+
      coord_cartesian(ylim = c(0,1))+
      labs(y = 'P for repro',
           x='Age  (years)')
  }, bg="transparent",type = "cairo-png")
  
  # output$efg <- renderPlot({
  #   
  #   alloc <- plot_data_growth_tO2()[['alloc']] %>%
  #     mutate(length=lw(Mass))
  #   
  #   ggplot(alloc) +
  #     geom_line(aes(x=length, y=efg,col=as.factor(Temperature))) +
  #     geom_hline(aes(yintercept=1))+
  #     #geom_point(data=growth_scenarios[seq.int(1,nrow(growth_scenarios),by = 5),],aes(x=Temperature, y=Temp, col=Scenario,shape=Scenario),alpha=0.5,size=2) +
  #     #theme_cowplot()+
  #     viridis::scale_colour_viridis(discrete = T,guide='none')+
  #     theme(
  #       panel.background = element_blank(),
  #       plot.background = element_blank())+
  #     #scale_color_discrete(guide='none')+
  #     #scale_shape_discrete(guide='none')+
  #     coord_cartesian()+
  #     labs(y = 'E for growth',
  #          x='Length (cm)')
  # }, bg="transparent",type = "cairo-png")
  # 
  output$ls <- renderPlot({
    
    alloc <- plot_data_growth_tO2()[['t_growth']]
    #browser()
    norm <- alloc %>% group_by(Temperature) %>% summarise(ts=t[ifelse(any(abs(allocs-0.5)<0.1),which.min(abs(allocs-0.5)),NA)-1],
                                                          t1=t[ifelse(any(abs(allocs-0.1)<0.1),which.min(abs(allocs-0.1)),NA)-1],
                                                          t3=t[ifelse(any(abs(allocs-0.9)<0.1),which.min(abs(allocs-0.9)),NA)-1],
                                                          m=t_length[ifelse(any(abs(allocs-0.5)<0.1),which.min(abs(allocs-0.5)),NA)-1],
                                                          
                                                          q1=t_length[ifelse(any(abs(allocs-0.1)<0.1),which.min(abs(allocs-0.1)),NA)-1],
                                                          q3=t_length[ifelse(any(abs(allocs-0.9)<0.1),which.min(abs(allocs-0.9)),NA)-1]) 
    
    ggplot() +
      geom_line(aes(x=t, y=t_length,col=as.factor(Temperature)),data=alloc) +
      geom_point(aes(x=ts, y=m), col='orange2',data=norm) +
      geom_point(aes(x=t1, y=q1),col='skyblue4',data=norm) +
      geom_point(aes(x=t3, y=q3),col='skyblue4',data=norm) +
      #geom_point(data=growth_scenarios[seq.int(1,nrow(growth_scenarios),by = 5),],aes(x=Temperature, y=Temp, col=Scenario,shape=Scenario),alpha=0.5,size=2) +
      #theme_cowplot()+
      viridis::scale_colour_viridis(discrete = T,guide='none')+
      theme(
        panel.background = element_blank(),
        plot.background = element_blank())+
      coord_cartesian()+
      labs(y = 'Length (cm)',
           x='Age  (years)')
  }, bg="transparent",type = "cairo-png")
  
  
  output$norm <- renderPlot({
    
    alloc <- plot_data_growth_tO2()[['t_growth']]
    #browser()
    norm <- alloc %>% group_by(Temperature) %>% summarise(ts=t[ifelse(any(abs(allocs-0.5)<0.05),which.min(abs(allocs-0.5)),NA)-1],
                                                          ls=t_length[ifelse(any(abs(allocs-0.5)<0.05),which.min(abs(allocs-0.5)),NA)-1],
                                                          l=max(t_length, na.rm=T),
                                                          alloc=allocs[ifelse(any(abs(allocs-0.5)<0.05),which.min(abs(allocs-0.5)),NA)-1]) 
    
    ggplot() +
      geom_line(aes(x=Temperature, y=ls),col='orange2',data=norm) +
      geom_line(aes(x=Temperature, y=l),col='skyblue4',data=norm) +
      #geom_point(data=growth_scenarios[seq.int(1,nrow(growth_scenarios),by = 5),],aes(x=Temperature, y=Temp, col=Scenario,shape=Scenario),alpha=0.5,size=2) +
      #theme_cowplot()+
      #viridis::scale_colour_viridis(discrete = T,guide='none')+
      theme(
        panel.background = element_blank(),
        plot.background = element_blank())+
      #scale_color_discrete(guide='none')+
      #scale_shape_discrete(guide='none')+
      coord_cartesian()+
      labs(y = 'Length (cm)',
           x='Temperature')
  }, bg="transparent",type = "cairo-png")
  
  output$la <- renderPlot({
    
    alloc <- plot_data_growth_tO2()[['tg_alloc']] 
    #browser()
    norm <- alloc %>% group_by(Temperature) %>% summarise(ts=t[ifelse(any(abs(allocs-0.5)<0.1),which.min(abs(allocs-0.5)),NA)-1],
                                                          tq=t[ifelse(any(abs(allocs-0.75)<0.1),which.min(abs(allocs-0.75)),NA)-1],
                                                          tp=t[ifelse(any(abs(allocs-0.25)<0.1),which.min(abs(allocs-0.25)),NA)-1],
                                                          m=ls[ifelse(any(abs(allocs-0.5)<0.1),which.min(abs(allocs-0.5)),NA)-1],
                                                          qs=ls[ifelse(any(abs(allocs-0.75)<0.1),which.min(abs(allocs-0.75)),NA)-1],
                                                          ps=ls[ifelse(any(abs(allocs-0.25)<0.1),which.min(abs(allocs-0.25)),NA)-1]) 
    m <- norm %>% select(-ts,-tq,-tp) %>% tidyr::gather(Length,mass,m,qs,ps)
    t <- norm %>% select(-m,-qs,-ps) %>% tidyr::gather(Age,age,ts,tq,tp)
    norm <- cbind(m,t%>% select(-Temperature))
    
    ggplot() +
      geom_line(aes(x=age, y=mass,col=as.factor(Temperature)),data=norm) +
      geom_point(aes(x=age, y=mass,col=as.factor(Temperature)),data=norm) +
      #geom_point(data=growth_scenarios[seq.int(1,nrow(growth_scenarios),by = 5),],aes(x=Temperature, y=Temp, col=Scenario,shape=Scenario),alpha=0.5,size=2) +
      #theme_cowplot()+
      viridis::scale_colour_viridis(discrete = T,guide='none')+
      theme(
        panel.background = element_blank(),
        plot.background = element_blank())+
      #scale_color_discrete(guide='none')+
      #scale_shape_discrete(guide='none')+
      coord_cartesian(xlim = c(min(norm$age, na.rm = T),max(norm$age, na.rm = T)/2))+
      labs(y = 'Length (cm)',
           x='Age (years)')
  }, bg="transparent",type = "cairo-png")
  
  output$MEplot <- renderPlot({
    #browser()
    plot_data_temp() %>%
      
      select(Cond,Cond_mort) %>%
      ggplot() + 
      geom_line(aes(y=Cond_mort, x=Cond)) +
      theme_cowplot()+
      theme(
        panel.background = element_blank(),
        plot.background = element_blank())+
      labs(x = 'Available Energy',
           y='M') 
    
  }, bg="transparent",type = "cairo-png")
  
}