
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
  
  plot_data_temp <- eventReactive(input$go,{
    
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
    
    data.frame(pd,
               o2frame,
               Temperature=v$temp,
               tau_temp = tau_max,
               tau_t = max_tau,
               M = tau^v$v*v$M*v$m^(v$q-1),
               tau_lim = tau,
               model_frame,
               po2=ifelse(scope<0.0001,0,scope),
               so2 = f*v$m^v$n,
               do2 = model_frame$Metabolism*v$omega,
               viable = as.numeric(tau_max>0.0001 & model_frame$e>0.0001))
    
  },ignoreNULL=F)
  
  
  plot_data_growth <- eventReactive(input$go,{

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
  },ignoreNULL=F)

  
  plot_data_growth_tO2 <- eventReactive(input$go,{
    
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
    
  },ignoreNULL=F)
  
  #######################################
  ############ plots ####################
  #######################################
  # 
  # reactive({  plot_data_temp %>%
  #     ggvis(x=~m, y=~tau_uc) %>%
  #     add_axis("x", title = 'Mass (g)', grid=F) %>%
  #     scale_numeric("x",trans='log',expand=0,nice=T) %>%
  #     add_axis("y", title = "\u03C4") %>%
  #     scale_numeric("y", domain = c(0, 1), clamp=T) %>%
  #     layer_lines()%>% 
  #     layer_points(size:=10)%>%
  #     set_options(height = 400, width = 460) })  %>%
  #   bind_shiny(plot_id = "ggvisplot")
  # 
  plot_data_temp %>%
    ggvis(x=~Temperature, y=~M) %>%
    add_axis("y", title = "Realised mortality m") %>%
    scale_numeric("y", domain = c(0, 1), clamp=T) %>%
    layer_lines()%>% 
    layer_points(size:=10)%>%
    set_options(height = 400, width = 460)%>%
    bind_shiny(plot_id = "ggvisO2plot")
  
  plot_data_temp %>%
    ggvis(x=~Temperature, y=~tau_temp) %>%
    add_axis("y", title = "\u03C4") %>%
    scale_numeric("y", domain = c(0, 1), clamp=T) %>%
    layer_points(fill=~viable,stroke=~viable,size:=10)%>%
    layer_lines(x=~Temperature, y=~tau_t, strokeDash:=6)%>% 
    layer_lines(x=~Temperature, y=~tau_lim, strokeDash:=2)%>% 
    hide_legend(c("fill","stroke"))%>%
    set_options(height = 400, width = 460)%>%
    bind_shiny(plot_id = "ggvisTempplot")
  
  plot_data_temp %>%
    ggvis(x=~Temperature, y=~po2) %>%
    add_axis("y", title = 'P\u1D0F\u2082 (g/yr)') %>%
    layer_lines()%>%
    layer_points(fill=~viable,stroke=~viable,size:=10,size:=10)%>%
    layer_lines(x=~Temperature, y=~do2, strokeDash:=6)%>% 
    layer_lines(x=~Temperature, y=~so2, strokeDash:=2)%>% 
    hide_legend(c("fill","stroke"))%>%
    set_options(height = 400, width = 460)%>%
    bind_shiny(plot_id = "ggvisMetplot")
  
  plot_data_temp %>%
    ggvis(x=~Temperature, y=~feeding) %>%
    add_axis("y", title = 'Feeding level (f)') %>%
    layer_lines()%>%
    layer_points(fill=~viable,stroke=~viable,size:=10)%>%
    hide_legend(c("fill","stroke"))%>%
    set_options(height = 400, width = 460)%>%
    bind_shiny(plot_id = "ggvisfeedplot")
  
  plot_data_temp %>%
    ggvis(x=~Temperature, y=~e) %>%
    add_axis("y", title = 'Available energy (P_C)') %>%
    layer_lines()%>%
    scale_numeric("y", domain = c(0, NA), clamp=T) %>%
    layer_points(fill=~viable,stroke=~viable,size:=10)%>%
    hide_legend(c("fill","stroke"))%>%
    set_options(height = 400, width = 460)%>%
    bind_shiny(plot_id = "ggvisEplot")
  
  plot_data_temp %>%
    ggvis(x=~Temperature, y=~Efficiency) %>%
    add_axis("y", title = 'Efficiency (P_C/f)') %>%
    layer_lines()%>%
    layer_points(fill=~viable,stroke=~viable,size:=10)%>%
    hide_legend(c("fill","stroke"))%>%
    set_options(height = 400, width = 460)%>%
    bind_shiny(plot_id = "ggvisEffplot")
  
  # plot_data_temp %>%
  #   ggvis(x=~Temperature, y=~predation) %>%
  #   add_axis("y", title = 'Predation rate (f/Θ)') %>%
  #   layer_lines()%>%
  #   layer_points(fill=~viable,stroke=~viable,size:=10)%>%
  #   hide_legend(c("fill","stroke"))%>%
  #   set_options(height = 400, width = 460)%>%
  #   bind_shiny(plot_id = "ggvispredplot")
  
  # plot_data_growth %>%
  #   ggvis(x=~m, y=~Rate,fill=~Term,size=3) %>%
  #   scale_numeric("x", trans="log",expand=0)%>%
  #   scale_numeric("y", trans="log",expand=0)%>%
  #   add_axis("y", grid=F) %>%
  #   add_axis("x", grid=F) %>%
  #   group_by(Limitation,Term) %>%
  #   layer_points(shape=~Limitation)%>%
  #   add_legend(scale="fill", properties = legend_props(legend = list(y = 200))) %>%
  #   #layer_paths()%>%
  #  
  #   set_options(height = 300, width = 720)%>%
  #   #scale_nominal("stroke",label='Limitation') %>%
  #   #scale_nominal("strokeDash",label='Term') %>%
  #   bind_shiny(plot_id = "ggvisGvis")%>%
  #   set_options(duration = 0)
  # 
  # plot_data_growth_tO2 %>%
  #   ggvis(x=~O2, y=~Oxygen,stroke:='orange') %>%
  #   layer_lines()%>%
  #   layer_lines(x=~O2, y=~Base,stroke:='black')%>%
  #   scale_numeric("y",label = 'm\u221E (g)',expand=0)%>%
  #   set_options(height = 300, width = 300)%>%
  #   bind_shiny(plot_id = "ggvisOGvis")
  
  plot_data_growth_tO2 %>%
    ggvis(x=~Temperature, y=~Temp,stroke:='green') %>%
    scale_numeric("y", trans="log",expand=1)%>%
    layer_lines()%>%
    #layer_lines(x=~Temperature, y=~Base,stroke:='black')%>%
    scale_numeric("y",label = 'm\u221E (g)',expand=0)%>%
    set_options(height = 500, width = 720)%>%
    bind_shiny(plot_id = "ggvisTGvis")
  
}