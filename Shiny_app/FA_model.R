require(dplyr)
require(purrr)
require(tidyr)
require(ggplot2)
require(cowplot)

lw <- function(w) (w/0.01)^(1/3)
wl <- function(l) 0.01*l^3
sc <- function(x) x/max(x)
logit <- function(p) log(p/(1-p))
inv_logit3 <- function(m,mstar,a,b,c) {
  z = (m-mstar)*cos(atan(b))-a*sin(atan(b))
  1/(1+exp(-(c*z)))
}

# Activity model

eval_tau_eq_temp <- function(Ea,
                             temp,
                             temp_ref=15,
                             gamma=50,
                             delta=2,
                             phi=10,
                             h=30,
                             beta=0.75,
                             k=2,
                             p=0.8,
                             q=0.9,
                             n=0.8,
                             m=100,
                             M=0.2,
                             v=1){
  
  tc = exp(Ea*((temp+(288.2-temp_ref))-288.2)/(8.6173324*10^(-5)*(temp+(288.2-temp_ref))*288.2))
  
  -(M*delta*h*k*m^(n - p + q)*tc - h*k*m^(n - p + q)*tc*v - sqrt(-((beta - 1)*h*k*m^(n - 2*p + 3*q) + h*k*m^(n - 2*p + 3*q)*phi)*tc^2*v^2 - ((beta - 1)*delta*gamma*k*m^(n - p + 2*q) + delta*gamma*k*m^(n - p + 2*q)*phi)*M^2*tc + (((beta - 1)*delta*h*k*m^(n - 2*p + 3*q) + delta*h*k*m^(n - 2*p + 3*q)*phi)*tc^2 + (gamma*h*m^(-p + 3*q)*phi^2 + (beta - 1)*gamma*k*m^(n - p + 2*q) + (beta^2 - 2*beta + 1)*gamma*h*m^(-p + 3*q) + (2*(beta - 1)*gamma*h*m^(-p + 3*q) + gamma*k*m^(n - p + 2*q))*phi)*tc)*M*v)*h)/(M*delta*gamma*k*m^n - ((beta - 1)*gamma*h*m^q + gamma*h*m^q*phi + gamma*k*m^n)*v)
  
# -(h*k*m^(n + p + q)*tc*v + (delta*h*k*m + h*k*m^n)*M*m^(p + q)*tc - sqrt(-((beta - 1)*h*k*m^(n + 2*p + 3*q) + h*k*m^(n + 2*p + 3*q)*phi)*tc^2*v^2 + ((beta - 1)*delta*gamma*k*m^(3*p + 2*q + 1) + (beta - 1)*gamma*k*m^(n + 3*p + 2*q) + (delta*gamma*k*m^(3*p + 2*q + 1) + gamma*k*m^(n + 3*p + 2*q))*phi)*M^2*tc - (((beta - 1)*delta*h*k*m^(2*p + 3*q + 1) + (beta - 1)*h*k*m^(n + 2*p + 3*q) + (delta*h*k*m^(2*p + 3*q + 1) + h*k*m^(n + 2*p + 3*q))*phi)*tc^2 - (gamma*h*m^(3*p + 3*q)*phi^2 + (beta - 1)*gamma*k*m^(n + 3*p + 2*q) + (beta^2 - 2*beta + 1)*gamma*h*m^(3*p + 3*q) + (2*(beta - 1)*gamma*h*m^(3*p + 3*q) + gamma*k*m^(n + 3*p + 2*q))*phi)*tc)*M*v)*h)/((delta*gamma*k*m + gamma*k*m^n)*M*m^(2*p) + ((beta - 1)*gamma*h*m^(2*p + q) + gamma*h*m^(2*p + q)*phi + gamma*k*m^(n + 2*p))*v)

  
  
}

model_out <- function(tau_max,
                      temp,
                      temp_ref=15,
                      Ea,
                      r=0.2,
                      gamma=50,
                      delta=2,
                      phi=0.15,
                      h=30,
                      beta=0.2,
                      k=2,
                      p=0.8,
                      q=0.9,
                      n=0.8,
                      m=100){
  
  tc = exp(Ea*((temp+(288.2-temp_ref))-288.2)/(8.6173324*10^(-5)*(temp+(288.2-temp_ref))*288.2))
  
  f <- tau_max*gamma*m^p/(tau_max*gamma*m^p+tc*h*m^q)
  inp <- (1-phi-beta)*f*tc*h*m^q
  out <- (1+tau_max*delta)*k*tc*m^n
  e <-  inp -out
  efficiency <- e/f
  efficiency[efficiency<0] <- 0
  predation_rate <- f/phi
  
  met = beta*f*tc*h*m^q+ ((1+tau_max*delta)*k*tc*m^n)
  
  #browser()
  
  data_frame(`Feeding level`= f, 
             Consumption = inp,
             `C used for Metabolism` = out,
             `C for growth` = e, 
             Efficiency = efficiency, 
             `Predation rate`=predation_rate, 
             Metabolism = met,
             Std = k*tc*m^n)
}
get_dPdm <- function(m=NULL,tau,
                     tc,
                     nu=-0.25,
                     v=0.2,
                     M=0.2,
                     gamma=50,
                     delta=2,
                     phi=0.15,
                     h=30,
                     beta=0.2,
                     k=2,
                     p=0.8,
                     q=0.9,
                     n=0.8) ((beta + phi - 1)*gamma*h*m^p*m^q*tau*tc/(gamma*m^p*tau + h*m^q*tc) + (delta*tau + 1)*k*m^n*tc)*m^(nu - 1)*nu/((tau*v + M)*m^(2*nu)) - ((beta + phi - 1)*gamma*h*m^(p - 1)*m^q*p*tau*tc/(gamma*m^p*tau + h*m^q*tc) + (beta + phi - 1)*gamma*h*m^p*m^(q - 1)*q*tau*tc/(gamma*m^p*tau + h*m^q*tc) - (gamma*m^(p - 1)*p*tau + h*m^(q - 1)*q*tc)*(beta + phi - 1)*gamma*h*m^p*m^q*tau*tc/(gamma*m^p*tau + h*m^q*tc)^2 + (delta*tau + 1)*k*m^(n - 1)*n*tc)/((tau*v + M)*m^nu)


model_out_growth <- function(temp,
                             temp_ref=15,
                             l,
                             lm,
                             Ea,
                             gamma=50,
                             delta=2,
                             phi=0.15,
                             h=30,
                             beta=0.2,
                             k=2,
                             p=0.8,
                             q=0.9,
                             n=0.8,
                             tmax=10,
                             slope=0.05,
                             tr = 1,
                             v=NULL,
                             dt=100){
  
 #browser()
  temps = length(temp)
  
  tc = exp(Ea*((temp+(288.2-temp_ref))-288.2)/(8.6173324*10^(-5)*(temp+(288.2-temp_ref))*288.2))
  
  ts <- seq(0,tmax,l=dt)
  
  withProgress(message = 'Calculating Winf', value = 0, {
    s <- array(0, c(temps,l,dt))
    allocs <- array(0, c(temps,l,dt))
    
    R0 <- array(0, c(temps,l,dt))
    surv <- array(0, c(temps,l,dt))
    
    s[,,1] <- min(wl(lm))
    
   
    dts <- (tmax/(dt-1))
    
    for(t in 2:dt) {
      tm1 <- get_taus(v,1,10,temp,s[,,t-1])
      Es <- (1-phi-beta)*(tm1*gamma*s[,,t-1]^p/(tm1*gamma*s[,,t-1]^p+tc*h*s[,,t-1]^q))*tc*h*s[,,t-1]^q-(1+ tm1*delta)*k*tc*s[,,t-1]^n 
      allocs[,,t] <- pmax(allocs[,,t-1],t(apply(lw(s[,,t-1]),1,inv_logit3,lm,ts[t],slope,tr)))
      
      s[,,t] <- s[,,t-1]+dts*(1-allocs[,,t])*Es
      
      surv[,,t] <- surv[,,t-1] + dts*(tm1*v$v+v$M)*s[,,t]^v$nu
      R0[,,t] <- R0[,,t-1] + dts*allocs[,,t]*Es*exp(-surv[,,t])
      
      
      incProgress(dts/tmax, detail = paste("Time", round(ts[t],2)))
    }
    ls <- lw(s)
    opt <- apply(R0[,,dt],1,function(x) ifelse(any(!is.nan(x)),which.max(x),NA))
    #browser()
    
    s=t(sapply(1:temps,function(x) s[x,opt[x],]))
    allocs=t(sapply(1:temps,function(x) allocs[x,opt[x],]))
    R00s=t(sapply(1:temps,function(x) R0[x,opt[x],]))
    
    lss <- reshape2::melt(s)
    colnames(lss) <- c('Temperature','t','size')
    lss$opt <- opt[lss$Temperature]
    lss$Temperature <- temp[lss$Temperature]
    lss$t <- ts[lss$t]
    
    
    alloc <- reshape2::melt(allocs)
    colnames(alloc) <- c('Temperature','t','allocs')
    alloc$opt <- opt[alloc$Temperature]
    alloc$Temperature <- temp[alloc$Temperature]
    alloc$t <- ts[alloc$t]
    
    R0s <- reshape2::melt(R00s)
    colnames(R0s) <- c('Temperature','t','R0')
    R0s$opt <- opt[R0s$Temperature]
    R0s$Temperature <- temp[R0s$Temperature]
    R0s$t <- ts[R0s$t]
    
    #browser()
    
    growth <- inner_join(inner_join(lss,alloc),R0s) %>% arrange(t,Temperature)
    
    
    winfs <- growth %>% group_by(Temperature) %>%
      summarise(l=size[ifelse(any(abs(allocs-0.5)<0.05),which.min(abs(allocs-0.5)),NA)-1],
                t=t[ifelse(any(abs(allocs-0.5)<0.05),which.min(abs(allocs-0.5)),NA)-1],
                opt=unique(opt))
    
    G <- rep(NA,length(unique(winfs$Temperature)))
    
    for(t in 2:(length(winfs$Temperature)-1)) {
      tau = winfs$Temperature[t]
      this.l <- winfs$l[t]
      if (is.na(this.l)) next
      tPM <- R0[t,opt[t],dt]
      lPM <- R0[t-1,opt[t],dt]
      nPM <- R0[t+1,opt[t],dt]
      
      sl <- abs((nPM-lPM)/(winfs$Temperature[t+1]-winfs$Temperature[t-1]))
      
      G[t] <- 0.04*this.l*(sl/tPM)
    }       
    
    
    #sG <- sign(G)
    winfs$G <- G/winfs$t ### why divide?
    winfs$L <- 10*(lw(winfs$l + winfs$G)-lw(winfs$l))
   
  })
  
  list(winfs=winfs, growth=growth)
  
}


model_out_growth_check <- function(temp,
                                   temp_ref=15,
                                   Ea,
                                   gamma=50,
                                   delta=2,
                                   phi=0.15,
                                   h=30,
                                   beta=0.2,
                                   k=2,
                                   p=0.8,
                                   q=0.9,
                                   n=0.8,
                                   mstar=1000,
                                   tmax=10,
                                   slope=0.05,
                                   tr = 1,
                                   v=NULL,
                                   dt=100,
                                   lm=NULL){
  

  temps = length(temp)
  
  tc = exp(Ea*((temp+(288.2-temp_ref))-288.2)/(8.6173324*10^(-5)*(temp+(288.2-temp_ref))*288.2))
  
  ts <- seq(0,tmax,l=dt)
  
  withProgress(message = 'Calculating Norms', value = 0, {
    s <- array(0, c(temps,temps,dt))
    allocs <- array(0, c(temps,temps,dt))
     
    s[,,1] <-  min(wl(lm))
    
    dts <- (tmax/(dt-1))
    
    for(t in 2:dt) {
      tm1 <- sapply(1:temps,function(g) {
        w <- v
        w$gamma <- gamma[g]
        get_taus(w,1,10,temp,s[,g,t-1])
      })
      f <- tm1*gamma*s[,,t-1]^p/(tm1*gamma*s[,,t-1]^p+tc*h*s[,,t-1]^q)
      inp <- (1-phi-beta)*f*tc*h*s[,,t-1]^q
      out <- (1+ tm1*delta)*k*tc*s[,,t-1]^n 
      Es=inp-out
      allocs[,,t] <- pmax(allocs[,,t-1],t(apply(lw(s[,,t-1]),1,inv_logit3,mstar,ts[t],slope,tr)))
        
      s[,,t] <- s[,,t-1]+dts*(1-allocs[,,t])*Es
      incProgress(dts/tmax, detail = paste("Time", round(ts[t],2)))
    }
    ls <- lw(s)
    
    #browser()
    
    ref = which.min(abs(temp-temp_ref))
    
    gls=t(sapply(1:temps,function(x) ls[ref,x,]))
    gallocs=t(sapply(1:temps,function(x) allocs[ref,x,]))
    
    tls=t(sapply(1:temps,function(x) ls[x,ref,]))
    tallocs=t(sapply(1:temps,function(x) allocs[x,ref,]))
     
    lss <- reshape2::melt(tls)
    colnames(lss) <- c('Temperature','t','t_length')
    lss$Temperature <- temp[lss$Temperature]
    lss$Gamma <- gamma[ref]
    lss$t <- ts[lss$t]
    
    alloc <- reshape2::melt(tallocs)
    colnames(alloc) <- c('Temperature','t','allocs')
    alloc$Temperature <- temp[alloc$Temperature]
    alloc$Gamma <- gamma[ref]
    alloc$t <- ts[alloc$t]
    
    
    lgs <- reshape2::melt(gls)
    colnames(lgs) <- c('Gamma','t','g_length')
    lgs$Temperature <- temp[ref]
    lgs$Gamma <- gamma[lgs$Gamma]
    lgs$t <- ts[lgs$t]
    
    galloc <- reshape2::melt(gallocs)
    colnames(galloc) <- c('Gamma','t','allocs')
    galloc$Temperature <- temp[ref]
    galloc$Gamma <- gamma[galloc$Gamma]
    galloc$t <- ts[galloc$t]
    
    #browser()
  
  })
  
  list(t_growth = inner_join(lss,alloc) %>% arrange(t,Temperature),
       g_growth = inner_join(lgs,galloc) %>% arrange(t,Temperature))
  
}


O2_supply <- function(O2 = 1:100,O2crit=20,P50 = 40,Tmax=30,Topt=15,T,omega=1.870,delta=1038,n=0.8){
  
  level <- delta*((Tmax-T)/(Tmax-Topt))^omega*exp(-omega*(Tmax-T)/(Tmax-Topt))/exp(-omega)
  365*24*level*(1-exp(-(O2-O2crit)/(-(P50-O2crit)/log(0.5))))/1000
  
}

# plot(O2_supply(level=200,P50 = 10),t='l',xlab='Disolved O2',ylab='02 supply')

eval_tau_max_o2 <- function(f=O2_supply(),
                            omega = 0.4,
                            gamma=50,
                            delta=2,
                            h=30,
                            phi=0.15,
                            beta=0.25,
                            k=2,
                            p=0.8,
                            q=0.9,
                            n=0.8,
                            m=100){
  
  
  tau_max<-(m^(-n-p)*(sqrt((gamma*(-f)*m^(n+p)+gamma*k*omega*m^(n+p)+delta*k*h*omega*m^(n+q)+beta*gamma*h*omega*m^(p+q))^2-4*gamma*delta*k*omega*m^(n+p)*(k*h*omega*m^(n+q)-f*h*m^(n+q)))+gamma*f*m^(n+p)-gamma*k*omega*m^(n+p)-delta*k*h*omega*m^(n+q)-beta*gamma*h*omega*m^(p+q)))/(2*gamma*delta*k*omega)
  
  tau_max
  
}

O2_fact <- function(temp,Tref=15){
  
  exp(-0.01851*(temp-Tref))
  
}

eval_tau_max_temp <- function(f=O2_supply(),
                              Ea = 0.52,
                              temp=seq(5,10,l=100),
                              temp_ref=15,
                              omega = 0.4,
                              gamma=50,
                              delta=2,
                              h=30,
                              phi=0.15,
                              beta=0.25,
                              k=2,
                              p=0.8,
                              q=0.9,
                              n=0.8,
                              m=100){
  
  tc = exp(Ea*((temp+(288.2-temp_ref))-288.2)/(8.6173324*10^(-5)*(temp+(288.2-temp_ref))*288.2))
  
  tau_max<-(m^(-n-p)*(sqrt((gamma*(-f)*m^(n+p)+gamma*tc*k*omega*m^(n+p)+delta*k*tc^2*h*omega*m^(n+q)+beta*tc*gamma*h*omega*m^(p+q))^2-4*gamma*delta*tc*k*omega*m^(n+p)*(k*h*tc^2*omega*m^(n+q)-f*h*tc*m^(n+q)))+gamma*f*m^(n+p)-gamma*k*omega*tc*m^(n+p)-delta*k*h*tc^2*omega*m^(n+q)-beta*gamma*h*tc*omega*m^(p+q)))/(2*gamma*delta*tc*k*omega)
  
  tau_max
  
}

get_taus <- function(v,tau_uc,O2_in,temp_in,m=10^seq(0,6,l=1000)){
  #browser()
 
  O2_tcor <- O2_fact(temp_in,5)
  O2 = O2_supply(O2=10*O2_tcor,Topt=v$Topt,O2crit=v$O2crit,Tmax=v$temp[length(v$temp)],T=temp_in,delta=v$lO,omega=v$shape,P50=v$P50,n=v$n)
  max_tau <- eval_tau_max_temp(f=O2,
                               temp=temp_in,
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
                               m=m)
  tau <- eval_tau_eq_temp(temp=temp_in,
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
                          m=m,
                          M=v$M,
                          v=v$v)
  
  tau_max = pmin(tau,max_tau)
  tau_max[tau_max<0] <- 0
  tau_max[tau_max>1] <- 1
  
  tau_max
  
}
