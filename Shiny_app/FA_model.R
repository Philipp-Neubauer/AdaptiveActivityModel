# Activity model


calc_feed <- function(gamma=calc_gamma()*10,
                      phi=10,
                      h=30,
                      m=100,
                      p=0.8,
                      q=0.9
){
  gamma *phi*m^p/(gamma *phi*m^p + h*m^q) 
}


calc_C02_demand <- function(feed,
                            h=30,
                            q=0.9,
                            alpha=0.4,
                            beta=0.6,
                            omega=0.4,
                            k=5,
                            n=0.8,
                            m=100
){
  omega*(alpha*beta*feed*h*m^q+k*m^n)
}

eval_dtau <- function(tau = seq(0,1,l=100),
                      gamma=50,
                      delta=2,
                      h=30,
                      beta=0.75,
                      omega=0.4,
                      k=2,
                      p=0.8,
                      q=0.9,
                      n=0.8,
                      m=100){
  
  delta*(-k)*m^n-((1-beta-phi)*gamma^2*tau*h*m^(2*p+q))/(gamma*tau*m^p+h*m^q)^2+((1-beta-phi)*gamma*h*m^(p+q))/(gamma*tau*m^p+h*m^q)
  
}

eval_dtau_temp <- function(Ea,
                           temp,
                           tau = seq(0,1,l=100),
                           gamma=50,
                           delta=2,
                           h=30,
                           beta=0.75,
                           omega=0.4,
                           k=2,
                           p=0.8,
                           q=0.9,
                           n=0.8,
                           m=100){
  
  tc = exp(Ea*((temp+273.2)-288.2)/(8.6173324*10^(-5)*(temp+273.2)*288.2))
  
  delta*(-k)*tc*m^n-((1-beta-phi)*gamma^2*tau*tc*h*m^(2*p+q))/(gamma*tau*m^p+tc*h*m^q)^2+((1-beta-phi)*gamma*tc*h*m^(p+q))/(gamma*tau*phi*m^p+tc*h*m^q)
  
}

eval_tau_eq <- function(gamma=50,
                        delta=2,
                        phi=10,
                        h=30,
                        beta=0.75,
                        k=2,
                        p=0.8,
                        q=0.9,
                        n=0.8,
                        m=100){
  
  (m^(-n-2*p)*(gamma*(-h)*delta*k*m^(n+p+q)+sqrt((1-beta-phi)*gamma^3*h^2*delta*k*m^(n+3*p+2*q))))/(gamma^2*delta*k)
  
}

eval_tau_eq_temp <- function(Ea,
                             temp,
                             gamma=50,
                             delta=2,
                             phi=10,
                             h=30,
                             beta=0.75,
                             k=2,
                             p=0.8,
                             q=0.9,
                             n=0.8,
                             m=100){
  
  tc = exp(Ea*((temp+273.2)-288.2)/(8.6173324*10^(-5)*(temp+273.2)*288.2))
  
  (m^(-n-2*p)*(gamma*(-h)*tc*delta*k*m^(n+p+q)+(sqrt(gamma^3*delta*k*tc*h^2*m^(n+3*p+2*q)*(1-beta-phi)))))/(gamma^2*delta*k)
  
}

model_out <- function(tau_max,
                      temp,
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
  
  tc = exp(Ea*((temp+273.2)-288.2)/(8.6173324*10^(-5)*(temp+273.2)*288.2))
  
  f <- tau_max*gamma*m^p/(tau_max*gamma*m^p+tc*h*m^q)
  inp <- (1-phi-beta)*f*tc*h*m^q
  out <- (1+tau_max*delta)*k*tc*m^n
  e <-  inp -out
  efficiency <- e/f
  efficiency[efficiency<0] <- 0
  predation_rate <- f/phi
  
  met = beta*f*tc*h*m^q+ ((1+tau_max*delta)*k*tc*m^n)
  data.frame(feeding = f, 
             input = inp,
             out = out,
             energy = e, 
             efficiency = efficiency, 
             predation=predation_rate, 
             metabolism = met)
}


model_out_par <- function(tau_max,
                          tau_uc,
                          tau_o2,
                          temp,
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
                          m=10^seq(1,5,l=100),
                          ret_winf=F){
  
  #browser()
  tc = exp(Ea*((temp+273.2)-288.2)/(8.6173324*10^(-5)*(temp+273.2)*288.2))
  
  f <- tau_max*gamma*m^p/(tau_max*gamma*m^p+tc*h*m^q)
  inp <- (1-phi-beta)*f*tc*h*m^q
  out <- (1+tau_max*delta)*k*tc*m^n + r*m
  winf <- min(m[out>inp],max(m))
  
  f <- tau_uc*gamma*m^p/(tau_uc*gamma*m^p+h*m^q)
  inp_uc <- (1-phi-beta)*f*h*m^q
  out_uc <- (1+tau_uc*delta)*k*m^n+ r*m
  winf_uc <- min(m[out_uc>inp_uc],max(m))
  
  f <- tau_o2*gamma*m^p/(tau_o2*gamma*m^p+h*m^q)
  inp_o2 <- (1-phi-beta)*f*h*m^q
  out_o2 <- (1+tau_o2*delta)*k*m^n+ r*m
  winf_o2 <- min(m[out_o2>inp_o2],max(m))
  
  winf <- c(winf,winf_o2,winf_uc)
  winf_range <- range(winf)
  
  if(ret_winf) return(winf)
  
  ix <- m<(winf_range[2]*2) & m>(winf_range[1]/2)
  ix <- which(ix)[round(seq(1,sum(ix),l=100))]
  
  data.frame(Carbon = c(inp_uc[ix],out_uc[ix],inp_o2[ix],out_o2[ix],inp[ix],out[ix]),
             Term = rep(rep(c('Gain','Loss'),each=100),3),
             Limitation = rep(c('None','O_2','Temp and O_2'),each=200),
             m=m[ix])
}


O2_supply <- function(O2 = 1:100,O2crit=20,P50 = 40, Tmax=30,Topt=15,T,omega=1.870,delta=1038,n=0.8){
  
  level <- delta*((Tmax-T)/(Tmax-Topt))^omega*exp(-omega*(Tmax-T)/(Tmax-Topt))
  365*24*level*(1-exp(-(O2-O2crit)/(-(P50-O2crit)/log(0.5))))/(1000*100^0.8)
  
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
  
  tc = exp(Ea*((temp+273.2)-288.2)/(8.6173324*10^(-5)*(temp+273.2)*288.2))
  
  tau_max<-(m^(-n-p)*(sqrt((gamma*(-f)*m^(n+p)+gamma*tc*k*omega*m^(n+p)+delta*k*tc^2*h*omega*m^(n+q)+beta*tc*gamma*h*omega*m^(p+q))^2-4*gamma*delta*tc*k*omega*m^(n+p)*(k*h*tc^2*omega*m^(n+q)-f*h*tc*m^(n+q)))+gamma*f*m^(n+p)-gamma*k*omega*tc*m^(n+p)-delta*k*h*tc^2*omega*m^(n+q)-beta*gamma*h*tc*omega*m^(p+q)))/(2*gamma*delta*tc*k*omega)
  
  tau_max
  
}

get_taus <- function(v,tau_uc,O2_in,temp_in){
  #browser()
  O2 = O2_supply(O2=O2_in,Topt=v$Topt,O2crit=v$O2crit,Tmax=v$temp[length(v$temp)],T=v$Topt,delta=v$lO,omega=1.8,P50=v$P50,n=v$n)
  tau_o2 = eval_tau_max_o2(f=O2,omega = v$omega,
                                  gamma=v$gamma,
                                  delta=v$delta,
                                  phi=v$phi,
                                  h=v$h,
                                  beta=v$beta,
                                  k=v$k,
                                  p=v$p,
                                  q=v$q,
                                  n=v$n,
                                  m=10^seq(1,5,l=10000))
  
  tau_o2 <- apply(cbind(tau_uc,tau_o2),1,min)
  
  tau_o2[tau_o2<0] <- 0
  tau_o2[tau_o2>1] <- 1
  
  O2_tcor <- O2_fact(temp_in,5)
  O2 = O2_supply(O2=10*O2_tcor,Topt=v$Topt,O2crit=v$O2crit,Tmax=v$temp[length(v$temp)],T=temp_in,delta=v$lO,omega=1.8,P50=v$P50,n=v$n)
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
                               m=10^seq(1,5,l=10000))
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
                          m=10^seq(1,5,l=10000))
  
  tau_max = apply(cbind(tau,max_tau),1,min)
  tau_max[tau_max<0] <- 0
  tau_max[tau_max>1] <- 1
  
  list(tau_max=tau_max, tau_o2=tau_o2)
  
}

