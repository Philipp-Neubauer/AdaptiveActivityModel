require(dplyr)
require(purrr)
require(tidyr)
require(ggplot2)
require(cowplot)

lw <- function(w) (w/0.0029)^(1/3.3365)
sc <- function(x) x/max(x)

O2_plotfun <- function(etas,Tops){
  grid <- expand.grid(eta = etas,
                      Topt = Tops)
  
  O2fun <- function(x) {
    ts = with(x,((30-seq(1,30,l=1000))/(30-Topt))^eta*exp(-eta*(30-seq(1,30,l=1000))/(30-Topt)))
    ts = ts/max(ts)
  }
  

  g <- grid %>% 
    by_row(~O2fun(.x), .collate='rows') %>%
    group_by(.row) %>%
    mutate(Temperature=seq(1,30,l=1000)) 
  
  ggplot(g) + 
    geom_line(aes(x = Temperature, 
                  y = .out,
                  linetype=factor(eta))) + 
    theme_bw() + 
    facet_grid(Topt~.) +
    scale_linetype_discrete(bquote(eta)) + 
    ylab(paste(expression(O_2, supply))) # + 
  #scale_color_discrete(bquote(T[opt]), h.start = 180)
  
}

plot_data_temp <- function(v){
  
  out = eval_tau_eq(gamma=v$gamma,
                    delta=v$delta,
                    phi=v$phi,
                    h=v$h,
                    beta=v$beta,
                    k=v$k,
                    p=v$p,
                    q=v$q,
                    n=v$n,
                    m=10^seq(-1,5,l=100),
                    M=v$M,
                    v=v$v)
  
  pd <- data_frame(m=10^seq(-1,5,l=100))
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
  
  o2frame <- data_frame(O2=seq(1,10,l=100),tau_o2 = tau_o2)
  
  O2_tcor <- 8.310212
  f=O2_supply(10*O2_tcor,Topt=v$Topt,O2crit=v$O2crit,Tmax=v$temp[length(v$temp)],T=v$temp,delta=v$lO,omega=v$shape,P50=v$P50,n=v$n)
  
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
}

plot_data_growth <- function(v){
  
  tau_uc = eval_tau_eq(gamma=v$gamma,
                       delta=v$delta,
                       phi=v$phi,
                       h=v$h,
                       beta=v$beta,
                       k=v$k,
                       p=v$p,
                       q=v$q,
                       n=v$n,
                       m=10^seq(1,5,l=1000),
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
                m = 10^seq(1,5,l=1000)
  )
}

plot_data_growth_tO2 <- function(v,lm=10^seq(0,6,l=1000),n_int=20){
  
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
  n_int <- n_int
  O2_range <- seq(1,10,l=n_int)
  winfs <- data.frame(matrix(NA,n_int,5))
  winfs[,1] <- O2_range
  winfs[,2] <- seq(min(v$temp),max(v$temp),l=n_int)
  
  winfs[,3:5] <- bind_rows(parallel::mclapply(1:n_int, function(i){#parallel::mclapply(1:n_int, function(i){
    taus <- get_taus(v,tau_uc,winfs[i,1],winfs[i,2],m=lm)
    
    #if(all(taus$tau_max<0.05)) return(data_frame(winf=0,winf_o2=0,winf_uc=0))
    model_out_par(tau_max = taus$tau_max,
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
  },mc.cores = 4))
  
  colnames(winfs) <- c('O2','Temperature','Temp','Oxygen','Base')
  winfs
  
}


# Activity model


calc_feed <- function(gamma=40,
                      phi=10,
                      h=30,
                      m=100,
                      p=0.8,
                      q=0.9
){
  (gamma *phi*m^p/(gamma *phi*m^p + h*m^q)) *h*m^q
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
  
  #delta*(-k)*tc*m^n-((1-beta-phi)*gamma^2*tau*tc*h*m^(2*p+q))/(gamma*tau*m^p+tc*h*m^q)^2+((1-beta-phi)*gamma*tc*h*m^(p+q))/(gamma*tau*phi*m^p+tc*h*m^q)
  
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
                        m=100,
                        M=0.2,
                        v=1){
  
  
  tc=1
  -(M*delta*h*k*m^(n - p + q)*tc - h*k*m^(n - p + q)*tc*v - sqrt(-((beta - 1)*h*k*m^(n - 2*p + 3*q) + h*k*m^(n - 2*p + 3*q)*phi)*tc^2*v^2 - ((beta - 1)*delta*gamma*k*m^(n - p + 2*q) + delta*gamma*k*m^(n - p + 2*q)*phi)*M^2*tc + (((beta - 1)*delta*h*k*m^(n - 2*p + 3*q) + delta*h*k*m^(n - 2*p + 3*q)*phi)*tc^2 + (gamma*h*m^(-p + 3*q)*phi^2 + (beta - 1)*gamma*k*m^(n - p + 2*q) + (beta^2 - 2*beta + 1)*gamma*h*m^(-p + 3*q) + (2*(beta - 1)*gamma*h*m^(-p + 3*q) + gamma*k*m^(n - p + 2*q))*phi)*tc)*M*v)*h)/(M*delta*gamma*k*m^n - ((beta - 1)*gamma*h*m^q + gamma*h*m^q*phi + gamma*k*m^n)*v)
  
  
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
                             m=100,
                             M=0.2,
                             v=1){
  
  tc = exp(Ea*((temp+273.2)-288.2)/(8.6173324*10^(-5)*(temp+273.2)*288.2))
  
  -(M*delta*h*k*m^(n - p + q)*tc - h*k*m^(n - p + q)*tc*v - sqrt(-((beta - 1)*h*k*m^(n - 2*p + 3*q) + h*k*m^(n - 2*p + 3*q)*phi)*tc^2*v^2 - ((beta - 1)*delta*gamma*k*m^(n - p + 2*q) + delta*gamma*k*m^(n - p + 2*q)*phi)*M^2*tc + (((beta - 1)*delta*h*k*m^(n - 2*p + 3*q) + delta*h*k*m^(n - 2*p + 3*q)*phi)*tc^2 + (gamma*h*m^(-p + 3*q)*phi^2 + (beta - 1)*gamma*k*m^(n - p + 2*q) + (beta^2 - 2*beta + 1)*gamma*h*m^(-p + 3*q) + (2*(beta - 1)*gamma*h*m^(-p + 3*q) + gamma*k*m^(n - p + 2*q))*phi)*tc)*M*v)*h)/(M*delta*gamma*k*m^n - ((beta - 1)*gamma*h*m^q + gamma*h*m^q*phi + gamma*k*m^n)*v)
  
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
                     r=0.2,
                     M=0.2,
                     gamma=50,
                     delta=2,
                     phi=0.15,
                     h=30,
                     beta=0.2,
                     k=2,
                     p=0.8,
                     q=0.9,
                     n=0.8) ((delta*tau + 1)*k*m^n*tc + (beta + phi - 1)*h*m^q*tau*tc/(h*m^(-p + q)*tc/gamma + tau))*m^(-q)*(q - 1)/(tau*r + M) - ((delta*tau + 1)*k*m^(n - 1)*n*tc + (beta + phi - 1)*h*m^(q - 1)*q*tau*tc/(h*m^(-p + q)*tc/gamma + tau) + (beta + phi - 1)*h^2*m^(-p + 2*q - 1)*(p - q)*tau*tc^2/((h*m^(-p + q)*tc/gamma + tau)^2*gamma))*m^(-q + 1)/(tau*r + M)

model_out_par <- function(tau_max,
                          tau,
                          tau_uc,
                          tau_o2,
                          temp,
                          Ea,
                          M=0.2,
                          v=1,
                          gamma=50,
                          delta=2,
                          phi=0.15,
                          h=30,
                          beta=0.2,
                          k=2,
                          p=0.8,
                          q=0.9,
                          n=0.8,
                          m=10^seq(0,6,l=1000),
                          ret_winf=F){
  
  #browser()
  tc = exp(Ea*((temp+273.2)-288.2)/(8.6173324*10^(-5)*(temp+273.2)*288.2))
  
  
  f <- tau_max*gamma*m^p/(tau_max*gamma*m^p+tc*h*m^q)
  inp <- (1-phi-beta)*f*tc*h*m^q
  out <- (1+tau_max*delta)*k*tc*m^n 
  
  #browser()
  mm <- get_dPdm(m=m,tau=tau_max,tc=tc,r=v,M=M,gamma=gamma,delta=delta,
                 phi=phi,h=h,beta=beta,k=k,p=p,q=q,n=n)
  mtemp <- tau_max^v*M*m^(n-1)
  
  test <- abs(mm-1)
  cond <- order(test[test<1])[1]
  winf <- ifelse(!is.na(cond), m[test<1][cond],min(m))
  
  # if(!(min(m)==winf | max(m)==winf)) stopifnot(numDeriv::grad(get_dPdm,x=winf,tau=tau_max[which.min(m[mtemp>mm])],
  
  #          tcc=tc,r=r,M=M,v=v,gamma=gamma,delta=delta,
  #          phi=phi,h=h,beta=beta,k=k,p=p,q=q,n=n,method='simple')<0)
  # 
  #f <- tau_uc*gamma*m^p/(tau_uc*gamma*m^p+h*m^q)
  #inp_uc <- (1-phi-beta)*f*h*m^q
  #out_uc <- (1+tau_uc*delta)*k*m^n+ r*m
  
  mm_uc <- get_dPdm(tau_uc,tc=1,r=v,M=M,gamma=gamma,delta=delta,
                    phi=phi,h=h,beta=beta,k=k,p=p,q=q,n=n,m=m)
  m_uc <- tau_uc^v*M*m^(n-1)
  winf_uc <- min(m[which.min(abs(mm_uc-1))],max(m))
  
  # if(!(min(m)==winf_uc | max(m)==winf_uc))  stopifnot(numDeriv::grad(get_dPdm,x=winf_uc,tau=tau_uc[which.min(m[m_uc>mm_uc])],
  #                tcc=1,r=r,M=M,v=v,gamma=gamma,delta=delta,
  #                phi=phi,h=h,beta=beta,k=k,p=p,q=q,n=n,method='simple')<0)
  # 
  #f <- tau_o2*gamma*m^p/(tau_o2*gamma*m^p+h*m^q)
  #inp_o2 <- (1-phi-beta)*f*h*m^q
  #out_o2 <- (1+tau_o2*delta)*k*m^n+ r*m
  
  mm_o2 <- get_dPdm(tau_o2,tc=1,r=v,M=M,gamma=gamma,delta=delta,
                    phi=phi,h=h,beta=beta,k=k,p=p,q=q,n=n,m=m)
  m_o2 <- tau_o2^v*M*m^(n-1)
  winf_o2 <- min(m[which.min(abs(mm_o2-1))],max(m))
  
  # if(!(min(m)==winf_o2 | max(m)==winf_o2)) stopifnot(min(m)==winf_o2 | max(m)==winf_o2| numDeriv::grad(get_dPdm,x=winf_o2 ,tau=tau_uc[which.min(m[m_o2>mm_o2])],
  #                tcc=1,r=r,M=M,v=v,gamma=gamma,delta=delta,
  #                phi=phi,h=h,beta=beta,k=k,p=p,q=q,n=n,method='simple')<0)
  # 
  winf <- data_frame(winf=winf,winf_o2=winf_o2,winf_uc=winf_uc)
  winf_range <- range(winf)
  
  if(ret_winf) return(winf)
  
  ix <- m<(winf_range[2]*2) & m>(winf_range[1]/2)
  ix <- which(ix)[round(seq(1,sum(ix),l=100))]
  
  data_frame(Rate = c(mm_uc[ix],m_uc[ix],mm_o2[ix],m_o2[ix],mm[ix],mtemp[ix]),
             Term = rep(rep(c('Production','Mortality'),each=100),3),
             Limitation = rep(c('None','O_2','Temp and O_2'),each=200),
             m=m[ix])
}


O2_supply <- function(O2 = 1:100,O2crit=20,P50 = 40, Tmax=30,Topt=15,T,omega=1.870,delta=1038,n=0.8){

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

get_taus <- function(v,tau_uc,O2_in,temp_in,m=10^seq(0,6,l=1000)){
  #browser()
  O2 = O2_supply(O2=O2_in,Topt=v$Topt,O2crit=v$O2crit,Tmax=v$temp[length(v$temp)],T=v$Topt,delta=v$lO,omega=v$shape,P50=v$P50,n=v$n)
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
                           m=m)
  
  tau_o2 <- apply(cbind(tau_uc,tau_o2),1,min)
  
  tau_o2[tau_o2<0] <- 0
  tau_o2[tau_o2>1] <- 1
  
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
  
  tau_max = apply(cbind(tau,max_tau),1,min)
  tau_max[tau_max<0] <- 0
  tau_max[tau_max>1] <- 1
  
  list(tau_max=tau_max, tau_o2=tau_o2, tau=tau)
  
}
