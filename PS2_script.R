#I608 PS2 Work 
#Lotka-Volterra Eqns
rm(list=ls())

library(readr)
library(readxl)
library(tidyverse)
library(ggplot2)
library(tidyverse)
library(haven)
library(lmtest)
library(sandwich)
library(ggquiver)
library(latex2exp)
library(rootSolve)
library(scales)


lv_solver=function(x0,y0,A,B,C,D,dt,iterations,conserve=F){
  
  equilib<-data.frame(number=c(1,2),
                      x=c(0,C/D),
                      y=c(0,A/B)) 
  
  x=x0
  y=y0
  pop_df<-c(0,x,y)
  for(i in 1:iterations){
    dx<-(A-B*y)*x
    dy<-(D*x-C)*y
    x_new<-x+dx*dt
    y_new<-y+dy*dt
    new_df=c(i*dt,x_new,y_new)
    pop_df=rbind(pop_df,new_df)
    x<-x_new
    y<-y_new
  }
  pop_df=as.data.frame(pop_df)
  
  consv_energy=exp(A*log(y0)-B*y0+C*log(x0)-D*x0)
  
  
  
  colnames(pop_df)<-c("time","Prey","Predator")
  rownames(pop_df)<-c()
  
  pop_df=pop_df%>%mutate(slope1=(-C+D*Prey)*Predator/(A-B*Predator)*Prey,
                         ene=A*log(Predator)-B*Predator+C*log(Prey)-D*Prey,
                         ene=exp(ene),
                         cons_ene=ene/consv_energy)
  if(conserve){
    pop_df=pop_df%>%
      mutate(Prey=Prey*cons_ene,
             Predator=Predator*cons_ene)
  }
  
  plot_sys=pop_df%>%
    gather(key="Species",value="pop",-c("time","ene","cons_ene","slope1"))%>%
    ggplot(aes(x=time,y=pop,color=Species))+
    geom_line()+
    theme_bw()+
    labs(x="Time",y="Population")+
    scale_color_manual(values = c("Red","Blue"))+
    scale_x_continuous(breaks = scales::pretty_breaks())
  return(list(pop_df=pop_df,equilib=equilib, plot=plot_sys))
}

num_deriv_fun=function(x,f_x){
  #browser()
  eps1=.Machine$double.eps
  fun1=approxfun(x,f_x)
  df_dx=(fun1(x+eps1)-fun1(x-eps1))/(2*eps1)
  return(df_dx)
}

lv_solver_feild=function(A,B,C,D){
  equilib<-data.frame(number=c(1,2),
                      x=c(0,C/D),
                      y=c(0,A/B)) 
  x_eq=C/D
  y_eq=A/B
  
  #lim1=max(x_eq,y_eq)
  
  x_seq=seq(0.1*x_eq,2*x_eq,length.out=25)
  y_seq=seq(0.1*y_eq,2*y_eq,length.out=25)
  
  graph_df=expand.grid(x=x_seq,y=y_seq)%>%
    as.data.frame()%>%mutate(dx=(A-B*y)*x,
                             dy=(D*x-C)*y,
                             xend=x+dx/sd(dx),
                             yend=y+dy/sd(dy),
                             K=y^A*exp(-B*y)*x^C*exp(-D*x))
  mes1=paste0("Note: (A=",A,"; B=",B,"; C=",C,"; D=",D,")")
  mes2=paste("(Equilibrium prey/predator(s): ",x_eq, " and ", y_eq,")")
  
  g1=graph_df%>%ggplot(aes(x,y,fill=K))+
    geom_raster(interpolate = TRUE)+
    geom_segment(aes(xend=xend,yend=yend),
                 color="black",
                 arrow = arrow(length = unit(0.01, "npc")))+
    scale_fill_viridis_c()+
    xlim(0,2*x_eq)+
    ylim(0,2*y_eq)+
    ggthemes::theme_tufte()+
    labs(x="Prey",y="Predator",fill="Energy",
         title="Lotka–Volterra Phase Diagram",
         subtitle = mes2,
         caption=mes1)+
    geom_point(x=x_eq,y=y_eq,color="red")
  
  
  return(g1)
}


################################################################################
# Testing Functions
################################################################################
if(FALSE){

dt1=0.001
time_max=54
iter_count=time_max/dt1

A1=1
B1=0.5
C1=2
D1=1

b1=lv_solver(x0=3,y0=5,A=A1,B=B1,C=C1,D=D1,iterations = iter_count,dt=dt1,conserve = F)

#b1$pop_df%>%ggplot(aes(x=time,y=slope1))+
#  geom_path()

#b1$pop_df$time[which.max(b1$pop_df$Predator)]

b1$pop_df%>%
  ggplot(aes(x=Prey,y=Predator))+
  geom_path()+
  geom_point(data=b1$equilib[2,],aes(x=x,y=y))

cycle_df=b1$pop_df%>%
  mutate(start_x=b1$equilib[2,"x"],
         start_y=b1$equilib[2,"y"],
         theta=atan((Prey-start_y)/(Predator-start_x)))


cycle_df=cycle_df%>%mutate(slope=(Predator-start_y)/(Prey-start_x))

slope_fun=approxfun(x=cycle_df$time,y=cycle_df$slope-cycle_df$slope[1])
rootSolve::uniroot.all(f=slope_fun,interval = c(min(cycle_df$time),max(cycle_df$time)))


theta_fun=approxfun(x=cycle_df$time,y=cycle_df$theta)


theta_roots=rootSolve::uniroot.all(theta_fun,interval = c(min(cycle_df$time),max(cycle_df$time)))
theta_df=data.frame(id=seq_along(theta_roots),time=theta_roots,theta=0)
theta_df$id<-(theta_df$id-1) %%4

theta_df$group_id=0
counter=0
for(i in 2:nrow(theta_df)){
  if(theta_df$id[i]==0){
    counter=counter+1
  }
  theta_df$group_id[i]<-counter
}



theta_df%>%group_by(id)%>%
  mutate(dif=time-lag(time))%>%
  filter(is.na(dif)==F)%>%ungroup()%>%
  group_by(group_id)%>%
  tally()

theta_dif=theta_df%>%group_by(id)%>%
  mutate(dif=time-lag(time))%>%
  filter(is.na(dif)==F)%>%ungroup()

theta_dif%>%
  ggplot(aes(x=time,y=dif))+
  geom_line()+
  geom_point(aes(color=as.factor(id)))+
  labs(x="Time",y="Estimated Period",color="Point")





period_finder=function(guess_T,C,D){
  r1=integrate(f=prey_function,lower = 0,upper = guess_T)$value
  r1=r1/guess_T
  avg_pop=D/C
  dif=abs(r1-avg_pop)
  return(dif)
}
period_finder=Vectorize(period_finder,vectorize.args = "guess_T")

uniroot(period_finder,lower=0.1,upper=5,C=C1,D=D1)

period_finder(guess_T =seq(0,5,by=0.01),C=C1,D=D1)%>%plot()

nlm(f=period_finder,p=1,C=C1,D=D1)


mod1=lm(data=theta_dif,dif~time)
summary(mod1)

pi_scales <- math_format(.x * pi, format = function(x) x / pi)

cycle_df%>%
  filter(time<5.2)%>%
  ggplot()+
  geom_line(aes(x=time,y=theta))+
  geom_point(data=filter(theta_df,time<5.2),
             aes(x=time,y=theta,color=as.factor(round(time,3))),
             size=3)+
  scale_y_continuous(labels = pi_scales,
                     breaks = seq(-pi / 2,  pi / 2, pi/4))+
  labs(x="Time",y="Theta",color="Root",title="Theta Time Series")+
  theme_minimal()


index_find=data.frame()
for(i in 1:length(theta_df$time)){
  index_find=rbind(index_find,
                   data.frame(id=which.min(abs(theta_df$time[i]-cycle_df$time))))
}


cycle_df%>%
  filter(time<5.2)%>%
  ggplot()+
  geom_path(aes(x=Prey,y=Predator))+
  geom_segment(data=cycle_df[index_find$id,]%>%filter(time<5.2),
               aes(x=start_x,xend=Prey,y=start_y,yend=Predator,
                   color=as.factor(time)))+
  geom_point(data=cycle_df[index_find$id,]%>%filter(time<5.2),
             aes(x=Prey,y=Predator,color=as.factor(time)),
             size=2)+
  labs(color="Root",title="Phase Diagram")+
  theme_minimal()+
  geom_point(data=b1$equilib[2,],aes(x=x,y=y),size=2)





plotly::ggplotly()

b1$pop_df%>%ggplot(aes(x=time))+
  geom_line(aes(y=prey),color="blue")+
  geom_line(aes(y=predator),color="red")



p1=cycle_df%>%
  filter(time>=theta_roots[2],time<=theta_roots[6])%>%
  ggplot()+
  geom_path(aes(x=Prey,y=Predator))+
  geom_segment(aes(x=start_x,xend=Prey,y=start_y,yend=Predator))+
  labs(title="Time: {frame_time}")+
  transition_time(time)
p1  

period_estimate=theta_roots[6]-theta_roots[2]

prey_function=approxfun(x=cycle_df$time,y=cycle_df$Prey)
pred_function=approxfun(x=cycle_df$time,y=cycle_df$Predator)

time_test=seq(0,period_estimate,length.out=100)

period_test=data.frame(time0=time_test,time1=time_test+period_estimate)

period_test%>%
  mutate(X0=prey_function(time0),
         X1=prey_function(time1),
         Y0=pred_function(time0),
         Y1=pred_function(time1))%>%
  gather(key="series",value=pop,-c(time0,time1))%>%
  mutate(species=ifelse(grepl("X",series),"Prey","Predator"))%>%
  ggplot(aes(x=time0,y=pop,color=species))+
  geom_line()+
  facet_wrap(~series,scales = "free")+
  labs(x="Time",y="Pop",color="Species",title="Periodic Populations")+
  scale_color_manual(values = c("Red","Blue"))+
  ggthemes::theme_clean()}


################################################################################
# PS Questions
################################################################################
dt1=0.001
time_max=12
iter_count=time_max/dt1

A1=1
B1=0.5
C1=2
D1=1

q1_a=lv_solver(x0=3,y0=5,A=A1,B=B1,C=C1,D=D1,iterations = iter_count,dt=dt1,conserve = F)

A2=1.5

q2_a=lv_solver(x0=3,y0=5,A=A2,B=B1,C=C1,D=D1,iterations = iter_count,dt=dt1,conserve = F)



# Populations x(12) and y(10)
q1_a$pop_df%>%select(time,Prey,Predator)%>%
  filter(time %in% c(10,12))%>%round(3)

q2_a$pop_df%>%select(time,Prey,Predator)%>%
  filter(time %in% c(10,12))%>%round(3)

#Plots 

q1_a$plot+scale_x_continuous(limits = c(0,10),breaks = seq(0,10,1))+
  scale_y_continuous(limits = c(0,6),breaks = seq(0,6,1))+
  labs(title = "Prey x(t) and Predator y(t) Populations",
       subtitle = "Note: a=1")+
  geom_hline(yintercept = 0,color="purple",linetype="dashed")

q2_a$plot+scale_x_continuous(limits = c(0,10),breaks = seq(0,10,1))+
  scale_y_continuous(limits = c(0,6),breaks = seq(0,6,1))+
  labs(title = "Prey x(t) and Predator y(t) Populations",
       subtitle = "Note: a=1.5")+
  geom_hline(yintercept = 0,color="purple",linetype="dashed")


#Period testing 

dt1=0.001
time_max=100
iter_count=time_max/dt1

tictoc::tic()
q2=lv_solver(x0=3,y0=5,A=A1,B=B1,C=C1,D=D1,iterations = iter_count,dt=dt1,conserve = F)
tictoc::toc()

q2$pop_df%>%
  ggplot(aes(x=time,y=cons_ene))+
  geom_line(color="darkgreen",size=1.2)+
  theme_bw()+
  labs(x="Time",y="K(t)/K(0)",title="Dissipating K over Time")


simple_time=lm(data=q2$pop_df,cons_ene~time)
summary(simple_time)


q2$pop_df%>%
  ggplot(aes(x=Prey,y=Predator,color=ene))+
  geom_path()+
  geom_point(aes(x=2,y=2),color="red",size=3)+
  viridis::scale_color_viridis() + theme_bw()+
  labs(x="Prey",y="Predator",color="Law of Motion (K)",
       caption="Note: Equilibrium prey/redator population is (2,2)",
       title="Phase-Space Diagram: Dissipating K")

q2$pop_df[1,]%>%round(3)


tictoc::tic()
q3=lv_solver(x0=3,y0=5,A=A2,B=B1,C=C1,D=D1,iterations = iter_count,dt=dt1,conserve = F)
tictoc::toc()



cycle_df_a1=q2$pop_df%>%
  mutate(start_x=q2$equilib[2,"x"],
         start_y=q2$equilib[2,"y"],
         theta=atan((Prey-start_y)/(Predator-start_x)))


cycle_df_a2=q3$pop_df%>%
  mutate(start_x=q3$equilib[2,"x"],
         start_y=q3$equilib[2,"y"],
         theta=atan((Prey-start_y)/(Predator-start_x)))



theta_fun1=approxfun(x=cycle_df_a1$time,y=cycle_df_a1$theta)
theta_fun2=approxfun(x=cycle_df_a2$time,y=cycle_df_a2$theta)


theta_roots1=rootSolve::uniroot.all(theta_fun1,interval = c(min(cycle_df_a1$time),max(cycle_df_a1$time)),n=1000)
theta_roots2=rootSolve::uniroot.all(theta_fun2,interval = c(min(cycle_df_a2$time),max(cycle_df_a2$time)),n=1000)


theta_df1=data.frame(id=seq_along(theta_roots1),time=theta_roots1,theta=0)
theta_df1$id<-(theta_df1$id-1) %%4

theta_df2=data.frame(id=seq_along(theta_roots2),time=theta_roots2,theta=0)
theta_df2$id<-(theta_df2$id-1)%%4


index_find=data.frame()
for(i in 1:length(theta_df1$time)){
  index_find=rbind(index_find,
                   data.frame(id=which.min(abs(theta_df1$time[i]-cycle_df_a1$time))))
}


cycle_df_a1%>%
  filter(time<5.2)%>%
  ggplot()+
  geom_path(aes(x=Prey,y=Predator))+
  geom_segment(data=cycle_df_a1[index_find$id,]%>%filter(time<5.2),
               aes(x=start_x,xend=Prey,y=start_y,yend=Predator,
                   color=as.factor(time)),size=1.2)+
  geom_point(data=cycle_df_a1[index_find$id,]%>%filter(time<5.2),
             aes(x=Prey,y=Predator,color=as.factor(time)),
             size=4)+
  labs(color="Root",title="Phase Diagram")+
  theme_minimal()+
  geom_point(data=q2$equilib[2,],aes(x=x,y=y),size=4)

pi_scales <- math_format(.x * pi, format = function(x) x / pi)

cycle_df_a1%>%
  filter(time<5.2)%>%
  ggplot(aes(x=time,y=theta))+
  geom_line()+
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_point(data=filter(theta_df1,time<5.2),
             aes(x=time,y=theta,color=as.factor(round(time,3))),
             size=3)+
  ggrepel::geom_label_repel(data=filter(theta_df1,time<5.2),
                            aes(x=time,y=theta,label=paste0("Root:",id)))+
  scale_y_continuous(labels = pi_scales,
                     breaks = seq(-pi / 2,  pi / 2, pi/4))+
  labs(x="Time",y=TeX("$\\Theta(t)$"),color="Root",title="Theta Time Series")+
  theme_minimal()



theta_dif1=theta_df1%>%
  group_by(id)%>%
  mutate(period=time-lag(time))%>%
  filter(is.na(period)==F)%>%
  ungroup()

theta_dif2=theta_df2%>%
  group_by(id)%>%
  mutate(period=time-lag(time))%>%
  filter(is.na(period)==F)%>%
  ungroup()


rbind(theta_dif1%>%mutate(model="a=1"),
      theta_dif2%>%mutate(model="a=1.5"))%>%
  ggplot(aes(x=time,y=period))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~model,scales = "free")+
  labs(x="Time",y="Estimated Period")+
  theme_bw()

period_mod1=lm(data=theta_dif1,period~time)
period_mod2=lm(data=theta_dif2,period~time)

texreg::texreg(list(period_mod1,period_mod2),digits = 4,
               stars = c(0.01,0.05,0.1))

#Double Time
(coef(period_mod1)[1]/coef(period_mod1)[2]%>%as.numeric())%>%round(0)
(coef(period_mod2)[1]/coef(period_mod2)[2]%>%as.numeric())%>%round(0)

#Max pop during first period
cycle_df_a1%>%filter(as.numeric(coef(period_mod1)[1])>time)%>%
  summarize(max_prey=max(Prey),
            max_pred=max(Predator),
            min_prey=min(Prey),
            min_pred=min(Predator))%>%
  round(3)



cycle_df_a2%>%filter(as.numeric(coef(period_mod2)[1])>time)%>%
  summarize(max_prey=max(Prey),
            max_pred=max(Predator),
            min_prey=min(Prey),
            min_pred=min(Predator))%>%
  round(3)





cycle_df_a1%>%select(time,Prey,Predator)%>%
  mutate(Predator=cummax(Predator),
         Prey=cummax(Prey))%>%
  filter(time>coef(period_mod1)[1])%>%
  gather(key=Species,value="pop",-1)%>%
  ggplot(aes(x=time,y=pop,color=Species))+
  geom_line(size=1.2)+
  labs(x="Time",y="Max Pop")+
  scale_color_manual(values = c("Red","Blue"))+
  theme_minimal()+
  theme(legend.position = "top")



cycle_df_a1%>%select(time,Prey,Predator)%>%
  mutate(Predator=cummin(Predator),
         Prey=cummin(Prey))%>%
  filter(time>coef(period_mod1)[1])%>%
  gather(key=Species,value="pop",-1)%>%
  ggplot(aes(x=time,y=pop,color=Species))+
  geom_line(size=1.2)+
  labs(x="Time",y="Min Pop")+
  scale_color_manual(values = c("Red","Blue"))+
  theme_minimal()+
  theme(legend.position = "top")


cycle_df_a2%>%select(time,Prey,Predator)%>%
  mutate(Predator=cummax(Predator),
         Prey=cummax(Prey))%>%
  filter(time>coef(period_mod2)[1])%>%
  gather(key=Species,value="pop",-1)%>%
  ggplot(aes(x=time,y=pop,color=Species))+
  geom_line(size=1.2)+
  labs(x="Time",y="Max Pop")+
  scale_color_manual(values = c("Red","Blue"))+
  theme_minimal()+
  theme(legend.position = "top")


cycle_df_a2%>%select(time,Prey,Predator)%>%
  mutate(Predator=cummin(Predator),
         Prey=cummin(Prey))%>%
  filter(time>coef(period_mod2)[1])%>%
  gather(key=Species,value="pop",-1)%>%
  ggplot(aes(x=time,y=pop,color=Species))+
  geom_line(size=1.2)+
  labs(x="Time",y="Min Pop")+
  scale_color_manual(values = c("Red","Blue"))+
  theme_minimal()+
  theme(legend.position = "top")


#Extinction Event

a1_min=cycle_df_a1%>%select(time,Prey,Predator)%>%
  mutate(Predator=cummin(Predator),
         Prey=cummin(Prey))%>%
  filter(time>coef(period_mod1)[1])
 
a2_min=cycle_df_a2%>%select(time,Prey,Predator)%>%
  mutate(Predator=cummin(Predator),
         Prey=cummin(Prey))%>%
  filter(time>coef(period_mod1)[1])


extinct_mod1=lm(data=a1_min,Predator~time)
extinct_mod2=lm(data=a2_min,Predator~time)
