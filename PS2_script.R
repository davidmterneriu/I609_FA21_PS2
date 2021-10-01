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


################################################################################
# Question 1: 
################################################################################
dt1=0.001
time_max=54
iter_count=time_max/dt1

b1=lv_solver(x0=3,y0=5,A=1,B=0.5,C=2,D=1,iterations = iter_count,dt=dt1,conserve = F)

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
         theta=atan((Prey-start_x)/(Predator-start_y)))


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

theta_roots[6]-theta_roots[2]
