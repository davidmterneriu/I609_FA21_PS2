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


lv_solver=function(x0,y0,A,B,C,D,dt,iterations){
  
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
  colnames(pop_df)<-c("time","Prey","Predator")
  rownames(pop_df)<-c()
  plot_sys=pop_df%>%
    gather(key="Species",value="pop",-1)%>%
    ggplot(aes(x=time,y=pop,color=Species))+
    geom_line()+
    theme_bw()+
    labs(x="Time",y="Population")+
    scale_color_manual(values = c("Red","Blue"))+
    scale_x_continuous(breaks = scales::pretty_breaks())
  return(list(pop_df=pop_df,equilib=equilib, plot=plot_sys))
}

################################################################################
# Question 1: 
################################################################################
dt1=0.001
iter_count=12/dt1

b1=lv_solver(x0=3,y0=5,A=10,B=0.5,C=2,D=1,iterations = iter_count,dt=dt1)

b1$plot
#plotly::ggplotly()

b1$pop_df%>%ggplot(aes(x=time))+
  geom_line(aes(y=prey),color="blue")+
  geom_line(aes(y=predator),color="red")
