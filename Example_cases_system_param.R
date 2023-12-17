
# Libraries necessary to run the project
library(RxODE)
library(reshape2)
library(gridExtra)
library(ggplot2)
library(stats)
library(MASS)
library(plyr)
library(dplyr)
library(grid)
library(Metrics)
library(ODEsensitivity)
library(data.table)   
library(tidyr)
library(Phxnlme)
library(forecast)



##################                                             ####################
##################      SA  and volume in the apical side      ####################
##################                                             ####################

ini.Conc=10
gut.liver.pop <-  RxODE({         # General model
  # Volume of the cells
  V5 = 1.7    # Typical enterocyte volume of 1 Mio of cells
  V6 = 3.9    # Typical enterocyte volume of 1 Mio of cells
  V2=V2
  Q =Med_Q
  ### Concentration/amount relationship for the substrate
  C2.S = A2.S / V2               # Substrate concentration in the AP
  C3.S = A3.S / V3               # Substrate concentration in the BS
  C4.S = A4.S  / V4               # Substrate concentration in the LMPS
  C5.S = A5.S  / (N3 * V5)       # Substrate concentration in the ent.
  C6.S = A6.S  / (N4 * V6)       # Substrate concentration in the hep.
  
  
  CL23 = SA23 * convfactor*Papp*exp(eta.Papp)
  
  
  
  # Substrate
  d/dt(A2.S) = - CL23 * NPERM * (C2.S * fum2 - C5.S) # AP. passive diffusion ent.
  d/dt(A3.S) = CL23 * NPERM * (C5.S- C3.S * fum)     # BS. passive diffusion
  + Q  * C4.S 
  - Q  * C3.S                 # Circulation flow out (BS -> LIVER MPS)     
  d/dt(A4.S) =  Q * C3.S                      # Circulation flow in (MC -> LMPS)
  - Q * C4.S                       # Circulation flow out (LMPS -> MC)        
  - CL56S * (C4.S * fum - C6.S)            # Passive diffusion hep.
  d/dt(A5.S) = CL23 * NPERM * (C2.S * fum2 - C5.S)   # ENT. passive diffusion .ent apical side
  - CL23 * NPERM *(C5.S - C3.S * fum)     # passive diffusion .ent basolateral side
  - CL3tot * N3 * C5.S                     # Intestinal metabolism
  d/dt(A6.S) = CL56S * (C4.S * fum - C6.S)            # HEP. passive diffusion hep.
  - CL4tot * N4 * C6.S                     # Metabolism in the liver cells
  
  
  C2ObsS = C2.S *exp(CEps.C2.S) 
  C3ObsS = C3.S *exp(CEps.C3.S) 
  C4ObsS = C3.S *exp(CEps.C4.S) 
  
})

SA23_i= c(0.143,0.33,1.12,4.67)
V2_i= c(75,100,500,1500)
#V2_i= c(325,325,325,325)

sim = data.frame()
for (i in 1:length(SA23_i)){
  
  #Parameters
  param=c( V2=V2_i[i]  ,V3=1506 , V4 =1394, fum =  1,fum2 = 1, fumM =  1,fum2M = 1, NPERM=2,Papp=0.4,CL3tot=0, CL4tot=0,CL56S=10000000, SA23=SA23_i[i],N3=0.5,N4=0.5,convfactor=0.006, Med_Q=60  )
  
  #Event table                           
  qd <- et(seq(0, 4320, by = 1))
  
  #Ini. cond                          
  # Initialization 
  initialization <- c(A2.S = ini.Conc* as.numeric(param["V2"]), A3.S = 0)
  
  #
  omega.matrix= lotri({c(eta.Papp ~  0.00001) })
  
  sigma.matrix <- lotri(CEps.C2.S ~0.000015^2, CEps.C3.S ~ 0.00000015^2, CEps.C4.S ~ 0.000015^2)
  ### Sensitivity analysis ####
  set.seed(120)
  sim_i <- rxSolve(object=gut.liver.pop, params =param,events =qd,inits=initialization, omega=omega.matrix, sigma=sigma.matrix, sub=1)
  
  sim_i = as.data.frame(sim_i)
  
  sim_i$SA23= SA23_i[i]
  
  sim= rbind(sim, sim_i)
}



sim_samples = subset(sim, time %in% c(0,30,60,120,240,480,1440,2880,4320))


p1=plot(ggplot(sim_samples,aes(time/60,C2ObsS)) +  geom_line(aes(color = as.factor(SA23)),size=1) +  geom_point(aes(color = as.factor(SA23)),size=3)+ylab("uM") + xlab("Time (h)") +ggtitle("Apical side") + theme(legend.position="bottom")+  theme(legend.key.size = unit(0.2, "cm")) +theme(legend.title=element_blank())+
          theme_bw()+ scale_y_continuous( trans="log10")+labs(color='Surface Area (cm2)')+ theme(legend.position=c(0.8, 0.8)) + theme(legend.text = element_text(size=15)) + theme(legend.key.size= unit(1, 'cm')) +theme(axis.text.x=element_text(size=15)) + theme(axis.text.y=element_text(size=15))) 

p2=plot(ggplot(sim_samples,aes(time/60,C3ObsS)) + geom_line(aes(color = as.factor(SA23)),size=1) +  geom_point(aes(color = as.factor(SA23)),size=3)+ylab("uM") + xlab("Time (h)") +ggtitle("Basolateral  side")+ theme(legend.position="bottom")+  theme(legend.key.size = unit(0.2, "cm")) +theme(legend.title=element_blank())+
          theme_bw()+  scale_y_continuous( trans="log10") + geom_hline(yintercept=0.001, linetype="dashed", color = "black")+labs(color='Surface Area (cm2)')+theme(legend.position = "none")+theme(axis.text.x=element_text(size=15))+ theme(axis.text.y=element_text(size=15))) 

p3=plot(ggplot(sim_samples,aes(time/60,C4ObsS)) + geom_line(aes(color = as.factor(SA23)),size=1) +  geom_point(aes(color = as.factor(SA23)),size=3)+ylab("uM") + xlab("Time (h)") +ggtitle("Liver  MPS")+ theme(legend.position="bottom")+  theme(legend.key.size = unit(0.2, "cm")) +theme(legend.title=element_blank())+
          theme_bw()+  scale_y_continuous( trans="log10")+ geom_hline(yintercept=0.001, linetype="dashed", color = "black")+labs(color='Surface Area (cm2)')+theme(legend.position = "none")+theme(axis.text.x=element_text(size=15))+ theme(axis.text.y=element_text(size=15))) 

grid.arrange(p1,p2,p3,nrow = 1)

p2zoom=plot(ggplot(sim_samples,aes(time/60,C3ObsS)) + geom_line(aes(color = as.factor(SA23)), size=1) +  geom_point(aes(color = as.factor(SA23)),size=3)+ylab("uM") + xlab("Time (h)") +ggtitle("Basolateral  side (0-0.01 uM)")+ theme(legend.position="bottom")+  theme(legend.key.size = unit(0.2, "cm")) +theme(legend.title=element_blank())+ylim(0,0.01)+
          theme_bw()+ geom_hline(yintercept=0.001, linetype="dashed", color = "black")+geom_text(aes(0,0.001,label = "BLQ", vjust = -1),color="black")+labs(color='Surface Area (cm2)')+theme(legend.position = "none")) 

grid.arrange(p1,p2,p3,nrow = 1)


##### Same model but the investigation here considers the ratio SA/V


ratio_i= c(0.33,0.5,1,2,3)

sim = data.frame()
for (i in 1:length(ratio_i)){
  
  #Parameters
  param=c( V2=500*ratio_i[i]  ,V3=1506 , V4 =1394, fum =  1,fum2 = 1, fumM =  1,fum2M = 1, NPERM=2,Papp=0.4,CL3tot=0, CL4tot=0,CL56S=10000000, SA23=1.12,N3=0.5,N4=0.5,convfactor=0.006, Med_Q=60  )
  
  #Event table                           
  qd <- et(seq(0, 2880, by = 1))
  
  #Ini. cond                          
  # Initialization 
  initialization <- c(A2.S = ini.Conc* as.numeric(param["V2"]), A3.S = 0)
  
  #
  omega.matrix= lotri({c(eta.Papp ~  0.00001) })
  
  sigma.matrix <- lotri(CEps.C2.S ~0.000015^2, CEps.C3.S ~ 0.00000015^2, CEps.C4.S ~ 0.000015^2)
  ### Sensitivity analysis ####
  set.seed(120)
  sim_i <- rxSolve(object=gut.liver.pop, params =param,events =qd,inits=initialization, omega=omega.matrix, sigma=sigma.matrix, sub=1)
  
  sim_i = as.data.frame(sim_i)
  
  sim_i$ratio_i= ratio_i[i]
  
  sim= rbind(sim, sim_i)
}



sim_samples = subset(sim, time %in% c(0,30,60,120,240,480,1440,2880,4320,5999))


p11=plot(ggplot(sim_samples,aes(time/60,C2ObsS)) +  geom_line(aes(color = as.factor(ratio_i)),size=1) +  geom_point(aes(color = as.factor(ratio_i)),size=3)+ylab("uM") + xlab("Time (h)") +ggtitle("Apical side") + theme(legend.position="bottom")+  theme(legend.key.size = unit(0.2, "cm")) +theme(legend.title=element_blank())+
          theme_bw()+ scale_y_continuous( trans="log10")+labs(color='Volume Apical Side x-fold')+ theme(legend.position=c(0.2, 0.3)) + theme(legend.text = element_text(size=15)) + theme(legend.key.size= unit(1, 'cm')) +theme(axis.text.x=element_text(size=15))+ theme(axis.text.y=element_text(size=15)))

p22=plot(ggplot(sim_samples,aes(time/60,C3ObsS)) + geom_line(aes(color = as.factor(ratio_i)),size=1) +  geom_point(aes(color = as.factor(ratio_i)),size=3)+ylab("uM") + xlab("Time (h)") +ggtitle("Basolateral  side")+ theme(legend.position="bottom")+  theme(legend.key.size = unit(0.2, "cm")) +theme(legend.title=element_blank())+
          theme_bw()+  scale_y_continuous( trans="log10") + geom_hline(yintercept=0.001, linetype="dashed", color = "black")+labs(color='Surface Area (cm2)')+theme(legend.position = "none")+theme(axis.text.x=element_text(size=15))+ theme(axis.text.y=element_text(size=15))) 

p33=plot(ggplot(sim_samples,aes(time/60,C4ObsS)) + geom_line(aes(color = as.factor(ratio_i)),size=1) +  geom_point(aes(color = as.factor(ratio_i)),size=3)+ylab("uM") + xlab("Time (h)") +ggtitle("Liver  MPS")+ theme(legend.position="bottom")+  theme(legend.key.size = unit(0.2, "cm")) +theme(legend.title=element_blank())+
          theme_bw()+  scale_y_continuous( trans="log10")+ geom_hline(yintercept=0.001, linetype="dashed", color = "black")+labs(color='Surface Area (cm2)')+theme(legend.position = "none")+theme(axis.text.x=element_text(size=15))+ theme(axis.text.y=element_text(size=15))) 

grid.arrange(p11,p22,p33,nrow = 1)


p11=plot(ggplot(sim_samples,aes(time/60,A2.S)) +  geom_line(aes(color = as.factor(ratio_i)),size=1) +  geom_point(aes(color = as.factor(ratio_i)),size=3)+ylab("pmol") + xlab("Time (h)") +ggtitle("Apical side") +theme(legend.position = "none")+
           theme_bw()+ scale_y_continuous( trans="log10")  +theme(axis.text.x=element_text(size=15))+ theme(axis.text.y=element_text(size=15))) 

p22=plot(ggplot(sim_samples,aes(time/60,A3.S)) + geom_line(aes(color = as.factor(ratio_i)),size=1) +  geom_point(aes(color = as.factor(ratio_i)),size=3)+ylab("pmol") + xlab("Time (h)") +ggtitle("Basolateral  side")+ theme(legend.position="bottom")+  theme(legend.key.size = unit(0.2, "cm")) +theme(legend.title=element_blank())+ 
           theme_bw()+  scale_y_continuous( trans="log10") + geom_hline(yintercept=0.001, linetype="dashed", color = "black")+labs(color='Surface Area (cm2)')+theme(legend.position = "none")+theme(axis.text.x=element_text(size=15))+ theme(axis.text.y=element_text(size=15))) 

p33=plot(ggplot(sim_samples,aes(time/60,A4.S)) + geom_line(aes(color = as.factor(ratio_i)),size=1) +  geom_point(aes(color = as.factor(ratio_i)),size=3)+ylab("pmol") + xlab("Time (h)") +ggtitle("Liver  MPS")+ theme(legend.position="bottom")+  theme(legend.key.size = unit(0.2, "cm")) +theme(legend.title=element_blank())+
           theme_bw()+  scale_y_continuous( trans="log10")+ geom_hline(yintercept=0.001, linetype="dashed", color = "black")+labs(color='Surface Area (cm2)')+theme(legend.position = "none")+theme(axis.text.x=element_text(size=15))+ theme(axis.text.y=element_text(size=15))) 


grid.arrange(p11,p22,p33,nrow = 1)


##################                                              ####################
##################  Flow rate between basol side and liver MPS  ####################
##################                                              ####################
ini.Conc=10

gut.liver.pop <-  RxODE({         # General model
  # Volume of the cells
  V5 = 1.7    # Typical enterocyte volume of 1 Mio of cells
  V6 = 3.9    # Typical enterocyte volume of 1 Mio of cells
  
  Q =Med_Q
  ### Concentration/amount relationship for the substrate
  C2.S = A2.S / V2               # Substrate concentration in the AP
  C3.S = A3.S / V3               # Substrate concentration in the BS
  C4.S = A4.S  / V4               # Substrate concentration in the LMPS
  C5.S = A5.S  / (N3 * V5)       # Substrate concentration in the ent.
  C6.S = A6.S  / (N4 * V6)       # Substrate concentration in the hep.
  
  
  CL23 = SA23 * convfactor*Papp*exp(eta.Papp)
  
  
  
  # Substrate
  d/dt(A2.S) = - CL23 * NPERM * (C2.S * fum2 - C5.S) # AP. passive diffusion ent.
  d/dt(A3.S) = CL23 * NPERM * (C5.S- C3.S * fum)     # BS. passive diffusion
  + Q  * C4.S 
  - Q  * C3.S                 # Circulation flow out (BS -> LIVER MPS)     
  d/dt(A4.S) =  Q * C3.S                      # Circulation flow in (MC -> LMPS)
  - Q * C4.S                       # Circulation flow out (LMPS -> MC)        
  - CL56S * (C4.S * fum - C6.S)            # Passive diffusion hep.
  d/dt(A5.S) = CL23 * NPERM * (C2.S * fum2 - C5.S)   # ENT. passive diffusion .ent apical side
  - CL23 * NPERM *(C5.S - C3.S * fum)     # passive diffusion .ent basolateral side
  - CL3tot * N3 * C5.S                     # Intestinal metabolism
  d/dt(A6.S) = CL56S * (C4.S * fum - C6.S)            # HEP. passive diffusion hep.
  - CL4tot * N4 * C6.S                     # Metabolism in the liver cells
  
  
  C2ObsS = C2.S *exp(CEps.C2.S) 
  C3ObsS = C3.S *exp(CEps.C3.S) 
  C4ObsS = C3.S *exp(CEps.C4.S) 
  
})

Q_i= c(10,30,50)
#Q_i= c(1000,5000,9000)

sim = data.frame()
for (i in 1:length(Q_i)){
  
  #Parameters
  param=c( V2=325  ,V3=1506 , V4 =1394, fum =  1,fum2 = 1, fumM =  1,fum2M = 1, NPERM=2,Papp=800,CL3tot=0, CL4tot=120,CL56S=10000000, SA23=0.33,N3=0.5,N4=0.5,convfactor=0.006, Med_Q=Q_i[i]  )
  
  #Event table                           
  qd <- et(seq(0, 2880, by = 1))
  
  #Ini. cond                          
  # Initialization 
  initialization <- c(A2.S = ini.Conc* as.numeric(param["V2"]), A3.S = 0)
  
  #
  omega.matrix= lotri({c(eta.Papp ~  0.00001) })
  
  sigma.matrix <- lotri(CEps.C2.S ~0.15^2, CEps.C3.S ~ 0.15^2, CEps.C4.S ~ 0.15^2)
  ### Sensitivity analysis ####
  set.seed(120)
  sim_i <- rxSolve(object=gut.liver.pop, params =param,events =qd,inits=initialization, omega=omega.matrix, sigma=sigma.matrix, sub=1)
  
  sim_i = as.data.frame(sim_i)
  

  sim= rbind(sim, sim_i)
}
sim1=sim
sim1$Qbl=30

#sim_samples = subset(sim, time %in% c(0,30,60,120,240,1440,2880,4320))
sim_samples=sim

p1=plot(ggplot(sim_samples,aes(time/60,C2.S)) +  geom_line(aes(color = as.factor(Q))) +ylab("uM") + xlab("Time (h)") +ggtitle("Apical side") + theme(legend.position="bottom")+  theme(legend.key.size = unit(0.2, "cm")) +theme(legend.title=element_blank())+
          theme_bw()+ scale_y_continuous(limits = c(0, 10.52))+labs(color='Q (uL/min)')+ theme(legend.position=c(0.75, 0.3))) 

p2=plot(ggplot(sim_samples,aes(time/60,C3.S)) + geom_line(aes(color = as.factor(Q))) +ylab("uM") + xlab("Time (h)") +ggtitle("Basolateral  side")+ theme(legend.position="bottom")+  theme(legend.key.size = unit(0.2, "cm")) +theme(legend.title=element_blank())+
          theme_bw()+labs(color='Surface Area (cm2)')+theme(legend.position = "none")) 

p3=plot(ggplot(sim_samples,aes(time/60,C4.S)) + geom_line(aes(color = as.factor(Q))) +ylab("uM") + xlab("Time (h)") +ggtitle("Liver  MPS")+ theme(legend.position="bottom")+  theme(legend.key.size = unit(0.2, "cm")) +theme(legend.title=element_blank())+
          theme_bw()+labs(color='Surface Area (cm2)')+theme(legend.position = "none")) 

grid.arrange(p1,p2,p3,nrow = 1)


######### Second set
Q_i=c(96,120,144)
#Q_i= c(1000,5000,9000)

sim = data.frame()
for (i in 1:length(Q_i)){
  
  #Parameters
   param=c( V2=325  ,V3=1594 , V4 =1394, fum =  1,fum2 = 1, fumM =  1,fum2M = 1, NPERM=2,Papp=500,CL3tot=0, CL4tot=120,CL56S=10000000, SA23=0.33,N3=0.5,N4=0.5,convfactor=0.006, Med_Q=Q_i[i]  )
  
  #Event table                           
  qd <- et(seq(0, 2880, by = 1))
  
  #Ini. cond                          
  # Initialization 
  initialization <- c(A2.S = ini.Conc* as.numeric(param["V2"]), A3.S = 0)
  
  #
  omega.matrix= lotri({c(eta.Papp ~  0.00001) })
  
  sigma.matrix <- lotri(CEps.C2.S ~0.15^2, CEps.C3.S ~ 0.15^2, CEps.C4.S ~ 0.15^2)
  ### Sensitivity analysis ####
  set.seed(120)
  sim_i <- rxSolve(object=gut.liver.pop, params =param,events =qd,inits=initialization, omega=omega.matrix, sigma=sigma.matrix, sub=1)
  
  sim_i = as.data.frame(sim_i)
  
  
  sim= rbind(sim, sim_i)
}



#sim_samples = subset(sim, time %in% c(0,30,60,120,240,1440,2880,4320))
#sim_samples = subset(sim, time %in% c(0,30,60,120,240,1440,2880,4320))
sim_samples=sim

p1=plot(ggplot(sim_samples,aes(time/60,C2.S)) +  geom_line(aes(color = as.factor(Q))) +ylab("uM") + xlab("Time (h)") +ggtitle("Apical side") + theme(legend.position="bottom")+  theme(legend.key.size = unit(0.2, "cm")) +theme(legend.title=element_blank())+
          theme_bw()+ scale_y_continuous(limits = c(0, 10.52))+labs(color='Q (uL/min)')+ theme(legend.position=c(0.75, 0.3))) 

p2=plot(ggplot(sim_samples,aes(time/60,C3.S)) + geom_line(aes(color = as.factor(Q))) +ylab("uM") + xlab("Time (h)") +ggtitle("Basolateral  side")+ theme(legend.position="bottom")+  theme(legend.key.size = unit(0.2, "cm")) +theme(legend.title=element_blank())+
          theme_bw()+labs(color='Surface Area (cm2)')+theme(legend.position = "none")) 

p3=plot(ggplot(sim_samples,aes(time/60,C4.S)) + geom_line(aes(color = as.factor(Q))) +ylab("uM") + xlab("Time (h)") +ggtitle("Liver  MPS")+ theme(legend.position="bottom")+  theme(legend.key.size = unit(0.2, "cm")) +theme(legend.title=element_blank())+
          theme_bw()+labs(color='Surface Area (cm2)')+theme(legend.position = "none")) 

grid.arrange(p1,p2,p3,nrow = 1)

sim2=sim
sim2$Qbl = 120

sim_all=rbind(sim1,sim2)

p1=plot(ggplot(sim_all,aes(time/60,C2.S)) +  geom_line(aes(color = as.factor(Qbl))) +ylab("uM") + xlab("Time (h)") +ggtitle("Apical side") + theme(legend.position="bottom")+  theme(legend.key.size = unit(0.2, "cm")) +theme(legend.title=element_blank())+
          theme_bw()+ scale_y_continuous(limits = c(0, 10.52))+labs(color='Q (uL/min)')+ theme(legend.position=c(0.75, 0.3))) 

p2=plot(ggplot(sim_all,aes(time/60,C3.S)) + geom_line(aes(color = as.factor(Qbl))) +ylab("uM") + xlab("Time (h)") +ggtitle("Basolateral  side")+ theme(legend.position="bottom")+  theme(legend.key.size = unit(0.2, "cm")) +theme(legend.title=element_blank())+
          theme_bw()+labs(color='Surface Area (cm2)')+theme(legend.position = "none")) 

p3=plot(ggplot(sim_all,aes(time/60,C4.S)) + geom_line(aes(color = as.factor(Qbl))) +ylab("uM") + xlab("Time (h)") +ggtitle("Liver  MPS")+ theme(legend.position="bottom")+  theme(legend.key.size = unit(0.2, "cm")) +theme(legend.title=element_blank())+
          theme_bw()+labs(color='Surface Area (cm2)')+theme(legend.position = "none")) 

grid.arrange(p1,p2,p3,nrow = 1)




##################                ####################
##################  Cell density  ####################
##################                ####################



gut.liver.MM <-  RxODE({         
  N3 = tN3 * exp(eta.N3)
  N4 = tN4 * exp(eta.N4)
  Qr = tQ * exp(eta.Q)
  kev = tkev * exp(eta.kev)
  
  d/dt(Vmsa) = 0
  #d/dt(Vmsb) = 0
  d/dt(Vmsl) = 0
  
  V2s = V2 - Vmsa
  V3evs =  V3-kev*time*(Ab/(Ab+Al)) 
  Ve=kev*time
  V4evs = V4-kev*time*(Al/(Ab+Al)) - Vmsl #- Vmsb
  
  #ter=eta.Q
  
  V2s = V2s -Vtab
  V4evs = V4evs+Vtab  
  d/dt(Vtab) =   ktab * V2s
  
  #C2MM = A2MM/V2s
  C2MPA = A2MPA/V2s
  C3MPA = A3MPA/V3evs
  C4MPA = A4MPA/V4evs
  C5MPA = A5MPA/(V5*N3)
  C2GMPA = A2GMPA/V2s
  C3GMPA = A3GMPA/V3evs
  C4GMPA = A4GMPA/V4evs
  C5GMPA = A5GMPA/(V5*N3)
  
  #  d/dt(timein) = 0    # Time delay after activation of interconnection flow
  #  Q = tQ * timein
  
  #d/dt(A2MM) =  - CL3MM * N3 * C2MM *fum
 # f(A2MM)    = (V2s - Sample.Va) / V2s    
  d/dt(A2MPA) =  PappMPA* SA23 * convfactor * NPERM * (Er*C5MPA- C2MPA * fum)  
  f(A2MPA)    = (V2s - Sample.Va_i) / V2s 
  d/dt(A3MPA) =  PappMPA* SA23 * convfactor * NPERM * (C5MPA- C3MPA * fumMPA)- Qr* C3MPA +  Qr  * C4MPA  
  #  f(A3MPA)    = (V3evs - Sample.V) / V3evs 
  d/dt(A4MPA) =  Qr  * C3MPA - Qr * C4MPA - C4MPA * N4 * fumMPA * CL4MPA
  f(A4MPA)    = (V4evs - Sample.V) / V4evs 
  d/dt(A5MPA) =  - PappMPA* SA23 * convfactor * NPERM *(C5MPA - C3MPA * fumMPA) - PappMPA* SA23 * convfactor * NPERM *(Er*C5MPA - C2MPA * fum) - CL3MPA * N3 * C5MPA
  
  d/dt(A2GMPA) =  PappGMPA* SA23 * convfactor * NPERM * (ErG*C5GMPA- C2GMPA * fum)  
  f(A2GMPA)    = (V2s - Sample.Va_i) / V2s 
  d/dt(A3GMPA) =  PappGMPA* SA23 * convfactor * NPERM * (C5GMPA- C3GMPA * fumGMPA)- Qr* C3GMPA +  Qr  * C4GMPA 
  #  f(A3GMPA)    = (V3evs - Sample.V) / V3evs 
  d/dt(A4GMPA) =  Qr  * C3GMPA - Qr * C4GMPA  +  C4MPA * N4 * fumMPA * CL4MPA 
  f(A4GMPA)    = (V4evs - Sample.V) / V4evs 
  d/dt(A5GMPA) =  - PappGMPA* SA23 * convfactor * NPERM *(C5GMPA - C3GMPA * fumGMPA) - PappGMPA* SA23 * convfactor * NPERM *(ErG*C5GMPA - C2GMPA * fum) + CL3MPA * N3 * C5MPA 
  
  # cell desnity in the apical side (enterocytes)
  cell_den = N3/V2s
  
#  C2ObsMM=C2MM *exp(CEps.C2.MM)
  C2ObsMPA=C2MPA *exp(CEps.C2.MPA)
  C3ObsMPA=C3MPA *exp(CEps.C3.MPA)
  C4ObsMPA=C4MPA *exp(CEps.C4.MPA)
  
  C2ObsGMPA=C2GMPA *exp(CEps.C2.GMPA)
  C3ObsGMPA=C3GMPA *exp(CEps.C3.GMPA)
  C4ObsGMPA=C4GMPA *exp(CEps.C4.GMPA)
  
})

# Inizialitzation flow 
timeQ = 1 # Start flow after x min
time.inc = 2880 # incubation time in min
Conc.ini=10 # uM Mofetil
o.target.AP <- c(30,60,120,480,2880)            #114
o.target.BS <- c(30,60,120,480,1440,2880)        #293
o.target.LM <- c(120,240,480,1440,2880)
Sample.V=25    # sampling volume uL bas and liver
Sample.Va=c(5,15,30,45)   # sampling volume uL apical

t12tab = 1000000000  # 2880   # t1/2 (first oreder process) in min



sim.rep.pop.all = data.frame()
for (i in 1:length(Sample.Va)){
  Sample.Va_i = Sample.Va[i]
  
  # paramters
  par<- c(tQ=90,                 
          V2=325, 
          V3=1506,             
          V4=1394,
          PappMPA=546,          # Evaluated from  exp (nm/s)
          SA23=0.33, 
          convfactor=0.006,
          NPERM=2,
          fum=1,
          fumMPA=0.38,
          fumGMPA=0.70,
          CL3MPA=17,         # Evaluated from  exp (uL/min/Mio cells)
          CL3MM = 13,         # Parameter from EXP5 all simul. fitting
          V5=2.6,
          tN3=0.45,           # From mean of gut-liver data Exp 5
          tN4=0.30,          # From mean of gut-liver data Exp 5
          Er=3.0,                # From Cn Bio model MG with input
          ErG=3.1,
          PappGMPA=0.35,
          CL4MPA= 26,
          Sample.V = Sample.V,
          Sample.Va_i = Sample.Va[i],
          tkev = 0,   # From exp 5 (mean of values)
          Al = 7,
          Ab = 0,
          ktab = log(2)/t12tab)
          
  sim.ev.MM <- et(seq(0,tail(o.target.AP,1), by = 1)) %>%
    et(cmt="C2ObsMM", time=0, evid = 1, amt=Conc.ini*unname(par["V2"])) %>%      
    et(cmt="C2ObsMPA", time=0, evid = 1, amt=0)  %>%    
    et(cmt="C2ObsGMPA", time=0, evid = 1, amt=0) %>%   
    et(cmt="C3ObsMPA", time=0, evid = 1, amt=0)  %>%    
    et(cmt="C3ObsGMPA", time=0, evid = 1, amt=0) %>%   
    et(cmt="C4ObsMPA", time=0, evid = 1, amt=0)  %>%    
    et(cmt="C4ObsGMPA", time=0, evid = 1, amt=0) %>%   
    
    # Sampling volume in the Apical 
    et(cmt="Vmsa", time=o.target.AP[1], evid = 1, amt=Sample.Va_i) %>%      
    et(cmt="Vmsa", time=o.target.AP[2], evid = 1, amt=Sample.Va_i) %>% 
    et(cmt="Vmsa", time=o.target.AP[3], evid = 1, amt=Sample.Va_i) %>% 
    et(cmt="Vmsa", time=o.target.AP[4], evid = 1, amt=Sample.Va_i) %>% 
    et(cmt="Vmsa", time=o.target.AP[5], evid = 1, amt=Sample.Va_i) %>% 
    
    # Sampling amount of substrate in Apical
    et(cmt="A2MM", time=o.target.AP[1], evid=6, amt=1) %>%
    et(cmt="A2MM", time=o.target.AP[2], evid=6, amt=1) %>%
    et(cmt="A2MM", time=o.target.AP[3], evid=6, amt=1) %>%
    et(cmt="A2MM", time=o.target.AP[4], evid=6, amt=1) %>%
    et(cmt="A2MM", time=o.target.AP[5], evid=6, amt=1) %>%
    
    # Sampling amount of MPA in Apical
    et(cmt="A2MPA", time=o.target.AP[1], evid=6, amt=1) %>%
    et(cmt="A2MPA", time=o.target.AP[2], evid=6, amt=1) %>%
    et(cmt="A2MPA", time=o.target.AP[3], evid=6, amt=1) %>%
    et(cmt="A2MPA", time=o.target.AP[4], evid=6, amt=1) %>%
    et(cmt="A2MPA", time=o.target.AP[5], evid=6, amt=1) %>%
    
    # Sampling amount of MPA in Apical
    et(cmt="A2GMPA", time=o.target.AP[1], evid=6, amt=1) %>%
    et(cmt="A2GMPA", time=o.target.AP[2], evid=6, amt=1) %>%
    et(cmt="A2GMPA", time=o.target.AP[3], evid=6, amt=1) %>%
    et(cmt="A2GMPA", time=o.target.AP[4], evid=6, amt=1) %>%
    et(cmt="A2GMPA", time=o.target.AP[5], evid=6, amt=1) %>%
    
    # # Sampling volume in the basolateral 
    # et(cmt="Vmsb", time=o.target.BS[1], evid = 1, amt=Sample.V) %>%      
    # et(cmt="Vmsb", time=o.target.BS[2], evid = 1, amt=Sample.V) %>% 
    # et(cmt="Vmsb", time=o.target.BS[3], evid = 1, amt=Sample.V) %>% 
    # et(cmt="Vmsb", time=o.target.BS[4], evid = 1, amt=Sample.V) %>% 
    # et(cmt="Vmsb", time=o.target.BS[5], evid = 1, amt=Sample.V) %>% 
    # et(cmt="Vmsb", time=o.target.BS[6], evid = 1, amt=Sample.V) %>% 
    
    # # Sampling amount of MPA in Basolateral 
    # et(cmt="A3MPA", time=o.target.BS[1], evid=6, amt=1) %>%
  # et(cmt="A3MPA", time=o.target.BS[2], evid=6, amt=1) %>%
  # et(cmt="A3MPA", time=o.target.BS[3], evid=6, amt=1) %>%
  # et(cmt="A3MPA", time=o.target.BS[4], evid=6, amt=1) %>%
  # et(cmt="A3MPA", time=o.target.BS[5], evid=6, amt=1) %>%
  # et(cmt="A3MPA", time=o.target.BS[6], evid=6, amt=1) %>%
  # 
  # # Sampling amount of GMPA in basolateral
  # et(cmt="A3GMPA", time=o.target.BS[1], evid=6, amt=1) %>%
  # et(cmt="A3GMPA", time=o.target.BS[2], evid=6, amt=1) %>%
  # et(cmt="A3GMPA", time=o.target.BS[3], evid=6, amt=1) %>%
  # et(cmt="A3GMPA", time=o.target.BS[4], evid=6, amt=1) %>%
  # et(cmt="A3GMPA", time=o.target.BS[5], evid=6, amt=1) %>%
  # et(cmt="A3GMPA", time=o.target.BS[6], evid=6, amt=1) %>%
  # 
  # Sampling volume in the liver 
  et(cmt="Vmsl", time=o.target.LM[1], evid = 1, amt=Sample.V) %>%      
    et(cmt="Vmsl", time=o.target.LM[2], evid = 1, amt=Sample.V) %>% 
    et(cmt="Vmsl", time=o.target.LM[3], evid = 1, amt=Sample.V) %>% 
    et(cmt="Vmsl", time=o.target.LM[4], evid = 1, amt=Sample.V) %>% 
    et(cmt="Vmsl", time=o.target.LM[5], evid = 1, amt=Sample.V) %>% 
    et(cmt="Vmsl", time=o.target.BS[1], evid = 1, amt=Sample.V) %>%
    et(cmt="Vmsl", time=o.target.BS[2], evid = 1, amt=Sample.V) %>%
    et(cmt="Vmsl", time=o.target.BS[3], evid = 1, amt=Sample.V) %>%
    et(cmt="Vmsl", time=o.target.BS[4], evid = 1, amt=Sample.V) %>%
    et(cmt="Vmsl", time=o.target.BS[5], evid = 1, amt=Sample.V) %>%
    et(cmt="Vmsl", time=o.target.BS[6], evid = 1, amt=Sample.V) %>%
    
    # Sampling amount of MPA in Liver 
    et(cmt="A4MPA", time=o.target.LM[1], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.LM[2], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.LM[3], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.LM[4], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.LM[5], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.BS[1], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.BS[2], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.BS[3], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.BS[4], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.BS[5], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.BS[6], evid=6, amt=1) %>%
    
    # Sampling amount of GMPA in Liver
    et(cmt="A4GMPA", time=o.target.LM[1], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.LM[2], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.LM[3], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.LM[4], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.LM[5], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.BS[1], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.BS[2], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.BS[3], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.BS[4], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.BS[5], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.BS[6], evid=6, amt=1) 
  
  # Omega and sigma matrix 
  
  omega.matrix.mod1.pop <- lotri({
    c(eta.N4 ~ 0.000015^2, eta.N3 ~ 0.000015^2, eta.Q ~ 0.00002^2, eta.kev ~ 0.000015^2) })
  
  
  sigma.matrix.mod1.pop <- lotri( CEps.C2.MPA ~ 0.15^2, CEps.C3.MPA ~ 0.15^2, CEps.C4.MPA ~ 0.15^2, 
                                  CEps.C2.GMPA ~ 0.15^2, CEps.C3.GMPA ~ 0.15^2, CEps.C4.GMPA ~ 0.15^2)
  
  initialization <- c(A2MPA=Conc.ini*unname(par["V2"]), A3MPA = 0, A4MPA = 0, A5MPA =0, A3GMPA = 0, A4GMPA = 0, A5GMPA =0, Vtab=0)
  
  #qd <- et(seq(0,time.inc, by = 1)) %>%
  #   et(cmt="A2MM", time=0, evid = 1, amt=Conc.ini*325) 
  
  qd <- data.frame(qd)
  
  
  set.seed(123)
  sim.rep.pop <- rxSolve(gut.liver.MM,par,sim.ev.MM,inits= initialization, omega=omega.matrix.mod1.pop, addDosing=TRUE,sigma=sigma.matrix.mod1.pop, nSub=2)
  
  sim.rep.pop<-data.frame(sim.rep.pop)
  sim.rep.pop.plot.2 <- sim.rep.pop %>% filter(sim.rep.pop$evid == 0)
  sim.rep.pop.plot.2 <- sim.rep.pop.plot.2 %>% filter(sim.rep.pop.plot.2$time >1)
  
 # sim.rep.pop.plot.2 <- sim.rep.pop.plot.2 %>% filter(sim.rep.pop.plot.2$sim.id < 1001 & !sim.rep.pop.plot.2$sim.id == 0)
  
  # define sampling volume columns
  sim.rep.pop.plot.2$VsamplingA = Sample.Va[i]
  sim.rep.pop.all= rbind(sim.rep.pop.all, sim.rep.pop.plot.2)
}

x_ax = theme(axis.text.x=element_text(size=15),axis.title=element_text(size=15))
y_ax = theme(axis.text.y=element_text(size=15),axis.title=element_text(size=15))


PV2<-ggplot(sim.rep.pop.all,aes(time/60,V2s)) + geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Volume Apical Side (uL)") + xlab("Time (h)") + theme_bw()  +xlim(0,25) + ylim(80,325) +theme(legend.position = "none") +geom_hline(yintercept=325-325*0.15, linetype="dashed", color = "red") +x_ax+y_ax
Pcell_d<-ggplot(sim.rep.pop.all,aes(time/60,cell_den*1000)) +  geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Cells density (Mio cells/mL)") + xlab("Time (h)")+xlim(0,25)  + theme_bw()  +  theme(legend.position = "none") +  geom_hline(yintercept=(sim.rep.pop.all[1,"cell_den"]+sim.rep.pop.all[1,"cell_den"]*0.15)*1000, linetype="dashed", color = "red")+x_ax+y_ax
C2<-ggplot(sim.rep.pop.all,aes(time/60,C2MPA)) +  geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Conc. Sub. Apical Side (uM)") + xlab("Time (h)") + theme_bw()  +xlim(0,25)    + labs(color='Sampling Vol. (uL)') + scale_y_continuous(trans='log10') +  theme(legend.position=c(0.25, 0.25))+x_ax+y_ax +theme(legend.text=element_text(size=rel(1.3)))+theme(legend.title =element_text(size=rel(1.3)))# + theme(legend.position = "none")
A2<-ggplot(sim.rep.pop.all,aes(time/60,A2MPA)) +  geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Amount Sub. Apical Side (pmol)") + xlab("Time (h)") + theme_bw() +xlim(0,25) + theme(legend.position = "none") + scale_y_continuous(trans='log10')+x_ax+y_ax

A2_zoom<-ggplot(sim.rep.pop.all,aes(time/60,A2MPA)) +  geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Amount Apical Side (pmol)") + xlab("Time (h)") + theme_bw() + xlim(0,6)

grid.arrange(PV2,Pcell_d,C2,A2)


C3<-ggplot(sim.rep.pop.all,aes(time/60,C3MPA)) + geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Conc. Sub. Basol. Side (uM)") + xlab("Time (h)") + theme_bw()   +theme(legend.position = "none") 
C4<-ggplot(sim.rep.pop.all,aes(time/60,C4MPA)) +  geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Conc. Sub. Liver MPS (uM)") + xlab("Time (h)") + theme_bw()  +  theme(legend.position = "none") 
#C2<-ggplot(sim.rep.pop.all,aes(time/60,C2MPA)) +  geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Conc. Apical Side (uM)") + xlab("Time (h)") + theme_bw()     +  theme(legend.position=c(0.75, 0.55))+ labs(color='Media Evap. (uL/min)')
#A2<-ggplot(sim.rep.pop.all,aes(time/60,A2MPA)) +  geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Amount Apical Side (pmol)") + xlab("Time (h)") + theme_bw() + theme(legend.position = "none")
C2M<-ggplot(sim.rep.pop.all,aes(time/60, C2GMPA)) + geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Conc. Met. Apical Side (uM)") + xlab("Time (h)") + theme_bw()  +theme(legend.position = "none") 
C3M<-ggplot(sim.rep.pop.all,aes(time/60, C3GMPA)) +  geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Conc. Met Basol. Side (uM)") + xlab("Time (h)") + theme_bw()  +  theme(legend.position = "none") 
C4M<-ggplot(sim.rep.pop.all,aes(time/60, C4GMPA)) +  geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Conc. Met. Liver MPS (uM)") + xlab("Time (h)") + theme_bw()  +  theme(legend.position = "none") 


#A2_zoom<-ggplot(sim.rep.pop.all,aes(time/60,A2MPA)) +  geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Amount Apical Side (pmol)") + xlab("Time (h)") + theme_bw() + xlim(0,6)

grid.arrange(PV2,Pcell_d,C2,A2)

grid.arrange(C2,C3,C4,C2M,C3M,C4M,ncol = 3)






gut.liver.MM <-  RxODE({         
  N3 = tN3 * exp(eta.N3)
  N4 = tN4 * exp(eta.N4)
  Qr = tQ * exp(eta.Q)
  kev = tkev * exp(eta.kev)
  
  d/dt(Vmsa) = 0
  #d/dt(Vmsb) = 0
  d/dt(Vmsl) = 0 
  
  V2s = V2 -kev.a*time- Vmsa
  V3evs =  V3-kev*time*(Ab/(Ab+Al)) 
  Ve=kev*time
  V4evs = V4-kev*time*(Al/(Ab+Al)) - Vmsl #- Vmsb
  
  #ter=eta.Q
  
  V2s = V2s -Vtab
  V4evs = V4evs+Vtab  
  d/dt(Vtab) =   ktab * V2s
 
  #C2MM = A2MM/V2s
  C2MPA = A2MPA/V2s
  C3MPA = A3MPA/V3evs
  C4MPA = A4MPA/V4evs
  C5MPA = A5MPA/(V5*N3)
  C2GMPA = A2GMPA/V2s
  C3GMPA = A3GMPA/V3evs
  C4GMPA = A4GMPA/V4evs
  C5GMPA = A5GMPA/(V5*N3)
  
  #  d/dt(timein) = 0    # Time delay after activation of interconnection flow
  #  Q = tQ * timein
  
  #d/dt(A2MM) =  - CL3MM * N3 * C2MM *fum
  # f(A2MM)    = (V2s - Sample.Va) / V2s    
  d/dt(A2MPA) =  PappMPA* SA23 * convfactor * NPERM * (Er*C5MPA- C2MPA * fum)  
  f(A2MPA)    = (V2s - Sample.Va) / V2s 
  d/dt(A3MPA) =  PappMPA* SA23 * convfactor * NPERM * (C5MPA- C3MPA * fumMPA)- Qr* C3MPA +  Qr  * C4MPA  
  #  f(A3MPA)    = (V3evs - Sample.V) / V3evs 
  d/dt(A4MPA) =  Qr  * C3MPA - Qr * C4MPA - C4MPA * N4 * fumMPA * CL4MPA
  f(A4MPA)    = (V4evs - Sample.V) / V4evs 
  d/dt(A5MPA) =  - PappMPA* SA23 * convfactor * NPERM *(C5MPA - C3MPA * fumMPA) - PappMPA* SA23 * convfactor * NPERM *(Er*C5MPA - C2MPA * fum) - CL3MPA * N3 * C5MPA
  
  d/dt(A2GMPA) =  PappGMPA* SA23 * convfactor * NPERM * (ErG*C5GMPA- C2GMPA * fum)  
  f(A2GMPA)    = (V2s - Sample.Va) / V2s 
  d/dt(A3GMPA) =  PappGMPA* SA23 * convfactor * NPERM * (C5GMPA- C3GMPA * fumGMPA)- Qr* C3GMPA +  Qr  * C4GMPA 
  #  f(A3GMPA)    = (V3evs - Sample.V) / V3evs 
  d/dt(A4GMPA) =  Qr  * C3GMPA - Qr * C4GMPA  +  C4MPA * N4 * fumMPA * CL4MPA 
  f(A4GMPA)    = (V4evs - Sample.V) / V4evs 
  d/dt(A5GMPA) =  - PappGMPA* SA23 * convfactor * NPERM *(C5GMPA - C3GMPA * fumGMPA) - PappGMPA* SA23 * convfactor * NPERM *(ErG*C5GMPA - C2GMPA * fum) + CL3MPA * N3 * C5MPA 
  
  # cell desnity in the apical side (enterocytes)
  cell_den = N3/V2s
  
  #  C2ObsMM=C2MM *exp(CEps.C2.MM)
  C2ObsMPA=C2MPA *exp(CEps.C2.MPA)
  C3ObsMPA=C3MPA *exp(CEps.C3.MPA)
  C4ObsMPA=C4MPA *exp(CEps.C4.MPA)
  
  C2ObsGMPA=C2GMPA *exp(CEps.C2.GMPA)
  C3ObsGMPA=C3GMPA *exp(CEps.C3.GMPA)
  C4ObsGMPA=C4GMPA *exp(CEps.C4.GMPA)
  
})

# Inizialitzation flow 
timeQ = 1 # Start flow after x min
time.inc = 2880 # incubation time in min
Conc.ini=10 # uM Mofetil
o.target.AP <- c(30,60,120,480,2880)            #114
o.target.BS <- c(30,60,120,480,1440,2880)        #293
o.target.LM <- c(120,240,480,1440,2880)
Sample.V=15    # sampling volume uL bas and liver
Sample.Va=c(15)   # sampling volume uL apical
tkev.a = c(0, 0.01, 0.03, 0.05)
tkev = c(0, 0.01, 0.03, 0.05)

t12tab = 1000000000  # 2880   # t1/2 (first oreder process) in min

#Event table                           
qd <- et(seq(0, 2880, by = 1))

sim.rep.pop.all = data.frame()
for (i in 1:length(tkev.a)){

  # paramters
  par<- c(tQ=90,                 
          V2=325, 
          V3=1506,             
          V4=1394,
          PappMPA=546,          # Evaluated from  exp (nm/s)
          SA23=0.33, 
          convfactor=0.006,
          NPERM=2,
          fum=1,
          fumMPA=0.38,
          fumGMPA=0.70,
          CL3MPA=17,         # Evaluated from  exp (uL/min/Mio cells)
          CL3MM = 13,         # Parameter from EXP5 all simul. fitting
          V5=2.6,
          tN3=0.45,           # From mean of gut-liver data Exp 5
          tN4=0.30,          # From mean of gut-liver data Exp 5
          Er=3.0,                # From Cn Bio model MG with input
          ErG=3.1,
          PappGMPA=0.35,
          CL4MPA= 26,
          Sample.V = Sample.V,
          Sample.Va = Sample.Va,
          kev.a= tkev.a[i],
          tkev = tkev[i],   # From exp 5 (mean of values)
          Al = 7,
          Ab = 0,
          ktab = log(2)/t12tab)
  
  sim.ev.MM <- et(seq(0,tail(o.target.AP,1), by = 1)) %>%
    et(cmt="C2ObsMM", time=0, evid = 1, amt=Conc.ini*unname(par["V2"])) %>%      
    et(cmt="C2ObsMPA", time=0, evid = 1, amt=0) %>%    
    et(cmt="C2ObsGMPA", time=0, evid = 1, amt=0) %>%   
    et(cmt="C3ObsMPA", time=0, evid = 1, amt=0) %>%    
    et(cmt="C3ObsGMPA", time=0, evid = 1, amt=0) %>%   
    et(cmt="C4ObsMPA", time=0, evid = 1, amt=0) %>%    
    et(cmt="C4ObsGMPA", time=0, evid = 1, amt=0) %>%   
    
    # Sampling volume in the Apical 
    et(cmt="Vmsa", time=o.target.AP[1], evid = 1, amt=Sample.Va) %>%      
    et(cmt="Vmsa", time=o.target.AP[2], evid = 1, amt=Sample.Va) %>% 
    et(cmt="Vmsa", time=o.target.AP[3], evid = 1, amt=Sample.Va) %>% 
    et(cmt="Vmsa", time=o.target.AP[4], evid = 1, amt=Sample.Va) %>% 
    et(cmt="Vmsa", time=o.target.AP[5], evid = 1, amt=Sample.Va) %>% 
    
    # Sampling amount of substrate in Apical
    et(cmt="A2MM", time=o.target.AP[1], evid=6, amt=1) %>%
    et(cmt="A2MM", time=o.target.AP[2], evid=6, amt=1) %>%
    et(cmt="A2MM", time=o.target.AP[3], evid=6, amt=1) %>%
    et(cmt="A2MM", time=o.target.AP[4], evid=6, amt=1) %>%
    et(cmt="A2MM", time=o.target.AP[5], evid=6, amt=1) %>%
    
    # Sampling amount of MPA in Apical
    et(cmt="A2MPA", time=o.target.AP[1], evid=6, amt=1) %>%
    et(cmt="A2MPA", time=o.target.AP[2], evid=6, amt=1) %>%
    et(cmt="A2MPA", time=o.target.AP[3], evid=6, amt=1) %>%
    et(cmt="A2MPA", time=o.target.AP[4], evid=6, amt=1) %>%
    et(cmt="A2MPA", time=o.target.AP[5], evid=6, amt=1) %>%
    
    # Sampling amount of MPA in Apical
    et(cmt="A2GMPA", time=o.target.AP[1], evid=6, amt=1) %>%
    et(cmt="A2GMPA", time=o.target.AP[2], evid=6, amt=1) %>%
    et(cmt="A2GMPA", time=o.target.AP[3], evid=6, amt=1) %>%
    et(cmt="A2GMPA", time=o.target.AP[4], evid=6, amt=1) %>%
    et(cmt="A2GMPA", time=o.target.AP[5], evid=6, amt=1) %>%
    
    # # Sampling volume in the basolateral 
    # et(cmt="Vmsb", time=o.target.BS[1], evid = 1, amt=Sample.V) %>%      
    # et(cmt="Vmsb", time=o.target.BS[2], evid = 1, amt=Sample.V) %>% 
    # et(cmt="Vmsb", time=o.target.BS[3], evid = 1, amt=Sample.V) %>% 
    # et(cmt="Vmsb", time=o.target.BS[4], evid = 1, amt=Sample.V) %>% 
    # et(cmt="Vmsb", time=o.target.BS[5], evid = 1, amt=Sample.V) %>% 
    # et(cmt="Vmsb", time=o.target.BS[6], evid = 1, amt=Sample.V) %>% 
    
    # # Sampling amount of MPA in Basolateral 
    # et(cmt="A3MPA", time=o.target.BS[1], evid=6, amt=1) %>%
  # et(cmt="A3MPA", time=o.target.BS[2], evid=6, amt=1) %>%
  # et(cmt="A3MPA", time=o.target.BS[3], evid=6, amt=1) %>%
  # et(cmt="A3MPA", time=o.target.BS[4], evid=6, amt=1) %>%
  # et(cmt="A3MPA", time=o.target.BS[5], evid=6, amt=1) %>%
  # et(cmt="A3MPA", time=o.target.BS[6], evid=6, amt=1) %>%
  # 
  # # Sampling amount of GMPA in basolateral
  # et(cmt="A3GMPA", time=o.target.BS[1], evid=6, amt=1) %>%
  # et(cmt="A3GMPA", time=o.target.BS[2], evid=6, amt=1) %>%
  # et(cmt="A3GMPA", time=o.target.BS[3], evid=6, amt=1) %>%
  # et(cmt="A3GMPA", time=o.target.BS[4], evid=6, amt=1) %>%
  # et(cmt="A3GMPA", time=o.target.BS[5], evid=6, amt=1) %>%
  # et(cmt="A3GMPA", time=o.target.BS[6], evid=6, amt=1) %>%
  # 
  # Sampling volume in the liver 
  et(cmt="Vmsl", time=o.target.LM[1], evid = 1, amt=Sample.V) %>%      
    et(cmt="Vmsl", time=o.target.LM[2], evid = 1, amt=Sample.V) %>% 
    et(cmt="Vmsl", time=o.target.LM[3], evid = 1, amt=Sample.V) %>% 
    et(cmt="Vmsl", time=o.target.LM[4], evid = 1, amt=Sample.V) %>% 
    et(cmt="Vmsl", time=o.target.LM[5], evid = 1, amt=Sample.V) %>% 
    et(cmt="Vmsl", time=o.target.BS[1], evid = 1, amt=Sample.V) %>%
    et(cmt="Vmsl", time=o.target.BS[2], evid = 1, amt=Sample.V) %>%
    et(cmt="Vmsl", time=o.target.BS[3], evid = 1, amt=Sample.V) %>%
    et(cmt="Vmsl", time=o.target.BS[4], evid = 1, amt=Sample.V) %>%
    et(cmt="Vmsl", time=o.target.BS[5], evid = 1, amt=Sample.V) %>%
    et(cmt="Vmsl", time=o.target.BS[6], evid = 1, amt=Sample.V) %>%
    
    # Sampling amount of MPA in Liver 
    et(cmt="A4MPA", time=o.target.LM[1], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.LM[2], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.LM[3], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.LM[4], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.LM[5], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.BS[1], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.BS[2], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.BS[3], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.BS[4], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.BS[5], evid=6, amt=1) %>%
    et(cmt="A4MPA", time=o.target.BS[6], evid=6, amt=1) %>%
    
    # Sampling amount of GMPA in Liver
    et(cmt="A4GMPA", time=o.target.LM[1], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.LM[2], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.LM[3], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.LM[4], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.LM[5], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.BS[1], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.BS[2], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.BS[3], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.BS[4], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.BS[5], evid=6, amt=1) %>%
    et(cmt="A4GMPA", time=o.target.BS[6], evid=6, amt=1) 
  
  # Omega and sigma matrix 
  
  omega.matrix.mod1.pop <- lotri({
    c(eta.N4 ~ 0.000015^2, eta.N3 ~ 0.000015^2, eta.Q ~ 0.00002^2, eta.kev ~ 0.000015^2) })
  
  
  sigma.matrix.mod1.pop <- lotri( CEps.C2.MPA ~ 0.15^2, CEps.C3.MPA ~ 0.15^2, CEps.C4.MPA ~ 0.15^2, 
                                  CEps.C2.GMPA ~ 0.15^2, CEps.C3.GMPA ~ 0.15^2, CEps.C4.GMPA ~ 0.15^2)
  
  initialization <- c(A2MPA=Conc.ini*unname(par["V2"]), A3MPA = 0, A4MPA = 0, A5MPA =0, A3GMPA = 0, A4GMPA = 0, A5GMPA =0, Vtab=0)
  
  #qd <- et(seq(0,time.inc, by = 1)) %>%
  #   et(cmt="A2MM", time=0, evid = 1, amt=Conc.ini*325) 
  
  qd <- data.frame(qd)
  
  
  set.seed(123)
  sim.rep.pop <- rxSolve(gut.liver.MM,par,sim.ev.MM,inits= initialization, omega=omega.matrix.mod1.pop, addDosing=TRUE,sigma=sigma.matrix.mod1.pop, nSub=2)
  
  sim.rep.pop<-data.frame(sim.rep.pop)
  sim.rep.pop.plot.2 <- sim.rep.pop %>% filter(sim.rep.pop$evid == 0)

  # define sampling volume columns
  sim.rep.pop.plot.2$VsamplingA = tkev.a[i]
  sim.rep.pop.all= rbind(sim.rep.pop.all, sim.rep.pop.plot.2)
}

sim.rep.pop.all= sim.rep.pop.all %>% filter(sim.rep.pop.all$time>2)
sim.rep.pop.all= sim.rep.pop.all %>% filter(sim.rep.pop.all$sim.id==1)

x_ax = theme(axis.text.x=element_text(size=15),axis.title=element_text(size=15))
y_ax = theme(axis.text.y=element_text(size=15),axis.title=element_text(size=15))


PV2<-ggplot(sim.rep.pop.all,aes(time/60,V2s)) + geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Volume Apical Side (uL)") + xlab("Time (h)") + theme_bw()  + ylim(30,325) +theme(legend.position = "none") + geom_hline(yintercept=325-325*0.15, linetype="dashed", color = "red") +x_ax + y_ax
Pcell_d<-ggplot(sim.rep.pop.all,aes(time/60,cell_den*1000)) +  geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Cells density (Mio cells/mL)") + xlab("Time (h)") + theme_bw()  +  theme(legend.position = "none") + geom_hline(yintercept=(sim.rep.pop.all[1,"cell_den"]+sim.rep.pop.all[1,"cell_den"]*0.15)*1000, linetype="dashed", color = "red") +x_ax + y_ax 
C2<-ggplot(sim.rep.pop.all,aes(time/60,C2MPA)) +  geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Conc. Sub. Apical Side (uM)") + xlab("Time (h)") + theme_bw()     +  theme(legend.position=c(0.75, 0.55))+ labs(color='Media Evap. (uL/min)') + scale_y_continuous(trans='log10') +x_ax + y_ax + theme(legend.text=element_text(size=rel(1.3))) +theme(legend.title =element_text(size=rel(1.3)))
A2<-ggplot(sim.rep.pop.all,aes(time/60,A2MPA)) +  geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Amount Sub. Apical Side (pmol)") + xlab("Time (h)") + theme_bw() + theme(legend.position = "none")+ scale_y_continuous(trans='log10') +x_ax + y_ax


C3<-ggplot(sim.rep.pop.all,aes(time/60,C3MPA)) + geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Conc. Sub. Basol. Side (uM)") + xlab("Time (h)") + theme_bw()   +theme(legend.position = "none") 
C4<-ggplot(sim.rep.pop.all,aes(time/60,C4MPA)) +  geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Conc. Sub. Liver MPS (uM)") + xlab("Time (h)") + theme_bw()  +  theme(legend.position = "none") 
#C2<-ggplot(sim.rep.pop.all,aes(time/60,C2MPA)) +  geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Conc. Apical Side (uM)") + xlab("Time (h)") + theme_bw()     +  theme(legend.position=c(0.75, 0.55))+ labs(color='Media Evap. (uL/min)')
#A2<-ggplot(sim.rep.pop.all,aes(time/60,A2MPA)) +  geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Amount Apical Side (pmol)") + xlab("Time (h)") + theme_bw() + theme(legend.position = "none")
C2M<-ggplot(sim.rep.pop.all,aes(time/60, C2GMPA)) + geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Conc. Met. Apical Side (uM)") + xlab("Time (h)") + theme_bw()  +theme(legend.position = "none") 
C3M<-ggplot(sim.rep.pop.all,aes(time/60, C3GMPA)) +  geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Conc. Met Basol. Side (uM)") + xlab("Time (h)") + theme_bw()  +  theme(legend.position = "none") 
C4M<-ggplot(sim.rep.pop.all,aes(time/60, C4GMPA)) +  geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Conc. Met. Liver MPS (uM)") + xlab("Time (h)") + theme_bw()  +  theme(legend.position = "none") 


#A2_zoom<-ggplot(sim.rep.pop.all,aes(time/60,A2MPA)) +  geom_line(aes(color=as.factor(VsamplingA)), size=1) + ylab("Amount Apical Side (pmol)") + xlab("Time (h)") + theme_bw() + xlim(0,6)

grid.arrange(PV2,Pcell_d,C2,A2)

grid.arrange(C2,C3,C4,C2M,C3M,C4M,ncol = 3)


max(sim.rep.pop.all$C2GMPA)


##### TEST A SMART WAY TO WRITE EV TABLE ######

prova <- et(seq(0,tail(o.target.AP,1), by = 1))
  
for (i in 1:5){
  prova <- prova %>%
  et(cmt="C2ObsMM", time=i-1, evid = 1, amt=Conc.ini*unname(par["V2"]),ii=0)      
 
}
