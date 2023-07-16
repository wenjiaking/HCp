library(ggplot2)
library(cowplot)
library(rlist)
prep.plots=function(orig.libs) {
  libs=sapply(orig.libs,function(x) unlist(x))
  mean.libs=apply(libs,1,mean)
  mse.libs=apply(libs,1,function(x) mean((log10(x)-log10(mean(x)))^2))
  cv.libs=apply(libs,1,function(x) sd(x)/mean(x))
  return(list(mean.libs,mse.libs,cv.libs))
}

#### Generate Figrue 2 #########
q.vals50=exp(seq(log(10^3),log(10^7),length.out = 200))
ZY500.HCorig=readRDS("ZY500.HCorig.rds")
XH500=readRDS("XH500.rds")
libsHC500=list.load("libs500.RData")
pvalsHC500=prep.plots(libsHC500)
ZY500=readRDS("ZY500.rds")
subHC=seq(1,200,length.out = 15)

compHC500=data.frame(x=rep(log(q.vals50[subHC]),8),
                     y=c(-log10(XH500)[subHC],-log10(abs(ZY500.HCorig))[subHC],-log10(XH500)[subHC],-log10(pvalsHC500[[1]])[subHC],
                         -log10(abs(ZY500.HCorig))[subHC],-log10(pvalsHC500[[1]])[subHC],-log10(abs(ZY500))[subHC],-log10(pvalsHC500[[1]])[subHC]),
                     sd=c(rep(NA,length(subHC)),rep(NA,length(subHC)),rep(NA,length(subHC)),sqrt(pvalsHC500[[2]][subHC]),
                          rep(NA,length(subHC)),sqrt(pvalsHC500[[2]][subHC]),rep(NA,length(subHC)),sqrt(pvalsHC500[[2]][subHC])),
                     group=factor(c(rep(c("LIN","SetTest"),each=length(subHC)),rep(c("LIN","IS"),each=length(subHC)),
                                    rep(c("SetTest","IS"),each=length(subHC)),rep(c("HC.mod","IS"),each=length(subHC))),
                                  levels = c("LIN","SetTest","IS","HC.mod")),
                     comp=factor(c(rep("LIN vs SetTest",2*length(subHC)),rep("LIN vs IS",2*length(subHC)),
                                   rep("SetTest vs IS",2*length(subHC)),rep("HC.mod vs IS",2*length(subHC))),
                                 levels=c("LIN vs SetTest","LIN vs IS","SetTest vs IS","HC.mod vs IS"),
                                 labels = c(expression("BL-SetTest"),expression(paste("IS"["CE"],"-BL")),expression(paste("IS"["CE"],"-SetTest")),expression(paste("IS"["CE"],"-MST")))))
f1a=ggplot(compHC500,aes(x=x,y=y,group=group))+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.5),alpha=0.6)+
  geom_point(aes(color=group,shape=group))+
  scale_color_manual(name="",values = c(LIN="#6B8E23",SetTest="#D2691E",IS="#808080",HC.mod="blue"),
                     labels=c("Barnett-Lin (BL)",expression("SetTest"),expression("*IS"["CE"]),"*MST"))+
  ylab("-log10 Estimated P-value") +
  xlab("log statistic") +
  facet_grid(~comp,labeller = "label_parsed")+
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"))+
  #scale_color_discrete(labels=c("Barnett-Lin",expression("SetTest.RT"^"HC"),"IS",expression("SetTest.RT"^"HC")))+
  scale_shape_manual(values = c(LIN=15,SetTest=17,IS=3,HC.mod=19),labels=c("Barnett-Lin (BL)",expression("SetTest"),expression("*IS"["CE"]),"*MST"))+
  theme(legend.text.align = 0,
        legend.text=element_text(size=8),
        legend.title = element_text(face="bold"),
        strip.text = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 8),
        axis.title = element_text(face="bold"))+
  ggtitle("(A)")

mixgradlibs_modHC500=list.load("mixgradlibs_modHC500.RData")
mixgrad_modHC500=prep.plots(mixgradlibs_modHC500)
ZY_modtrunc500=readRDS("ZY_modtrunc500.rds")
ZY_trunc500=readRDS("ZY_trunc500.rds")

qvals=exp(seq(log(30),17,length=100))
subTHC=seq(1,100,length.out = 15)
compTHC500=data.frame(x=rep(log(qvals)[subTHC],4),y=c(-log10(abs(ZY_trunc500))[subTHC],-log10(mixgrad_modHC500[[1]])[subTHC],-log10(abs(ZY_modtrunc500))[subTHC],-log10(mixgrad_modHC500[[1]])[subTHC]),
                      sd=c(rep(NA,length(subTHC)),sqrt(mixgrad_modHC500[[2]][subTHC]),rep(NA,length(subTHC)),sqrt(mixgrad_modHC500[[2]][subTHC])),
                      group=factor(c(rep(c("SetTest","IS"),each=length(subTHC)),rep(c("HC.mod","IS"),each=length(subTHC))),
                                   levels = c("SetTest","IS","HC.mod")),
                      comp=factor(c(rep("SetTest vs IS",2*length(subTHC)),rep("HC.mod vs IS",2*length(subTHC))),
                                  levels=c("SetTest vs IS","HC.mod vs IS"),
                                  labels = c(expression(paste("IS"["CE"],"-SetTest")),expression(paste("IS"["CE"],"-MST")))))

f1b=ggplot(compTHC500,aes(x=x,y=y,group=group))+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.5),alpha=0.6)+
  geom_point(aes(color=group,shape=group))+
  scale_color_manual(name="",values = c(SetTest="#D2691E",IS="#808080",HC.mod="blue"),
                     labels=c(expression("SetTest"),expression("*IS"["CE"]),"*MST"))+
  ylab("-log10 Estimated P-value") +
  xlab("log statistic") +
  facet_grid(~comp,labeller = "label_parsed")+
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"))+
  #scale_color_discrete(labels=c("BL",expression("SetTest.RT"^"HC"),"IS",expression("SetTest.RT"^"HC")))+
  scale_shape_manual(values = c(SetTest=17,IS=3,HC.mod=19),labels=c(expression("SetTest"),expression("*IS"["CE"]),"*MST"))+
  theme(legend.text.align = 0,
        legend.position = "bottom",
        legend.text=element_text(size=6),
        legend.title = element_text(face="bold"),
        strip.text = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 8),
        axis.title = element_text(face="bold"))+
  ggtitle("(B)")

mse=function(x) {
  return(mean((log10(x)-log10(mean(x,na.rm = T)))^2,na.rm = T))
}

mixgrad5libs50_MHC100=readRDS("mixgrad5libs50N7_MHC100.rds")
ZY_MHC100=readRDS("ZY_MHC100.rds")
mixgrad5libs50_TMHC100=readRDS("mixgrad5libs50N7_TMHC100.rds")
ZY_TMHC100=readRDS("ZY_TMHC100.rds")
subInd=seq(1,100,length.out = 15)

MHC100prop_data=data.frame(q=rep(seq(log(3),log(13),length=100)[subInd],2),
                           p=c(-log10(apply(mixgrad5libs50_MHC100,1,mean)),-log10(ZY_MHC100)[subInd]),
                           sd=c(sqrt(apply(mixgrad5libs50_MHC100, 1, mse)),rep(NA,length(subInd))),
                           method=factor(c(rep("IS",length(subInd)),rep("Analytic",length(subInd))),levels = c("IS","Analytic")),
                           comp=factor(rep("Analytic vs IS",2*length(subInd)),
                                       labels = c(expression(paste("IS"["CE"],"-MST")))))
TMHC100prop_data=data.frame(q=rep(seq(log(3),log(13),length=100)[subInd],2),
                            p=c(-log10(apply(mixgrad5libs50_TMHC100,1,mean)),-log10(ZY_TMHC100)[subInd]),
                            sd=c(sqrt(apply(mixgrad5libs50_TMHC100, 1, mse)),rep(NA,length(subInd))),
                            method=factor(c(rep("IS",length(subInd)),rep("Analytic",length(subInd))),levels = c("IS","Analytic")),
                            comp=factor(rep("Analytic vs IS",2*length(subInd)),
                                        labels = c(expression(paste("IS"["CE"],"-MST")))))

f1c=ggplot(MHC100prop_data)+
  geom_point(aes(x=q,y=p,group=method,color=method,shape=method))+
  geom_errorbar(aes(x=q,ymin=p-sd,ymax=p+sd,group=method,color=method))+
  #geom_line(aes(x=q,y=p,group=method,color=method))+
  scale_color_manual(name="",values = c(IS="#808080",Analytic="blue"),
                     labels=c(expression("*IS"["CE"]),"*MST"))+
  scale_shape_manual(values = c(IS=3,Analytic=19),
                     labels=c(expression("*IS"["CE"]),"*MST"))+
  guides(col=guide_legend(title=""),
         shape=guide_legend(title=""))+
  ylab("-log10 Estimated P-value") +
  xlab("log statistic") +
  facet_grid(~comp,labeller = "label_parsed")+
  theme(legend.text.align = 0,
        legend.position = "bottom",
        legend.text=element_text(size=6),
        legend.title = element_text(face="bold"),
        strip.text = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 8),
        axis.title = element_text(face="bold"))+
  ggtitle("(C)")


f1d=ggplot(TMHC100prop_data)+
  geom_point(aes(x=q,y=p,group=method,color=method,shape=method))+
  geom_errorbar(aes(x=q,ymin=p-sd,ymax=p+sd,group=method,color=method))+
  #geom_line(aes(x=q,y=p,group=method,color=method))+
  scale_color_manual(name="",values = c(IS="#808080",Analytic="blue"),
                     labels=c(expression("*IS"["CE"]),"*MST"))+
  scale_shape_manual(values = c(IS=3,Analytic=19),
                     labels=c(expression("*IS"["CE"]),"*MST"))+
  guides(col=guide_legend(title=""),
         shape=guide_legend(title=""))+
  ylab("-log10 Estimated P-value") +
  xlab("log statistic") +
  facet_grid(~comp,labeller = "label_parsed")+
  theme(legend.text.align = 0,
        legend.position = "bottom",
        legend.text=element_text(size=6),
        legend.title = element_text(face="bold"),
        strip.text = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 8),
        axis.title = element_text(face="bold"))+
  ggtitle("(D)")

pdf("Figure2.pdf")
plot_grid(f1a,plot_grid(f1b,f1c,f1d,rel_widths = c(3,2,2),nrow = 1),nrow = 2)
dev.off()
##### Generate Figure 3 #######
q_sm=exp(seq(log(0.5),log(3.7),length.out = 50))
mixgradSplib_MHC100=readRDS("mixgradSplib_MHC100.rds") #10^7 sampling size
MCLpMHC100lib=readRDS("MCLpMHC100lib.rds")
LSapproSpMHC100=readRDS("LSapproSpMHC100.rds")
LSapproLpMHC100=readRDS("LSapproLpMHC100.rds")
ZY_SpMHC100=readRDS("ZY_SpMHC100.rds")
ZY_LpMHC100=readRDS("ZY_LpMHC100.rds")
qs=c(exp(seq(log(3.8),log(13),length.out = 150))[seq(80,150,length.out = 12)],exp(seq(log(13.5),log(14.5),length.out=3)))[1:11]
ql=c(q_sm[seq(17,42,length.out = 10)],exp(seq(log(3),log(13),length=100)[seq(1,100,length.out = 15)][1:5]))
MHC100Sp_dat=data.frame(q=c(qs,qs-0.1,qs+0.1),
                        p=-log10(c(ZY_SpMHC100[1:11],LSapproSpMHC100[1:11],apply(mixgradSplib_MHC100,1,mean)[1:11])),
                        sd=c(rep(NA,22),sqrt(apply(mixgradSplib_MHC100,1,mse))[1:11]),
                        method=factor(c(rep("MST",11),rep("Li-Siegmund",11),rep("IS",11)),levels = c("MST","Li-Siegmund","IS")))
MHC100Lp_dat=data.frame(q=c(ql,ql-0.1,ql+0.1),
                        p=c(ZY_LpMHC100,LSapproLpMHC100,apply(MCLpMHC100lib,1,mean)),
                        sd=c(rep(NA,30),sqrt(apply(MCLpMHC100lib,1,function(x) mean((x-mean(x))^2,na.rm = T)))),
                        method=factor(c(rep("MST",15),rep("Li-Siegmund",15),rep("MC",15)),levels = c("MST","Li-Siegmund","MC")))

mixgradSplib_TMHC100=readRDS("mixgradSplib_TMHC100.rds")
MCLpTMHC100lib=readRDS("MCLpTMHC100lib.rds")
mod_SpTMHC100=readRDS("mod_SpTMHC100.rds")
mod_LpTMHC100=readRDS("mod_LpTMHC100.rds")
LSapproSpTMHC100=readRDS("LSapproSpTMHC100.rds")
LSapproLpTMHC100=readRDS("LSapproLpTMHC100.rds")

TMHC100Sp_dat=data.frame(q=c(qs,qs+0.1),
                         p=-log10(c(LSapproSpTMHC100[1:11],apply(mixgradSplib_TMHC100,1,mean)[1:11])),
                         sd=c(rep(NA,11),sqrt(apply(mixgradSplib_TMHC100,1,mse))[1:11]),
                         method=factor(c(rep("Li-Siegmund",11),rep("IS",11)),levels = c("MST","Li-Siegmund","IS")))
TMHC100Lp_dat=data.frame(q=c(ql,ql+0.1),
                         p=c(LSapproLpTMHC100,apply(MCLpTMHC100lib,1,mean)),
                         sd=c(rep(NA,15),sqrt(apply(MCLpTMHC100lib,1,function(x) mean((x-mean(x))^2,na.rm = T)))),
                         method=factor(c(rep("Li-Siegmund",15),rep("MC",15)),levels = c("MST","Li-Siegmund","MC")))

MHC100Lpeval.plot=ggplot(MHC100Lp_dat)+
  geom_point(aes(x=q,y=p,color=method,shape=method),alpha=0.7)+
  geom_errorbar(aes(x=q,ymin=p-sd,ymax=p+sd,color=method))+
  scale_colour_manual("",
                      breaks = c("MST", "Li-Siegmund", "MC"),
                      values = c("red","#D2691E","blue"),
                      labels=c("MST", "Li-Siegmund","Monte Carlo"))+
  scale_shape_manual("",
                     breaks = c("MST", "Li-Siegmund", "MC"),
                     values = c(18,17,19),
                     labels=c("MST", "Li-Siegmund","Monte Carlo")) +
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"),
         linetype=guide_legend(title="Method"))+
  ylab("Estimated p-value") +
  xlab("statistic") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "none",
        legend.text = element_text(face="bold",size = 8),
        legend.title = element_text(face="bold"),
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 8),
        axis.title = element_text(face="bold"))+
  ggtitle("(A)")


MHC100Speval.plot=ggplot(MHC100Sp_dat)+
  geom_point(aes(x=q,y=p,color=method,shape=method),alpha=0.7)+
  geom_errorbar(aes(x=q,ymin=p-sd,ymax=p+sd,color=method))+
  scale_colour_manual("",
                      breaks = c("MST", "Li-Siegmund", "IS"),
                      values = c("red","#D2691E","blue"),
                      labels=c("MST", "Li-Siegmund",expression("IS"["CE"])))+
  scale_shape_manual("",
                     breaks = c("MST", "Li-Siegmund", "IS"),
                     values = c(18,17,19),
                     labels=c("MST", "Li-Siegmund",expression("IS"["CE"]))) +
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"),
         linetype=guide_legend(title="Method"))+
  ylab("-log10 estimated p-value") +
  xlab("statistic") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "none",
        legend.text = element_text(face="bold",size = 8),
        legend.title = element_text(face="bold"),
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 8),
        axis.title = element_text(face="bold"))+
  ggtitle("")

TMHC100Lpeval.plot=ggplot(TMHC100Lp_dat)+
  geom_point(aes(x=q,y=p,color=method,shape=method),alpha=0.7)+
  geom_errorbar(aes(x=q,ymin=p-sd,ymax=p+sd,color=method))+
  scale_colour_manual("",
                      breaks = c("Li-Siegmund", "MC"),
                      values = c("#D2691E","blue"),
                      labels=c( "Li-Siegmund","Monte Carlo"))+
  scale_shape_manual("",
                     breaks = c( "Li-Siegmund", "MC"),
                     values = c(17,19),
                     labels=c( "Li-Siegmund","Monte Carlo")) +
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"),
         linetype=guide_legend(title="Method"))+
  ylab("Estimated p-value") +
  xlab("statistic") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "bottom",
        legend.text = element_text(face="bold",size = 8),
        legend.title = element_text(face="bold"),
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 8),
        axis.title = element_text(face="bold"))+
  ggtitle("(C)")


TMHC100Speval.plot=ggplot(TMHC100Sp_dat)+
  geom_point(aes(x=q,y=p,color=method,shape=method),alpha=0.7)+
  geom_errorbar(aes(x=q,ymin=p-sd,ymax=p+sd,color=method))+
  scale_colour_manual("",
                      breaks = c("Li-Siegmund", "IS"),
                      values = c("#D2691E","blue"),
                      labels=c("Li-Siegmund",expression("IS"["CE"])))+
  scale_shape_manual("",
                     breaks = c("Li-Siegmund", "IS"),
                     values = c(17,19),
                     labels=c("Li-Siegmund",expression("IS"["CE"]))) +
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"),
         linetype=guide_legend(title="Method"))+
  ylab("-log10 estimated p-value") +
  xlab("statistic") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "bottom",
        legend.text = element_text(face="bold",size = 8),
        legend.title = element_text(face="bold"),
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 8),
        axis.title = element_text(face="bold"))+
  ggtitle("(D)")


pdf("Figure3.pdf")
plot_grid(THC2000Lpeval.plot,THC2000Speval.plot,TMHC100Lpeval.plot,TMHC100Speval.plot,ncol = 2)
dev.off()

######Generate Figure 4 (time complexity LIB HC)###############
Kset=c(30,100,500,1000,2000)
nset=c(1,10,100,1000,10000,100000)
HCres_libs=list.load("HCres_libs.RData")
HCres_IS.sub=list.load("HCres_IS.sub.RData")
HCres_IS.lg=list.load("HCres_IS.lg.RData")
HCres_IS.sub[[6]]=HCres_IS.lg
HCres_mod=list.load("HCres_mod.RData")
HCres_ZY=list.load("HCres_ZY.RData")
HCres_LS=list.load("HCres_LS.RData")
HCK500_libs=sapply(1:6, function(i) as.numeric(HCres_libs[[i]][[3]][[1]],units="secs"))
HCK500_libs.fit=lm(log(HCK500_libs)~log(nset))
HCK500_ISsub=sapply(1:6, function(i) as.numeric(HCres_IS.sub[[i]][[3]][[1]],units="secs"))
HCK500_ISsub.fit=lm(log(HCK500_ISsub)~log(nset))
HCK500_mod=sapply(1:6, function(i) as.numeric(HCres_mod[[i]][[3]][[1]],units="secs"))
HCK500_mod.fit=lm(log(HCK500_mod)~log(nset))
HCK500_ZY=sapply(1:6, function(i) as.numeric(HCres_ZY[[i]][[3]][[1]],units="secs"))
HCK500_ZY.fit=lm(log(HCK500_ZY)~log(nset))
HCK500_LS=sapply(1:6, function(i) as.numeric(HCres_LS[[i]][[3]][[1]],units="secs"))
HCK500_LS.fit=lm(log(HCK500_LS)~log(nset))

Kset2=c(50,200,300,400,600,700,800,900,seq(1100,1900,by=100))
Ksetall=c(30,100,Kset2)
HCres2_IS=list.load("HCres2_IS.RData")
HCres2_libs=list.load("HCres2_libs.RData")
HCres2_mod=list.load("HCres2_mod.RData")
HCres2_ZY=list.load("HCres2_ZY.RData")
HCres2_LIN=list.load("HCres2_LIN.RData")
HCres2_LS=list.load("HCres2_LS.RData")

HCres2_IS.t=sapply(HCres2_IS,function(x) as.numeric(x[[1]],units="secs"))
HC_IS.allt=c(sapply(1:2, function(i) as.numeric(HCres_IS.sub[[3]][[i]][[1]],units="secs")),HCres2_IS.t)
HCres2_libs.t=sapply(HCres2_libs,function(x) as.numeric(x[[1]],units="secs"))
HC_libs.allt=c(sapply(1:2, function(i) as.numeric(HCres_libs[[3]][[i]][[1]],units="secs")),HCres2_libs.t)
HCres2_mod.t=sapply(HCres2_mod,function(x) as.numeric(x[[1]],units="secs"))
HC_mod.allt=c(sapply(1:2, function(i) as.numeric(HCres_mod[[3]][[i]][[1]],units="secs")),HCres2_mod.t)
HCres2_ZY.t=sapply(HCres2_ZY,function(x) as.numeric(x[[1]],units="secs"))
HC_ZY.allt=c(sapply(1:2, function(i) as.numeric(HCres_ZY[[3]][[i]][[1]],units="secs")),HCres2_ZY.t)
HCres2_LIN.t=sapply(HCres2_LIN,function(x) as.numeric(x[[1]],units="secs"))
HCres2_LS.t=sapply(HCres2_LS,function(x) as.numeric(x[[1]],units="secs"))
HC_LS.allt=c(sapply(1:2, function(i) as.numeric(HCres_LS[[3]][[i]][[1]],units="secs")),HCres2_LS.t)

polyd1_IS <- lm(formula = HC_IS.allt~poly(Ksetall,1))
summary(polyd1_IS)
polyd2_IS <- lm(formula = HC_IS.allt~poly(Ksetall,2))
summary(polyd2_IS) #second order term is not significant
polyd0_libs <- lm(formula = HC_libs.allt~1)
summary(polyd0_libs)
polyd1_libs <- lm(formula = HC_libs.allt~poly(Ksetall,1))
summary(polyd1_libs) #first order term is just *
polyd2_libs <- lm(formula = HC_libs.allt~poly(Ksetall,2))
summary(polyd2_libs)
polyd2_mod <- lm(formula = HC_mod.allt~poly(Ksetall,2))
summary(polyd2_mod)
polyd3_mod <- lm(formula = HC_mod.allt~poly(Ksetall,3))
summary(polyd3_mod) #second order term
polyd2_ZY <- lm(formula = HC_ZY.allt~poly(Ksetall,2))
summary(polyd2_ZY)
polyd3_ZY <- lm(formula = HC_ZY.allt~poly(Ksetall,3))
summary(polyd3_ZY)


polyd3_LIN <- lm(formula = HCres2_LIN.t~poly(Kset2,3))
summary(polyd3_LIN) #all term significant ***
polyd4_LIN <- lm(formula = HCres2_LIN.t~poly(Kset2,4))
summary(polyd4_LIN) #the 4th order is not signigicant

polyd1_LS <- lm(formula = HC_LS.allt~poly(Ksetall,1))
summary(polyd1_LS)
polyd2_LS <- lm(formula = HC_LS.allt~poly(Ksetall,2))
summary(polyd2_LS)

HCtime_K.plot=ggplot()+
  geom_point(aes(Ksetall, HC_libs.allt,col = "UFI",shape="UFI"), cex = 2)+
  stat_smooth(method = "lm", formula = y~1, aes(Ksetall, polyd0_libs$fitted.values, col = "UFI",lty="UFI")) +
  geom_point(aes(Ksetall, HC_IS.allt,col = "CEIS",shape="CEIS"), cex = 2)+
  stat_smooth(method = "lm", formula = y~poly(x,1), aes(Ksetall, polyd1_IS$fitted.values, col = "CEIS",lty="CEIS")) +
  geom_point(aes(Ksetall, HC_ZY.allt,col = "SetTest",shape="SetTest"), cex = 2)+
  stat_smooth(method = "lm", formula = y~poly(x,2), aes(Ksetall, polyd2_ZY$fitted.values, col = "SetTest",lty="SetTest")) +
  geom_point(aes(Ksetall, HC_mod.allt,col = "mod",shape="mod"), cex = 2)+
  stat_smooth(method = "lm", formula = y~poly(x,2), aes(Ksetall, polyd2_mod$fitted.values, col = "mod",lty="mod")) +
  geom_point(aes(Ksetall, HC_LS.allt,col = "LS",shape="LS"), cex = 2)+
  stat_smooth(method = "lm", formula = y~poly(x,2), aes(Ksetall, polyd1_LS$fitted.values, col = "LS",lty="LS")) +
  scale_colour_manual("",
                      breaks = c("UFI", "CEIS", "SetTest","mod","LS"),
                      values = c("red","#808080","#D2691E","blue","pink"),
                      labels=c(expression(paste("UFI: ","t=",0.00037)),expression(paste("IS"["CE"],": ","t=",2357.94*K+889.03)),
                               expression(paste("SetTest: ", "t=",32.68*K^2+294.13*K+89.07)),expression(paste("MST: ","t=",58.90*K^2+281.67*K+85.57)),
                               expression(paste("Li-Siegmund: ","t=",0.050*K+0.017)))) +
  scale_shape_manual("",
                     breaks = c("UFI", "CEIS", "SetTest","mod","LS"),
                     values = c(18,3,17,19,4),
                     labels=c(expression(paste("UFI: ","t=",0.00037)),expression(paste("IS"["CE"],": ","t=",2357.94*K+889.03)),
                              expression(paste("SetTest: ", "t=",32.68*K^2+294.13*K+89.07)),expression(paste("MST: ","t=",58.90*K^2+281.67*K+85.57)),
                              expression(paste("Li-Siegmund: ","t=",0.050*K+0.017)))) +
  scale_linetype_manual("",
                        breaks = c("UFI", "CEIS", "SetTest","mod","LS"),
                        values = c(6,2,3,5,4),
                        labels=c(expression(paste("UFI: ","t=",0.00037)),expression(paste("IS"["CE"],": ","t=",2357.94*K+889.03)),
                                 expression(paste("SetTest: ", "t=",32.68*K^2+294.13*K+89.07)),expression(paste("MST: ","t=",58.90*K^2+281.67*K+85.57)),
                                 expression(paste("Li-Siegmund: ","t=",0.050*K+0.017)))) +
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"),
         linetype=guide_legend(title="Method"))+
  ylab("Time (secs)") +
  xlab("K") +
  ggtitle("(B)")+
  theme(legend.text.align = 0,
        legend.position = c(0.3,0.83),
        #legend.position = c(0.27,0.74),
        legend.background = element_rect(fill='transparent'),
        legend.key = element_rect(fill="transparent"),
        legend.text=element_text(size=8,face = "bold"),
        legend.title = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 8),
        axis.title = element_text(face="bold",size=10),
        axis.text = element_text(face="bold",size=8))

HCtime_K.fullplot=ggplot()+
  geom_point(aes(Ksetall[-(1:5)], HC_libs.allt[-(1:5)],col = "UFI",shape="UFI"), cex = 2)+
  stat_smooth(method = "lm", formula = y~1, aes(Ksetall[-(1:5)], polyd0_libs$fitted.values[-(1:5)], col = "UFI",lty="UFI")) +
  geom_point(aes(Ksetall[-(1:5)], HC_IS.allt[-(1:5)],col = "CEIS",shape="CEIS"), cex = 2)+
  stat_smooth(method = "lm", formula = y~poly(x,1), aes(Ksetall[-(1:5)], polyd1_IS$fitted.values[-(1:5)], col = "CEIS",lty="CEIS")) +
  geom_point(aes(Ksetall[-(1:5)], HC_ZY.allt[-(1:5)],col = "SetTest",shape="SetTest"), cex = 2)+
  stat_smooth(method = "lm", formula = y~poly(x,2), aes(Ksetall[-(1:5)], polyd2_ZY$fitted.values[-(1:5)], col = "SetTest",lty="SetTest")) +
  geom_point(aes(Ksetall[-(1:5)], HC_mod.allt[-(1:5)],col = "mod",shape="mod"), cex = 2)+
  stat_smooth(method = "lm", formula = y~poly(x,2), aes(Ksetall[-(1:5)], polyd2_mod$fitted.values[-(1:5)], col = "mod",lty="mod")) +
  geom_point(aes(Kset2[-(1:3)], HCres2_LIN.t[-(1:3)],col = "BL",shape="BL"), cex = 2)+
  stat_smooth(method = "lm", formula = y~poly(x,2), aes(Kset2[-(1:3)], polyd3_LIN$fitted.values[-(1:3)], col = "BL",lty="BL")) +
  geom_point(aes(Ksetall[-(1:5)], HC_LS.allt[-(1:5)],col = "LS",shape="LS"), cex = 2)+
  stat_smooth(method = "lm", formula = y~poly(x,2), aes(Ksetall[-(1:5)], polyd1_LS$fitted.values[-(1:5)], col = "LS",lty="LS")) +
  scale_colour_manual("",
                      breaks = c("UFI", "CEIS", "SetTest","mod","BL","LS"),
                      values = c("red","#808080","#D2691E","blue","#52854C","pink"),
                      labels=c(expression(paste("UFI: ","t=",0.00037)),expression(paste("IS"["CE"],": ","t=",2357.94*K+889.03)),
                               expression(paste("SetTest: ", "t=",32.68*K^2+294.13*K+89.07)),expression(paste("MST: ","t=",58.90*K^2+281.67*K+85.57)),
                               expression(paste("Barnett-Lin: ", "t=",25370.68*K^3+159128.50*K^2+427810.00*K+117663.06)),
                               expression(paste("Li-Siegmund: ","t=",0.050*K+0.017)))) +
  scale_shape_manual("",
                     breaks = c("UFI", "CEIS", "SetTest","mod","BL","LS"),
                     values = c(18,3,17,19,7,4),
                     labels=c(expression(paste("UFI: ","t=",0.00037)),expression(paste("IS"["CE"],": ","t=",2357.94*K+889.03)),
                              expression(paste("SetTest: ", "t=",32.68*K^2+294.13*K+89.07)),expression(paste("MST: ","t=",58.90*K^2+281.67*K+85.57)),
                              expression(paste("Barnett-Lin: ", "t=",25370.68*K^3+159128.50*K^2+427810.00*K+117663.06)),
                              expression(paste("Li-Siegmund: ","t=",0.050*K+0.017)))) +
  scale_linetype_manual("",
                        breaks = c("UFI", "CEIS", "SetTest","mod","BL","LS"),
                        values = c(6,2,3,5,1,4),
                        labels=c(expression(paste("UFI: ","t=",0.00037)),expression(paste("IS"["CE"],": ","t=",2357.94*K+889.03)),
                                 expression(paste("SetTest: ", "t=",32.68*K^2+294.13*K+89.07)),expression(paste("MST: ","t=",58.90*K^2+281.67*K+85.57)),
                                 expression(paste("Barnett-Lin: ", "t=",25370.68*K^3+159128.50*K^2+427810.00*K+117663.06)),
                                 expression(paste("Li-Siegmund: ","t=",0.050*K+0.017)))) +
  guides(col=guide_legend(title="Method",override.aes=list(fill=NA)),
         shape=guide_legend(title="Method"),
         linetype=guide_legend(title="Method"))+
  ylab("Time (secs)") +
  ggtitle("(A)")+
  xlab("K") +
  theme(legend.text.align = 0,
        legend.position = c(0.4,0.8),
        #legend.position = c(0.3,0.7),
        legend.background = element_rect(fill='transparent'),
        legend.key = element_rect(fill="transparent"),
        legend.text=element_text(size=8,face = "bold"),
        legend.title = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 8),
        axis.title = element_text(face="bold",size=10),
        axis.text = element_text(face="bold",size=8))


HCtime_repeat.plot=ggplot()+
  geom_point(aes(log(nset), log(HCK500_libs),col = "UFT",shape="UFT"), cex = 2)+
  geom_line(aes(log(nset), log(HCK500_libs),col = "UFT",lty="UFT"))+
  geom_point(aes(log(nset), log(HCK500_ISsub),col = "CEIS",shape="CEIS"), cex = 2)+
  geom_line(aes(log(nset), log(HCK500_ISsub),col = "CEIS",lty="CEIS"))+
  geom_point(aes(log(nset), log(HCK500_ZY),col = "SetTest",shape="SetTest"), cex = 2)+
  geom_line(aes(log(nset), log(HCK500_ZY),col = "SetTest",lty="SetTest"))+
  geom_point(aes(log(nset), log(HCK500_mod),col = "mod",shape="mod"), cex = 2)+
  geom_line(aes(log(nset), log(HCK500_mod),col = "mod",lty="mod"))+
  geom_abline(aes(intercept = 0,slope = 1,col="baseline",lty="baseline"),size=0.5)+
  geom_abline(aes(intercept=2.1,slope=1,col="baseline",lty="baseline"),size=0.5)+
  geom_abline(aes(intercept=-11,slope=1,col="baseline",lty="baseline"),size=0.5)+
  geom_point(aes(log(nset), log(HCK500_mod),shape="baseline"),col="white",alpha=0)+
  scale_colour_manual("",
                      breaks = c("UFT", "CEIS", "SetTest","mod","baseline"),
                      values = c("red","#808080","#D2691E","blue","black"),
                      labels=c("UFI", expression("IS"["CE"]), "SetTest","MST","baseline (slope=1)"))  +
  scale_shape_manual("",
                     breaks = c("UFT", "CEIS", "SetTest","mod","baseline"),
                     values = c(18,3,17,19,7),
                     labels=c("UFI", expression("IS"["CE"]), "SetTest","MST","baseline (slope=1)")) +
  scale_linetype_manual("",
                        breaks = c("UFT", "CEIS", "SetTest","mod","baseline"),
                        values = c(6,2,3,5,1),
                        labels=c("UFI", expression("IS"["CE"]), "SetTest","MST","baseline (slope=1)")) +
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"),
         linetype=guide_legend(title="Method"))+
  ylim(-15,13)+
  ylab("log time (secs)") +
  xlab("log number of repeats") +
  ggtitle("(B)")+
  theme(legend.text.align = 0,
        legend.text=element_text(size=8),
        legend.title = element_text(face="bold"),
        legend.position = c(0.8, 0.2),
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 8),
        axis.title = element_text(face="bold",size=10),
        axis.text = element_text(face="bold",size=8))

HCvalues=c()
HClib_values=c()
for (i in 1:5) {
  temp.ind=which(-log10(HCres_libs[[5]][[i]][[2]])<12)[1:100]
  temp.val=HCres_mod[[5]][[i]][[2]][temp.ind]
  #temp.val=sapply(temp.val,function(x) ifelse(x<=0,10^(-13),x))
  HCvalues=c(HCvalues,temp.val)
  HClib_values=c(HClib_values,HCres_libs[[5]][[i]][[2]][temp.ind])
}

HCaccuracy.dat=data.frame(x=HCvalues,value=HClib_values,
                          K=rep(Kset,each=100))
HClib_accu.plot=ggplot(HCaccuracy.dat)+
  geom_point(aes(x=-log10(x),y=-log10(value),
                 color=factor(paste0("K=",K),levels =paste0("K=",Kset)),shape=factor(paste0("K=",K),levels =paste0("K=",Kset))),cex=2)+
  geom_abline(intercept=0,slope = 1)+
  #geom_line(aes(x=-log10(x),y=-log10(value),
  #               color=factor(methods,levels =c("HC.mod","SetTest","IS")),lty=factor(methods,levels =c("HC.mod","SetTest","IS"))))+

  guides(col=guide_legend(title="Dimension"),
         shape=guide_legend(title="Dimension"))+
  ylab("-log10 p-value by UFI") +
  xlab("-log10 analytic p-value (truth)") +
  ggtitle("(C)")+
  theme(legend.text.align = 0,
        legend.position = c(0.85, 0.25),
        legend.text=element_text(size=8,face = "bold"),
        legend.background = element_rect(fill='transparent'),
        legend.title = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 6),
        axis.title = element_text(face="bold",size=10),
        axis.text = element_text(face="bold",size=8))

pdf("Figure4.pdf",height=10,width = 10)
plot_grid(plot_grid(HCtime_K.fullplot,HCtime_K.plot,nrow = 1),HClib_accu.plot,nrow = 2)
dev.off()

############Real application: COVID19 desease survilence##############
#data is Provided by Johns Hopkins Center for Systems Science and Engineering (CSSE)
UScovid=read.csv("data/time_series_covid19_confirmed_US.csv",header = F)
Day=as.Date(as.character(UScovid[1,12:815]),format = "%m/%d/%y")
UScovid=UScovid[-1,-c(1,2,3,4,5,8)] #3342 counties in total and 59 states, 804 days.
colnames(UScovid)=c("county","state","Lat","Long","location",as.character(Day))
saveRDS(UScovid,"UScovid.rds")

#On "2022-03-10" (750) the number of significant counties reach the minimum,
# compared with "2021-03-10"(385), "2020-03-10" (20)
# match with PA part: "2020-07-11" (143),"2021-07-11"(508),"2022-04-02"(773)
library(maps)
library(dplyr)
#install.packages("tidyverse_1.3.1.tar.gz", repos = NULL, type ="source")
#library(socviz)
HCpvals=readRDS("US.HCpvals.rds")
HCresUS=matrix(HCpvals,nrow=3342)
UScovid=readRDS("UScovid.rds")
us_states <- map_data("state")
us_counties <- map_data("county")
colnames(us_counties)=c("long","lat","group","order","state","county" )
#"2021-03-22", "2021-07-20","2022-02-15"
countiesP=cbind(UScovid[,1:5],HCresUS[,c(3,11,14)]) #number of outbreak counties is 3297 3059 2441 by <0.05, and 2378 1594 1735 by <0.01
countiesP[is.na(countiesP)]=1
countiesP$state=tolower(countiesP$state)
countiesP$county=tolower(countiesP$county)
us_counties$county[is.na(match(us_counties$county,countiesP$county))]="unassigned"
us_countiesP=left_join(us_counties,countiesP)
map1=ggplot(data = us_countiesP,aes(x = long, y = lat,group = group,fill=-log10(as.numeric(`1`)))) +
  geom_polygon(color = "gray90", size = 0.05) + coord_equal()+
  scale_fill_gradient(low = "white", high = "#CB454A")+
  labs(fill = "-log10 p-value")+
  labs(title = "03/22/2021 - 04/20/2021")
map2=ggplot(data = us_countiesP,aes(x = long, y = lat,group = group,fill=-log10(as.numeric(`2`)))) +
  geom_polygon(color = "gray90", size = 0.05) + coord_equal()+
  scale_fill_gradient(low = "white", high = "#CB454A")+
  labs(fill = "-log10 p-value")+
  labs(title = "11/17/2021 - 12/16/2021")
map3=ggplot(data = us_countiesP,aes(x = long, y = lat,group = group,fill=-log10(as.numeric(`3`)))) +
  geom_polygon(color = "gray90", size = 0.05) + coord_equal()+
  scale_fill_gradient(low = "white", high = "#CB454A")+
  labs(fill = "-log10 p-value")+
  labs(title = "02/15/2022 - 03/16/2022")

pdf("Figure5.pdf")
plot_grid(map1,map2,map3,nrow = 3)
dev.off()

##### Generate Supplementary Figure A1 #####
unbias_plotdat=function(dat) {
  tempdat=lapply(dat,unlist)
  tempout=sapply(tempdat,function(x) {
    temp_p=x
    mean_p=mean(temp_p,na.rm=TRUE)
    mse_p=mean((temp_p-mean_p)^2,na.rm=TRUE)
    cv_p=sd(temp_p,na.rm = TRUE)/mean(temp_p,na.rm = TRUE)
    return(c(mean_p,mse_p,cv_p))
  })
  return(tempout)
}
idxHC100N4_fix=list.load("idxHC100N4_fix.RData")
idxHC100N5_fix=list.load("idxHC100N5_fix.RData")
idxHC100N6_fix=list.load("idxHC100N6_fix.RData")
idxTHC100N4_fix=list.load("idxTHC100N4_fix.RData")
idxTHC100N5_fix=list.load("idxTHC100N5_fix.RData")
idxTHC100N6_fix=list.load("idxTHC100N6_fix.RData")
idxHC100N4_fixdat=unbias_plotdat(idxHC100N4_fix)
idxHC100N5_fixdat=unbias_plotdat(idxHC100N5_fix)
idxHC100N6_fixdat=unbias_plotdat(idxHC100N6_fix)
idxTHC100N4_fixdat=unbias_plotdat(idxTHC100N4_fix)
idxTHC100N5_fixdat=unbias_plotdat(idxTHC100N5_fix)
idxTHC100N6_fixdat=unbias_plotdat(idxTHC100N6_fix)

idxHC100.fixdat=as.data.frame(rbind(t(idxHC100N4_fixdat),t(idxHC100N5_fixdat),t(idxHC100N6_fixdat)))
colnames(idxHC100.fixdat)=c("mean_p","mse_p","cv_p")
idxHC100.fixdat$idx=factor(c(0.3,0.5,1,2),levels = c("0.3","0.5","1","2"))
idxHC100.fixdat$SampleSize=factor(rep(c("10^4","10^5","10^6"),each=4),levels =c("10^4","10^5","10^6"))
idxHC100.fixdat$truth=rep(readRDS("idxHC100analy.rds"),12)

idxTHC100.fixdat=as.data.frame(rbind(t(idxTHC100N4_fixdat),t(idxTHC100N5_fixdat),t(idxTHC100N6_fixdat)))
colnames(idxTHC100.fixdat)=c("mean_p","mse_p","cv_p")
idxTHC100.fixdat$idx=factor(c(0.3,0.5,1,2),levels = c("0.3","0.5","1","2"))
idxTHC100.fixdat$SampleSize=factor(rep(c("10^4","10^5","10^6"),each=4),levels =c("10^4","10^5","10^6"))
idxTHC100.fixdat$truth=rep(readRDS("idxTHC100analy.rds"),12) #1.014893e-10

idxHC100.fixplot=ggplot(idxHC100.fixdat) +
  geom_point(position=position_dodge(width=0.3),aes(x=idx,y=mean_p,color=SampleSize))+
  geom_errorbar(position=position_dodge(width=0.3),aes(x=idx,ymax=mean_p+sqrt(mse_p),ymin=mean_p-sqrt(mse_p),width = 0.2,color=SampleSize))+
  geom_point(aes(x=idx,y=truth),size=2,shape=17)+
  ylab("Estimated P-value") +
  xlab("gradient degree (d)") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "none",
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 8),
        axis.title = element_text(face="bold"))+
  ggtitle("(B)")

idxTHC100.fixplot=ggplot(idxTHC100.fixdat) +
  geom_point(position=position_dodge(width=0.3),aes(x=idx,y=mean_p,color=SampleSize))+
  geom_errorbar(position=position_dodge(width=0.3),aes(x=idx,ymax=mean_p+sqrt(mse_p),ymin=mean_p-sqrt(mse_p),width = 0.2,color=SampleSize))+
  geom_point(aes(x=idx,y=truth),size=2,shape=17)+
  ylab("Estimated P-value") +
  xlab("gradient degree (d)") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "none",
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 8),
        axis.title = element_text(face="bold"))+
  ggtitle("(C)")

propMHC100N4_fix=list.load("propMHC100N4_fix.RData")
propMHC100N5_fix=list.load("propMHC100N5_fix.RData")
propMHC100N6_fix=list.load("propMHC100N6_fix.RData")
propTMHC100N4_fix=list.load("propTMHC100N4_fix.RData")
propTMHC100N5_fix=list.load("propTMHC100N5_fix.RData")
propTMHC100N6_fix=list.load("propTMHC100N6_fix.RData")
propMHC100N4_fixdat=unbias_plotdat(propMHC100N4_fix)
propMHC100N5_fixdat=unbias_plotdat(propMHC100N5_fix)
propMHC100N6_fixdat=unbias_plotdat(propMHC100N6_fix)
propTMHC100N4_fixdat=unbias_plotdat(propTMHC100N4_fix)
propTMHC100N5_fixdat=unbias_plotdat(propTMHC100N5_fix)
propTMHC100N6_fixdat=unbias_plotdat(propTMHC100N6_fix)
ZY_MHC100=readRDS("ZY_MHC100.rds")
ZY_TMHC100=readRDS("ZY_TMHC100.rds")
subInd=seq(1,100,length.out = 15)
propMHC100.fixdat=as.data.frame(rbind(t(propMHC100N4_fixdat),t(propMHC100N5_fixdat),t(propMHC100N6_fixdat)))
colnames(propMHC100.fixdat)=c("mean_p","mse_p","cv_p")
propMHC100.fixdat$prop=factor(c("-Inf","0.01","0.1","0.2","0.5"),levels = c("-Inf","0.01","0.1","0.2","0.5"))
propMHC100.fixdat$SampleSize=factor(rep(c("10^4","10^5","10^6"),each=5),levels =c("10^4","10^5","10^6"))
propMHC100.fixdat$truth=rep(ZY_MHC100[subInd][14],15)

propTMHC100.fixdat=as.data.frame(rbind(t(propTMHC100N4_fixdat),t(propTMHC100N5_fixdat),t(propTMHC100N6_fixdat)))
colnames(propTMHC100.fixdat)=c("mean_p","mse_p","cv_p")
propTMHC100.fixdat$prop=factor(c("-Inf","0.01","0.1","0.2","0.5"),levels = c("-Inf","0.01","0.1","0.2","0.5"))
propTMHC100.fixdat$SampleSize=factor(rep(c("10^4","10^5","10^6"),each=5),levels =c("10^4","10^5","10^6"))
propTMHC100.fixdat$truth=rep(ZY_TMHC100[subInd][14],15)
propMHC100.fixplot=ggplot(propMHC100.fixdat) +
  geom_point(position=position_dodge(width=0.3),aes(x=prop,y=mean_p,color=SampleSize))+
  geom_errorbar(position=position_dodge(width=0.3),aes(x=prop,ymax=mean_p+sqrt(mse_p),ymin=mean_p-sqrt(mse_p),width = 0.2,color=SampleSize))+
  geom_point(aes(x=prop,y=truth),size=2,shape=17)+
  ylab("Estimated P-value") +
  xlab("Mixture Proportion (prop)") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "bottom",
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 8),
        axis.title = element_text(face="bold"))+
  ggtitle("(D)")

propTMHC100.fixplot=ggplot(propTMHC100.fixdat) +
  geom_point(position=position_dodge(width=0.3),aes(x=prop,y=mean_p,color=SampleSize))+
  geom_errorbar(position=position_dodge(width=0.3),aes(x=prop,ymax=mean_p+sqrt(mse_p),ymin=mean_p-sqrt(mse_p),width = 0.2,color=SampleSize))+
  geom_point(aes(x=prop,y=truth),size=2,shape=17)+
  ylab("Estimated P-value") +
  xlab("Mixture Proportion (prop)") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "bottom",
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 8),
        axis.title = element_text(face="bold"))+
  ggtitle("(E)")
#calcualte the HC statistics for many sets of pvalues
libs.FHC50.Beta=list.load("libs.FHC50.Beta.RData")
libs.FHC50.Gauss=list.load("libs.FHC50.Gauss.RData")
libs.FHC50.mixGrad=list.load("libs.FHC50.mixGrad.RData")
libs.FHC50.Weibull_s5=list.load("libs.FHC50.Weibull_s5.RData")
libs.FHC50.Weibull_s1=list.load("libs.FHC50.Weibull_s1.RData")
libs.FHC50.Exp=list.load("libs.FHC50.Exp.RData")
FHC50.Beta=prep.plots(libs.FHC50.Beta)
FHC50.Gauss=prep.plots(libs.FHC50.Gauss)
FHC50.mixGrad=prep.plots(libs.FHC50.mixGrad)
FHC50.Weibull_s5=prep.plots(libs.FHC50.Weibull_s5)
FHC50.Weibull_s1=prep.plots(libs.FHC50.Weibull_s1)
FHC50.Exp=prep.plots(libs.FHC50.Exp)
q.vals50=exp(seq(log(10^3),log(10^7),length.out = 200))
analyticP=readRDS("FHC50.analyticP.rds")
familycomp.dat=data.frame(est=c(FHC50.Weibull_s1[[1]],FHC50.Weibull_s5[[1]],FHC50.Beta[[1]],FHC50.Exp[[1]],FHC50.Gauss[[1]],FHC50.mixGrad[[1]]),
                          truth=rep(analyticP,6),
                          q=rep(q.vals50[c(50,80,100,150)],6),
                          mse=c(FHC50.Weibull_s1[[2]],FHC50.Weibull_s5[[2]],FHC50.Beta[[2]],FHC50.Exp[[2]],FHC50.Gauss[[2]],FHC50.mixGrad[[2]]),
                          Family=rep(c("Weibull (shape=1)","Weibull (shape=5)","Beta","Exponential","Gaussian","Gaussion Mixture (gradient)"),each=4))
familycomp.plot=ggplot(familycomp.dat) +
  geom_point(position=position_dodge(width=0.3),aes(x=log(q),y=-log10(est),color=Family,shape=Family))+
  geom_errorbar(position=position_dodge(width=0.3),aes(x=log(q),ymax=-log10(est)+sqrt(mse),ymin=-log10(est)-sqrt(mse),width = 0.2,color=Family))+
  geom_point(aes(x=log(q)-0.05,y=-log10(truth)),shape=17)+
  ylab("-log10 Estimated P-value") +
  xlab("log statistic") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "right",
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 8),
        axis.title = element_text(face="bold"))+
  ggtitle("(A)")
pdf("FigureA1.pdf")
plot_grid(familycomp.plot,plot_grid(idxHC100.fixplot,idxTHC100.fixplot,ncol=2),plot_grid(propMHC100.fixplot,propTMHC100.fixplot,ncol=2),
          nrow = 3,rel_heights=c(4,4,5))
dev.off()


#######Generate Supplementary Figure A2##########
prep.plots=function(orig.libs) {
  libs=sapply(orig.libs,function(x) unlist(x))
  mean.libs=apply(libs,1,mean)
  mse.libs=apply(libs,1,function(x) mean((log10(x)-log10(mean(x)))^2))
  cv.libs=apply(libs,1,function(x) sd(x)/mean(x))
  return(list(mean.libs,mse.libs,cv.libs))
}
libsHC500=list.load("data/libs500.RData")
pvalsHC500=prep.plots(libsHC500)
ZY500.HCorig123=readRDS("data/ZY500.HCorig123.rds")
ZY500.HCorig124=readRDS("data/ZY500.HCorig124.rds")
ZY500.HCorig234=readRDS("data/ZY500.HCorig234.rds")
#orig134 NaN
ZY500.HCorig23=readRDS("data/ZY500.HCorig23.rds")
ZY500.HCorig24=readRDS("data/ZY500.HCorig24.rds")
ZY500.HCorig12=readRDS("data/ZY500.HCorig12.rds")
ZY500.HCorig13=readRDS("data/ZY500.HCorig13.rds")
#orig14/orig34 NaN and orig13 negative
#orig4 dimension will generate NaN
ZY500.HCorig3=readRDS("data/ZY500.HCorig3.rds")
ZY500.HCorig2=readRDS("data/ZY500.HCorig2.rds")
ZY500.HCorig1=readRDS("data/ZY500.HCorig1.rds")
#M1 is necessary: without M1 (orig123/234/134/23/34/13/3) will generate negative value
q.vals50=exp(seq(log(10^3),log(10^7),length.out = 200))

compZYmod1_CE500=data.frame(x=rep(log(q.vals50),2),y=c(-log10(abs(ZY500.HCorig234)),-log10(pvalsHC500[[1]])),
                            sd=c(rep(NA,length(q.vals50)),pvalsHC500[[2]]),group=as.factor(rep(c("M2","IS"),each=length(q.vals50))))
plotZYmod1_CE500=ggplot(data =compZYmod1_CE500,aes(x=x,group=group))+
  geom_point(aes(y=y,color=group,shape=group),alpha=0.5)+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.15),alpha=0.2)+
  ylab("-log10 Estimated P-value") +
  xlab("log statistics") +
  ggtitle("HC K=500: M2 vs IS") +
  theme_bw() +
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.25,.8),plot.title = element_text(hjust = 0.5))


#compZYmod4_CE500=data.frame(x=rep(log(q.vals50),2),y=c(sapply(ZY500.HCorig123,function(x) ifelse(x<=0,5,-log10(x))),-log10(pvalsHC500[[1]])),
#                            sd=c(rep(NA,length(q.vals50)),pvalsHC500[[2]]),group=as.factor(rep(c("M4","CE"),each=length(q.vals50))))
compZYmod4_CE500=data.frame(x=rep(log(q.vals50),2),y=c(-log10(abs(ZY500.HCorig123)),-log10(pvalsHC500[[1]])),
                            sd=c(rep(NA,length(q.vals50)),pvalsHC500[[2]]),group=as.factor(rep(c("M4","IS"),each=length(q.vals50))))
plotZYmod4_CE500=ggplot(data =compZYmod4_CE500,aes(x=x,group=group))+
  geom_point(aes(y=y,color=group,shape=group),alpha=0.5)+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.15),alpha=0.2)+
  ylab("-log10 Estimated P-value") +
  xlab("log statistics") +
  ggtitle("M4 vs IS") +
  theme_bw() +
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.25,.8),plot.title = element_text(hjust = 0.5))

#compZYmod14_CE500=data.frame(x=rep(log(q.vals50),2),y=c(sapply(ZY500.HCorig23,function(x) ifelse(x<=0,5,-log10(x))),-log10(pvalsHC500[[1]])),
#                            sd=c(rep(NA,length(q.vals50)),pvalsHC500[[2]]),group=as.factor(rep(c("M1&M4","CE"),each=length(q.vals50))))
compZYmod14_CE500=data.frame(x=rep(log(q.vals50),2),y=c(-log10(abs(ZY500.HCorig23)),-log10(pvalsHC500[[1]])),
                             sd=c(rep(NA,length(q.vals50)),pvalsHC500[[2]]),group=as.factor(rep(c("M2&M4","IS"),each=length(q.vals50))))
plotZYmod14_CE500=ggplot(data =compZYmod14_CE500,aes(x=x,group=group))+
  geom_point(aes(y=y,color=group,shape=group),alpha=0.5)+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.15),alpha=0.2)+
  ylab("-log10 Estimated P-value") +
  xlab("log statistics") +
  ggtitle("M2&M4 vs IS") +
  theme_bw() +
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.35,.8),plot.title = element_text(hjust = 0.5))

#compZYmod24_CE500=data.frame(x=rep(log(q.vals50),2),y=c(sapply(ZY500.HCorig13,function(x) ifelse(x<=0,5,-log10(x))),-log10(pvalsHC500[[1]])),
#                             sd=c(rep(NA,length(q.vals50)),pvalsHC500[[2]]),group=as.factor(rep(c("M2&M4","CE"),each=length(q.vals50))))
compZYmod24_CE500=data.frame(x=rep(log(q.vals50),2),y=c(-log10(abs(ZY500.HCorig13)),-log10(pvalsHC500[[1]])),
                             sd=c(rep(NA,length(q.vals50)),pvalsHC500[[2]]),group=as.factor(rep(c("M3&M4","IS"),each=length(q.vals50))))
plotZYmod24_CE500=ggplot(data =compZYmod24_CE500,aes(x=x,group=group))+
  geom_point(aes(y=y,color=group,shape=group),alpha=0.5)+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.15),alpha=0.2)+
  ylab("-log10 Estimated P-value") +
  xlab("log statistics") +
  ggtitle("M3&M4 vs IS") +
  theme_bw() +
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.35,.8),plot.title = element_text(hjust = 0.5))

#compZYmod124_CE500=data.frame(x=rep(log(q.vals50),2),y=c(sapply(ZY500.HCorig3,function(x) ifelse(x<=0,5,-log10(x))),-log10(pvalsHC500[[1]])),
#                             sd=c(rep(NA,length(q.vals50)),pvalsHC500[[2]]),group=as.factor(rep(c("M1&M2&M4","CE"),each=length(q.vals50))))
compZYmod124_CE500=data.frame(x=rep(log(q.vals50),2),y=c(-log10(abs(ZY500.HCorig3)),-log10(pvalsHC500[[1]])),
                              sd=c(rep(NA,length(q.vals50)),pvalsHC500[[2]]),group=as.factor(rep(c("M2&M3&M4","IS"),each=length(q.vals50))))
plotZYmod124_CE500=ggplot(data =compZYmod124_CE500,aes(x=x,group=group))+
  geom_point(aes(y=y,color=group,shape=group),alpha=0.5)+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.15),alpha=0.2)+
  ylab("-log10 Estimated P-value") +
  xlab("log statistics") +
  ggtitle("M2&M3&M4 vs IS") +
  theme_bw() +
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.35,.8),plot.title = element_text(hjust = 0.5))
pdf("FigureA2_1.pdf",width = 12,height=5)
grid.arrange(plotZYmod1_CE500, plotZYmod4_CE500,plotZYmod14_CE500, plotZYmod24_CE500,plotZYmod124_CE500,nrow = 1)
dev.off()

#With M3, M1 is necessary: without M1 (orig124/14/12) will be 0
#compZYmod3_CE500=data.frame(x=rep(log(q.vals50),2),y=c(sapply(ZY500.HCorig124,function(x) ifelse(x<=0,5,-log10(x))),-log10(pvalsHC500[[1]])),
#                              sd=c(rep(NA,length(q.vals50)),pvalsHC500[[2]]),group=as.factor(rep(c("M3","CE"),each=length(q.vals50))))
compZYmod3_CE500=data.frame(x=rep(log(q.vals50),2),y=c(-log10(abs(ZY500.HCorig124)),-log10(pvalsHC500[[1]])),
                            sd=c(rep(NA,length(q.vals50)),pvalsHC500[[2]]),group=as.factor(rep(c("M1","IS"),each=length(q.vals50))))

plotZYmod3_CE500=ggplot(data =compZYmod3_CE500,aes(x=x,group=group))+
  geom_point(aes(y=y,color=group,shape=group),alpha=0.5)+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.15),alpha=0.2)+
  ylab("-log10 Estimated P-value") +
  xlab("log statistics") +
  ggtitle("M1 vs IS") +
  theme_bw() +
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.25,.8),plot.title = element_text(hjust = 0.5))

#compZYmod34_CE500=data.frame(x=rep(log(q.vals50),2),y=c(sapply(ZY500.HCorig12,function(x) ifelse(x<=0,5,-log10(x))),-log10(pvalsHC500[[1]])),
#                            sd=c(rep(NA,length(q.vals50)),pvalsHC500[[2]]),group=as.factor(rep(c("M3&M4","CE"),each=length(q.vals50))))
compZYmod34_CE500=data.frame(x=rep(log(q.vals50),2),y=c(-log10(abs(ZY500.HCorig12)),-log10(pvalsHC500[[1]])),
                             sd=c(rep(NA,length(q.vals50)),pvalsHC500[[2]]),group=as.factor(rep(c("M1&M4","IS"),each=length(q.vals50))))
plotZYmod34_CE500=ggplot(data =compZYmod34_CE500,aes(x=x,group=group))+
  geom_point(aes(y=y,color=group,shape=group),alpha=0.5)+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.15),alpha=0.2)+
  ylab("-log10 Estimated P-value") +
  xlab("log statistics") +
  ggtitle("M1&M4 vs IS") +
  theme_bw() +
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.25,.8),plot.title = element_text(hjust = 0.5))
#compZYmod234_CE500=data.frame(x=rep(log(q.vals50),2),y=c(sapply(ZY500.HCorig1,function(x) ifelse(x<=0,5,-log10(x))),-log10(pvalsHC500[[1]])),
#                             sd=c(rep(NA,length(q.vals50)),pvalsHC500[[2]]),group=as.factor(rep(c("M2&M3&M4","CE"),each=length(q.vals50))))
compZYmod234_CE500=data.frame(x=rep(log(q.vals50),2),y=c(-log10(abs(ZY500.HCorig1)),-log10(pvalsHC500[[1]])),
                              sd=c(rep(NA,length(q.vals50)),pvalsHC500[[2]]),group=as.factor(rep(c("M1&M3&M4","IS"),each=length(q.vals50))))
plotZYmod234_CE500=ggplot(data =compZYmod234_CE500,aes(x=x,group=group))+
  geom_point(aes(y=y,color=group,shape=group),alpha=0.5)+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.15),alpha=0.2)+
  ylab("-log10 Estimated P-value") +
  xlab("log statistics") +
  ggtitle("M1&M3&M4 vs IS") +
  theme_bw() +
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.35,.8),plot.title = element_text(hjust = 0.5))
#M3&M2 (orig14) are NaN


#with M1, M3, M4 and M2 are not necessary: orig24 works

#compZYmod13_CE500=data.frame(x=rep(log(q.vals50),2),y=c(sapply(ZY500.HCorig24,function(x) ifelse(x<=0,5,-log10(x))),-log10(pvalsHC500[[1]])),
#                             sd=c(rep(NA,length(q.vals50)),pvalsHC500[[2]]),group=as.factor(rep(c("M1&M3","CE"),each=length(q.vals50))))
compZYmod13_CE500=data.frame(x=rep(log(q.vals50),2),y=c(-log10(abs(ZY500.HCorig24)),-log10(pvalsHC500[[1]])),
                             sd=c(rep(NA,length(q.vals50)),pvalsHC500[[2]]),group=as.factor(rep(c("M1&M2","IS"),each=length(q.vals50))))
plotZYmod13_CE500=ggplot(data =compZYmod13_CE500,aes(x=x,group=group))+
  geom_point(aes(y=y,color=group,shape=group),alpha=0.5)+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.15),alpha=0.2)+
  ylab("-log10 Estimated P-value") +
  xlab("log statistics") +
  ggtitle("M1&M2 vs IS") +
  theme_bw() +
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.25,.8),plot.title = element_text(hjust = 0.5))
#with M1,M3 and M2, M4 is necessary: without M4 (orig4) will be NaN
pdf("FigureA2_2.pdf",height = 5,width=10)
grid.arrange(plotZYmod3_CE500, plotZYmod34_CE500,plotZYmod234_CE500,plotZYmod13_CE500,nrow = 1)
dev.off()

#######Generate Supplementary Figure A3##########
mixgradlibs_modHC500=list.load("data/mixgradlibs_modHC500.RData")
mixgrad_modHC500=prep.plots(mixgradlibs_modHC500)
ZYorig123_MHC500=readRDS("data/ZYorig123_MHC500.rds")
ZYorig124_MHC500=readRDS("data/ZYorig124_MHC500.rds")
ZYorig234_MHC500=readRDS("data/ZYorig234_MHC500.rds")
#orig134 NaN
ZYorig23_MHC500=readRDS("data/ZYorig23_MHC500.rds")
ZYorig24_MHC500=readRDS("data/ZYorig24_MHC500.rds")
ZYorig12_MHC500=readRDS("data/ZYorig12_MHC500.rds")
ZYorig13_MHC500=readRDS("data/ZYorig13_MHC500.rds")
#orig14/orig34 NaN and orig13 negative
#orig4 dimension will generate NaN
ZYorig3_MHC500=readRDS("data/ZYorig3_MHC500.rds")
ZYorig2_MHC500=readRDS("data/ZYorig2_MHC500.rds")
ZYorig1_MHC500=readRDS("data/ZYorig1_MHC500.rds")
#M3 is necessary: without M3 (orig123/234/134/23/34/13/3) will generate negative value

qvals=exp(seq(log(30),17,length=100))
MHCZYmod1_CE500=data.frame(x=rep(log(qvals),2),y=c(-log10(abs(ZYorig234_MHC500)),-log10(mixgrad_modHC500[[1]])),
                           sd=c(rep(NA,length(qvals)),mixgrad_modHC500[[2]]),group=as.factor(rep(c("M2","IS"),each=length(qvals))))
plotMHCZYmod1_CE500=ggplot(data =MHCZYmod1_CE500,aes(x=x,group=group))+
  geom_point(aes(y=y,color=group,shape=group),alpha=0.5)+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.15),alpha=0.2)+
  ylab("-log10 Estimated P-value") +
  xlab("log statistics") +
  ggtitle("THC K=500: M2 vs IS") +
  theme_bw() +
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.25,.8),plot.title = element_text(hjust = 0.5))

MHCZYmod4_CE500=data.frame(x=rep(log(qvals),2),y=c(-log10(abs(ZYorig123_MHC500)),-log10(mixgrad_modHC500[[1]])),
                           sd=c(rep(NA,length(qvals)),mixgrad_modHC500[[2]]),group=as.factor(rep(c("M4","IS"),each=length(qvals))))
plotMHCZYmod4_CE500=ggplot(data =MHCZYmod4_CE500,aes(x=x,group=group))+
  geom_point(aes(y=y,color=group,shape=group),alpha=0.5)+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.15),alpha=0.2)+
  ylab("-log10 Estimated P-value") +
  xlab("log statistics") +
  ggtitle("M4 vs IS") +
  theme_bw() +
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.25,.8),plot.title = element_text(hjust = 0.5))

MHCZYmod14_CE500=data.frame(x=rep(log(qvals),2),y=c(-log10(abs(ZYorig23_MHC500)),-log10(mixgrad_modHC500[[1]])),
                            sd=c(rep(NA,length(qvals)),mixgrad_modHC500[[2]]),group=as.factor(rep(c("M2&M4","IS"),each=length(qvals))))
plotMHCZYmod14_CE500=ggplot(data =MHCZYmod14_CE500,aes(x=x,group=group))+
  geom_point(aes(y=y,color=group,shape=group),alpha=0.5)+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.15),alpha=0.2)+
  ylab("-log10 Estimated P-value") +
  xlab("log statistics") +
  ggtitle("M2&M4 vs IS") +
  theme_bw() +
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.35,.8),plot.title = element_text(hjust = 0.5))

MHCZYmod24_CE500=data.frame(x=rep(log(qvals),2),y=c(-log10(abs(ZYorig13_MHC500)),-log10(mixgrad_modHC500[[1]])),
                            sd=c(rep(NA,length(qvals)),mixgrad_modHC500[[2]]),group=as.factor(rep(c("M3&M4","IS"),each=length(qvals))))
plotMHCZYmod24_CE500=ggplot(data =MHCZYmod24_CE500,aes(x=x,group=group))+
  geom_point(aes(y=y,color=group,shape=group),alpha=0.5)+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.15),alpha=0.2)+
  ylab("-log10 Estimated P-value") +
  xlab("log statistics") +
  ggtitle("M3&M4 vs IS") +
  theme_bw() +
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.35,.8),plot.title = element_text(hjust = 0.5))

MHCZYmod124_CE500=data.frame(x=rep(log(qvals),2),y=c(-log10(abs(ZYorig3_MHC500)),-log10(mixgrad_modHC500[[1]])),
                             sd=c(rep(NA,length(qvals)),mixgrad_modHC500[[2]]),group=as.factor(rep(c("M2&M3&M4","IS"),each=length(qvals))))
plotMHCZYmod124_CE500=ggplot(data =MHCZYmod124_CE500,aes(x=x,group=group))+
  geom_point(aes(y=y,color=group,shape=group),alpha=0.5)+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.15),alpha=0.2)+
  ylab("-log10 Estimated P-value") +
  xlab("log statistics") +
  ggtitle("M2&M3&M4 vs IS") +
  theme_bw() +
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.35,.8),plot.title = element_text(hjust = 0.5))
pdf("FigureA3_1.pdf",width = 12,height=5)
grid.arrange(plotMHCZYmod1_CE500, plotMHCZYmod4_CE500,plotMHCZYmod14_CE500, plotMHCZYmod24_CE500,plotMHCZYmod124_CE500,nrow = 1)
dev.off()

#With M3, M1 is necessary: without M1 (orig124/14/12/orig1) will be 0
MHCZYmod3_CE500=data.frame(x=rep(log(qvals),2),y=c(-log10(abs(ZYorig124_MHC500)),-log10(mixgrad_modHC500[[1]])),
                           sd=c(rep(NA,length(qvals)),mixgrad_modHC500[[2]]),group=as.factor(rep(c("M1","IS"),each=length(qvals))))
plotMHCZYmod3_CE500=ggplot(data =MHCZYmod3_CE500,aes(x=x,group=group))+
  geom_point(aes(y=y,color=group,shape=group),alpha=0.5)+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.15),alpha=0.2)+
  ylab("-log10 Estimated P-value") +
  xlab("log statistics") +
  ggtitle("M1 vs IS") +
  theme_bw() +
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.25,.8),plot.title = element_text(hjust = 0.5))

MHCZYmod34_CE500=data.frame(x=rep(log(qvals),2),y=c(-log10(abs(ZYorig12_MHC500)),-log10(mixgrad_modHC500[[1]])),
                            sd=c(rep(NA,length(qvals)),mixgrad_modHC500[[2]]),group=as.factor(rep(c("M1&M4","IS"),each=length(qvals))))
plotMHCZYmod34_CE500=ggplot(data =MHCZYmod34_CE500,aes(x=x,group=group))+
  geom_point(aes(y=y,color=group,shape=group),alpha=0.5)+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.15),alpha=0.2)+
  ylab("-log10 Estimated P-value") +
  xlab("log statistics") +
  ggtitle("M1&M4 vs IS") +
  theme_bw() +
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.25,.8),plot.title = element_text(hjust = 0.5))

#M3&M2 (orig14) are NaN


#with M1, M3, M4 and M2 are not necessary: orig24 works

MHCZYmod13_CE500=data.frame(x=rep(log(qvals),2),y=c(-log10(abs(ZYorig24_MHC500)),-log10(mixgrad_modHC500[[1]])),
                            sd=c(rep(NA,length(qvals)),mixgrad_modHC500[[2]]),group=as.factor(rep(c("M1&M2","IS"),each=length(qvals))))
plotMHCZYmod13_CE500=ggplot(data =MHCZYmod13_CE500,aes(x=x,group=group))+
  geom_point(aes(y=y,color=group,shape=group),alpha=0.5)+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.15),alpha=0.2)+
  ylab("-log10 Estimated P-value") +
  xlab("log statistics") +
  ggtitle("M1&M2 vs IS") +
  theme_bw() +
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.25,.8),plot.title = element_text(hjust = 0.5))

MHCZYmod234_CE500=data.frame(x=rep(log(qvals),2),y=c(-log10(abs(ZYorig1_MHC500)),-log10(mixgrad_modHC500[[1]])),
                             sd=c(rep(NA,length(qvals)),mixgrad_modHC500[[2]]),group=as.factor(rep(c("M1&M3&M4","IS"),each=length(qvals))))
plotMHCZYmod234_CE500=ggplot(data =MHCZYmod234_CE500,aes(x=x,group=group))+
  geom_point(aes(y=y,color=group,shape=group),alpha=0.5)+
  geom_errorbar(aes(ymin=y-sd,
                    ymax=y+sd,
                    color=group,width=0.15),alpha=0.2)+
  ylab("-log10 Estimated P-value") +
  xlab("log statistics") +
  ggtitle("M1&M3&M4 vs IS") +
  theme_bw() +
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position=c(.3,.8),plot.title = element_text(hjust = 0.5))
#with M1,M3 and M2, M4 is necessary: without M4 (orig4) will be NaN
pdf("FigureA3_2.pdf",height = 5,width=10)
grid.arrange(plotMHCZYmod3_CE500, plotMHCZYmod34_CE500,plotMHCZYmod234_CE500,plotMHCZYmod13_CE500,nrow = 1)
dev.off()

##########Generate Supplementary Figure A4#########

Kset=c(30,100,500,1000,2000)
HCvalues=c()
HClib_values=c()
HClinearlib_value=c()
HCquadraticlib_value=c()
HCquartlib_value=c()
HCres_libs=list.load("HCres_libs.RData")
HCres_mod=list.load("HCres_mod.RData")
HCquadraticfit=list.load("HCquadraticfit.RData")
HClinearfit=list.load("HClinearfit.RData")
HCquartfit=list.load("HCquartlm.RData")
for (i in 1:5) {
  temp.ind=which(-log10(HCres_libs[[5]][[i]][[2]])<12)[1:100]
  temp.val=HCres_mod[[5]][[i]][[2]][temp.ind]
  #temp.val=sapply(temp.val,function(x) ifelse(x<=0,10^(-13),x))
  HCvalues=c(HCvalues,temp.val)
  HClib_values=c(HClib_values,HCres_libs[[5]][[i]][[2]][temp.ind])
  HClinearlib_value=c(HClinearlib_value,HClinearfit[[5]][temp.ind])
  HCquadraticlib_value=c(HCquadraticlib_value,HCquadraticfit[[5]][temp.ind])
  HCquartlib_value=c(HCquartlib_value,HCquartfit[[5]][temp.ind])
}

cubic.plot=ggplot()+
  geom_violin(aes(x="Cubic spline",y=log(HClib_values/HCvalues)),color="black",fill="red")+
  labs(y="log (interpolation/truth)", x= "")+
  geom_hline(aes(yintercept=0),color="black",linetype="dashed")+
  theme_classic()+
  ggtitle("(C)")+
  theme(legend.text.align = 0,
        legend.position ="bottom",
        legend.text=element_text(size=8,face = "bold"),
        legend.background = element_rect(fill='transparent'),
        legend.title = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 6),
        axis.title = element_text(face="bold",size=10),
        axis.text = element_text(face="bold",size=8))

linear.plot=ggplot()+
  geom_violin(aes(x="Linear regression",y=log(HClinearlib_value/HCvalues)),fill="#E69F00")+
  labs(y="log (interpolation/truth)", x= "")+
  geom_hline(aes(yintercept=0),color="black",linetype="dashed")+
  theme_classic()+
  ggtitle("(A)")+
  theme(legend.text.align = 0,
        legend.position ="bottom",
        legend.text=element_text(size=8,face = "bold"),
        legend.background = element_rect(fill='transparent'),
        legend.title = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 6),
        axis.title = element_text(face="bold",size=10),
        axis.text = element_text(face="bold",size=8))

quadratic.plot=ggplot()+
  geom_violin(aes(x="Quadratic regression",y=log(HCquadraticlib_value/HCvalues)),fill="#56B4E9")+
  labs(y="log (interpolation/truth)", x= "")+
  geom_hline(aes(yintercept=0),color="black",linetype="dashed")+
  theme_classic()+
  ggtitle("(B)")+
  theme(legend.text.align = 0,
        legend.position ="bottom",
        legend.text=element_text(size=8,face = "bold"),
        legend.background = element_rect(fill='transparent'),
        legend.title = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 6),
        axis.title = element_text(face="bold",size=10),
        axis.text = element_text(face="bold",size=8))
forthpoly.plot=ggplot()+
  geom_violin(aes(x="Bi-quadratic regression",y=log(HCquartlib_value/HCvalues)),fill="green")+
  labs(y="log (interpolation/truth)", x= "")+
  geom_hline(aes(yintercept=0),color="black",linetype="dashed")+
  theme_classic()+
  ggtitle("(D)")+
  theme(legend.text.align = 0,
        legend.position ="bottom",
        legend.text=element_text(size=8,face = "bold"),
        legend.background = element_rect(fill='transparent'),
        legend.title = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0,size = 14, face = "bold"),
        text = element_text(size = 6),
        axis.title = element_text(face="bold",size=10),
        axis.text = element_text(face="bold",size=8))

pdf("AppendFigureIV.pdf")
plot_grid(linear.plot,quadratic.plot,cubic.plot,forthpoly.plot,nrow = 2)
dev.off()





##########Generate Supplementary Figure A5#########

MCLpFHC2000lib=readRDS("MCLpFHC2000lib.rds")
MCLpFHC1000lib=readRDS("MCLpFHC1000lib.rds")
mixgradSplibs_FHC2000=readRDS("mixgradSplibs_FHC2000.rds")
mixgradSplibs_FHC1000=readRDS("mixgradSplibs_FHC1000.rds")

q.vals.HC=c(exp(seq(log(1),log(9*10^2),length.out = 200)),exp(seq(log(10^3),log(10^7),length.out = 200)))

LSapproSpFHC2000=readRDS("LSapproSpFHC2000.rds")
LSapproSpFHC1000=readRDS("LSapproSpFHC1000.rds")

FHC2000Lp_dat=data.frame(q=c(log(q.vals.HC[seq(1,100,length.out = 15)]),log(q.vals.HC[seq(1,100,length.out = 15)])+0.1),
                         p=c(LSapproSpFHC2000[1:15],apply(MCLpFHC2000lib,1,mean)),
                         sd=c(rep(NA,15),sqrt(apply(MCLpFHC2000lib,1,function(x) mean((x-mean(x))^2,na.rm = T)))),
                         method=factor(c(rep("Li-Siegmund",15),rep("MC",15)),levels = c("Li-Siegmund","MC")))
FHC2000Sp_dat=data.frame(q=log(c(q.vals.HC[seq(110,400,length.out = 25)][1:22],q.vals.HC[seq(110,400,length.out = 25)][1:22]+0.2)),
                         p=-log10(c(LSapproSpFHC2000[16:37],apply(mixgradSplibs_FHC2000,1,mean))),
                         sd=c(rep(NA,22),sqrt(apply(mixgradSplibs_FHC2000,1,mse))),
                         method=factor(c(rep("Li-Siegmund",22),rep("IS",22)),levels = c("Li-Siegmund","IS")))

FHC1000Lp_dat=data.frame(q=c(log(q.vals.HC[seq(1,100,length.out = 15)]),log(q.vals.HC[seq(1,100,length.out = 15)])+0.1),
                         p=c(LSapproSpFHC1000[1:15],apply(MCLpFHC1000lib,1,mean)),
                         sd=c(rep(NA,15),sqrt(apply(MCLpFHC1000lib,1,function(x) mean((x-mean(x))^2,na.rm = T)))),
                         method=factor(c(rep("Li-Siegmund",15),rep("MC",15)),levels = c("Li-Siegmund","MC")))
FHC1000Sp_dat=data.frame(q=log(c(q.vals.HC[seq(110,400,length.out = 25)][1:22],q.vals.HC[seq(110,400,length.out = 25)][1:22]+0.2)),
                         p=-log10(c(LSapproSpFHC1000[16:37],apply(mixgradSplibs_FHC1000,1,mean)[1:22])),
                         sd=c(rep(NA,22),sqrt(apply(mixgradSplibs_FHC1000,1,mse))[1:22]),
                         method=factor(c(rep("Li-Siegmund",22),rep("IS",22)),levels = c("Li-Siegmund","IS")))

FHC1000Lpeval.plot=ggplot(FHC1000Lp_dat)+
  geom_point(aes(x=q,y=p,color=method,shape=method),alpha=0.7,size=5)+
  geom_errorbar(aes(x=q,ymin=p-sd,ymax=p+sd,color=method))+
  scale_colour_manual("",
                      breaks = c("Li-Siegmund", "MC"),
                      values = c("#D2691E","red"),
                      labels=c("Li-Siegmund","Monte Carlo"))+
  scale_shape_manual("",
                     breaks = c("Li-Siegmund", "MC"),
                     values = c(17,18),
                     labels=c("Li-Siegmund","Monte Carlo")) +
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"),
         linetype=guide_legend(title="Method"))+
  ylab("Estimated p-value") +
  xlab("log statistic") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "none",
        legend.text = element_text(face="bold",size = 8),
        legend.title = element_text(face="bold"),
        plot.title = element_text(hjust = 0,size = 30, face = "bold"),
        axis.title = element_text(face="bold",size=20),
        axis.text = element_text(face="bold",size=20),
        text = element_text(size = 20))+
  ggtitle("K=1000")

FHC1000Speval.plot=ggplot(FHC1000Sp_dat)+
  geom_point(aes(x=q,y=p,color=method,shape=method),alpha=0.7,size=5)+
  geom_errorbar(aes(x=q,ymin=p-sd,ymax=p+sd,color=method))+
  scale_colour_manual("",
                      breaks = c("Li-Siegmund","IS"),
                      values = c("#D2691E","blue"),
                      labels=c("Li-Siegmund","IS"))+
  scale_shape_manual("",
                     breaks = c("Li-Siegmund","IS"),
                     values = c(17,19),
                     labels=c("Li-Siegmund","IS")) +
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"),
         linetype=guide_legend(title="Method"))+
  ylab("-log10 Estimated p-value") +
  xlab("log statistic") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "none",
        legend.text = element_text(face="bold",size = 8),
        legend.title = element_text(face="bold"),
        plot.title = element_text(hjust = 0,size = 30, face = "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold",size=20),
        axis.text = element_text(face="bold",size=20))+
  ggtitle("")
FHC1000=ggarrange(FHC1000Lpeval.plot,FHC1000Speval.plot,nrow = 1)

FHC2000Lpeval.plot=ggplot(FHC2000Lp_dat)+
  geom_point(aes(x=q,y=p,color=method,shape=method),alpha=0.7,size=5)+
  geom_errorbar(aes(x=q,ymin=p-sd,ymax=p+sd,color=method))+
  scale_colour_manual("",
                      breaks = c("Li-Siegmund", "MC"),
                      values = c("#D2691E","red"),
                      labels=c("Li-Siegmund","Monte Carlo"))+
  scale_shape_manual("",
                     breaks = c("Li-Siegmund", "MC"),
                     values = c(17,18),
                     labels=c("Li-Siegmund","Monte Carlo")) +
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"),
         linetype=guide_legend(title="Method"))+
  ylab("Estimated p-value") +
  xlab("log statistic") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "none",
        legend.text = element_text(face="bold",size = 8),
        legend.title = element_text(face="bold"),
        plot.title = element_text(hjust = 0,size = 30, face = "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold",size=20),
        axis.text = element_text(face="bold",size=20))+
  ggtitle("K=2000")

FHC2000Speval.plot=ggplot(FHC2000Sp_dat)+
  geom_point(aes(x=q,y=p,color=method,shape=method),alpha=0.7,size=5)+
  geom_errorbar(aes(x=q,ymin=p-sd,ymax=p+sd,color=method))+
  scale_colour_manual("",
                      breaks = c("Li-Siegmund","IS"),
                      values = c("#D2691E","blue"),
                      labels=c("Li-Siegmund","IS"))+
  scale_shape_manual("",
                     breaks = c("Li-Siegmund","IS"),
                     values = c(17,19),
                     labels=c("Li-Siegmund","IS")) +
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"),
         linetype=guide_legend(title="Method"))+
  ylab("-log10 Estimated p-value") +
  xlab("log statistic") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "none",
        legend.text = element_text(face="bold",size = 8),
        legend.title = element_text(face="bold"),
        plot.title = element_text(hjust = 0,size = 30, face = "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold",size=20),
        axis.text = element_text(face="bold",size=20))+
  ggtitle("")

FHC2000=ggarrange(FHC2000Lpeval.plot,FHC2000Speval.plot,nrow = 1)


# q.lg=exp(seq(log(5),log(10),length.out = 100))
# q.sm=exp(seq(log(0.5),log(5),length.out = 100))
# q_lg=exp(seq(log(3.8),log(13),length.out = 150))
q_sm=exp(seq(log(0.5),log(3.7),length.out = 50))

MHCSpevalN107_K300=readRDS("MHCSpevalN107_K300.rds")
MHCSpevalN107_K500=readRDS("MHCSpevalN107_K500.rds")
MHCSpevalN107_K1000=readRDS("MHCSpevalN107_K1000.rds")
MHCSpevalN107_K2000=readRDS("MHCSpevalN107_K2000.rds")

LSapproMHC_K300=readRDS("LSapproMHC_K300.rds")
LSapproMHC_K500=readRDS("LSapproMHC_K500.rds")
LSapproMHC_K1000=readRDS("LSapproMHC_K1000.rds")
LSapproMHC_K2000=readRDS("LSapproMHC_K2000.rds")

MHCLpevalN107.MC_K300=readRDS("MHCLpevalN107.MC_K300.rds")
MHCLpevalN107.MC_K500=readRDS("MHCLpevalN107.MC_K500.rds")
MHCLpevalN107.MC_K1000=readRDS("MHCLpevalN107.MC_K1000.rds")
MHCLpevalN107.MC_K2000=readRDS("MHCLpevalN107.MC_K2000.rds")


mse=function(x) {
  return(mean((log10(x)-log10(mean(x,na.rm = T)))^2,na.rm = T))
}


MHC300Lp_dat=data.frame(q=rep(c(q_sm[seq(17,42,length.out = 10)],exp(seq(log(3),log(13),length=100)[seq(1,100,length.out = 15)][1:5])),2),
                        p=c(apply(MHCLpevalN107.MC_K300,1,mean),LSapproMHC_K300[1:15]),
                        sd=c(sqrt(apply(MHCLpevalN107.MC_K300,1,function(x) sum((x-mean(x))^2))),rep(NA,15)),
                        method=factor(c(rep("MC",15),rep("Li-Siegmund",15)),levels = c("MC","Li-Siegmund")))
MHC300Sp_dat=data.frame(q=rep(exp(seq(log(5),log(11),length.out = 15))[1:12],2),
                        p=-log10(c(apply(MHCSpevalN107_K300,1,mean),LSapproMHC_K300[16:27])),
                        sd=c(sqrt(apply(MHCSpevalN107_K300,1,function(x) mse(x))),rep(NA,12)),
                        method=factor(c(rep("IS",12),rep("Li-Siegmund",12)),levels = c("IS","Li-Siegmund")))

MHC300Sp=ggplot(MHC300Sp_dat)+
  geom_point(aes(x=q,y=p,color=method,shape=method),alpha=0.7,size=5)+
  geom_errorbar(aes(x=q,ymin=p-sd,ymax=p+sd,color=method))+
  scale_colour_manual("",
                      breaks = c("IS",  "Li-Siegmund"),
                      values = c("blue","#D2691E"),
                      labels=c(expression("IS"["CE"]),"Li-Siegmund"))+
  scale_shape_manual("",
                     breaks = c("IS", "Li-Siegmund"),
                     values = c(19,17),
                     labels=c(expression("IS"["CE"]), "Li-Siegmund"))+
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"),
         linetype=guide_legend(title="Method"))+
  ylab("-log10 Estimated P-value") +
  xlab("statistic") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "none",
        legend.text = element_text(face="bold",size = 8),
        legend.title = element_text(face="bold"),
        plot.title = element_text(hjust = 0,size = 30, face = "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold",size=20),
        axis.text = element_text(face="bold",size=20))+
  ggtitle("")

MHC300Lp=ggplot(MHC300Lp_dat)+
  geom_point(aes(x=q,y=p,color=method,shape=method),alpha=0.7,size=5)+
  geom_errorbar(aes(x=q,ymin=p-sd,ymax=p+sd,color=method))+
  scale_colour_manual("",
                      breaks = c("MC", "Li-Siegmund"),
                      values = c("red","#D2691E"),
                      labels=c("Monte Carlo","Li-Siegmund"))+
  scale_shape_manual("",
                     breaks = c("MC", "Li-Siegmund"),
                     values = c(18,17),
                     labels=c("Monte Carlo","Li-Siegmund"))+
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"),
         linetype=guide_legend(title="Method"))+
  ylab("Estimated P-value") +
  xlab("statistic") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "none",
        legend.text = element_text(face="bold",size = 8),
        legend.title = element_text(face="bold"),
        plot.title = element_text(hjust = 0,size = 30, face = "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold",size=20),
        axis.text = element_text(face="bold",size=20))+
  ggtitle("K=300")
MHC300=ggarrange(MHC300Lp,MHC300Sp,nrow = 1)

MHC500Lp_dat=data.frame(q=rep(c(q_sm[seq(17,42,length.out = 10)],exp(seq(log(3),log(13),length=100)[seq(1,100,length.out = 15)][1:5])),2),
                        p=c(apply(MHCLpevalN107.MC_K500,1,mean),LSapproMHC_K500[1:15]),
                        sd=c(sqrt(apply(MHCLpevalN107.MC_K500,1,function(x) sum((x-mean(x))^2))),rep(NA,15)),
                        method=factor(c(rep("MC",15),rep("Li-Siegmund",15)),levels = c("MC","Li-Siegmund")))
MHC500Sp_dat=data.frame(q=rep(exp(seq(log(5),log(11),length.out = 15))[1:12],2),
                        p=-log10(c(apply(MHCSpevalN107_K500,1,mean),LSapproMHC_K500[16:27])),
                        sd=c(sqrt(apply(MHCSpevalN107_K500,1,function(x) mse(x))),rep(NA,12)),
                        method=factor(c(rep("IS",12),rep("Li-Siegmund",12)),levels = c("IS","Li-Siegmund")))

MHC500Sp=ggplot(MHC500Sp_dat)+
  geom_point(aes(x=q,y=p,color=method,shape=method),alpha=0.7,size=5)+
  geom_errorbar(aes(x=q,ymin=p-sd,ymax=p+sd,color=method))+
  scale_colour_manual("",
                      breaks = c("IS",  "Li-Siegmund"),
                      values = c("blue","#D2691E"),
                      labels=c(expression("IS"["CE"]),"Li-Siegmund"))+
  scale_shape_manual("",
                     breaks = c("IS", "Li-Siegmund"),
                     values = c(19,17),
                     labels=c(expression("IS"["CE"]), "Li-Siegmund"))+
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"),
         linetype=guide_legend(title="Method"))+
  ylab("-log10 Estimated P-value") +
  xlab("statistic") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "none",
        legend.text = element_text(face="bold",size = 8),
        legend.title = element_text(face="bold"),
        plot.title = element_text(hjust = 0,size = 30, face = "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold",size=20),
        axis.text = element_text(face="bold",size=20))+
  ggtitle("")

MHC500Lp=ggplot(MHC500Lp_dat)+
  geom_point(aes(x=q,y=p,color=method,shape=method),alpha=0.7,size=5)+
  geom_errorbar(aes(x=q,ymin=p-sd,ymax=p+sd,color=method))+
  scale_colour_manual("",
                      breaks = c( "MC", "Li-Siegmund"),
                      values = c("red","#D2691E"),
                      labels=c("Monte Carlo","Li-Siegmund"))+
  scale_shape_manual("",
                     breaks = c("MC", "Li-Siegmund"),
                     values = c(18,17),
                     labels=c("Monte Carlo","Li-Siegmund"))+
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"),
         linetype=guide_legend(title="Method"))+
  ylab("Estimated P-value") +
  xlab("statistic") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "none",
        legend.text = element_text(face="bold",size = 8),
        legend.title = element_text(face="bold"),
        plot.title = element_text(hjust = 0,size = 30, face = "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold",size=20),
        axis.text = element_text(face="bold",size=20))+
  ggtitle("K=500")
MHC500=ggarrange(MHC500Lp,MHC500Sp,nrow = 1)

MHC1000Lp_dat=data.frame(q=rep(c(q_sm[seq(17,42,length.out = 10)],exp(seq(log(3),log(13),length=100)[seq(1,100,length.out = 15)][1:5])),2),
                        p=c(apply(MHCLpevalN107.MC_K1000,1,mean),LSapproMHC_K1000[1:15]),
                        sd=c(sqrt(apply(MHCLpevalN107.MC_K1000,1,function(x) sum((x-mean(x))^2))),rep(NA,15)),
                        method=factor(c(rep("MC",15),rep("Li-Siegmund",15)),levels = c("MC","Li-Siegmund")))
MHC1000Sp_dat=data.frame(q=rep(exp(seq(log(5),log(11),length.out = 15))[1:12],2),
                        p=-log10(c(apply(MHCSpevalN107_K1000,1,mean),LSapproMHC_K1000[16:27])),
                        sd=c(sqrt(apply(MHCSpevalN107_K1000,1,function(x) mse(x))),rep(NA,12)),
                        method=factor(c(rep("IS",12),rep("Li-Siegmund",12)),levels = c("IS","Li-Siegmund")))

MHC1000Sp=ggplot(MHC1000Sp_dat)+
  geom_point(aes(x=q,y=p,color=method,shape=method),alpha=0.7,size=5)+
  geom_errorbar(aes(x=q,ymin=p-sd,ymax=p+sd,color=method))+
  scale_colour_manual("",
                      breaks = c("IS",  "Li-Siegmund"),
                      values = c("blue","#D2691E"),
                      labels=c(expression("IS"["CE"]),"Li-Siegmund"))+
  scale_shape_manual("",
                     breaks = c("IS", "Li-Siegmund"),
                     values = c(19,17),
                     labels=c(expression("IS"["CE"]), "Li-Siegmund"))+
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"),
         linetype=guide_legend(title="Method"))+
  ylab("-log10 Estimated P-value") +
  xlab("statistic") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "none",
        legend.text = element_text(face="bold",size = 8),
        legend.title = element_text(face="bold"),
        plot.title = element_text(hjust = 0,size = 30, face = "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold",size=20),
        axis.text = element_text(face="bold",size=20))+
  ggtitle("")

MHC1000Lp=ggplot(MHC1000Lp_dat)+
  geom_point(aes(x=q,y=p,color=method,shape=method),alpha=0.7,size=5)+
  geom_errorbar(aes(x=q,ymin=p-sd,ymax=p+sd,color=method))+
  scale_colour_manual("",
                      breaks = c( "MC", "Li-Siegmund"),
                      values = c("red","#D2691E"),
                      labels=c("Monte Carlo","Li-Siegmund"))+
  scale_shape_manual("",
                     breaks = c("MC", "Li-Siegmund"),
                     values = c(18,17),
                     labels=c("Monte Carlo","Li-Siegmund"))+
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"),
         linetype=guide_legend(title="Method"))+
  ylab("Estimated P-value") +
  xlab("statistic") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "none",
        legend.text = element_text(face="bold",size = 8),
        legend.title = element_text(face="bold"),
        plot.title = element_text(hjust = 0,size = 30, face = "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold",size=20),
        axis.text = element_text(face="bold",size=20))+
  ggtitle("K=1000")
MHC1000=ggarrange(MHC1000Lp,MHC1000Sp,nrow = 1)

MHC2000Lp_dat=data.frame(q=rep(c(q_sm[seq(17,42,length.out = 10)],exp(seq(log(3),log(13),length=100)[seq(1,100,length.out = 15)][1:5])),2),
                         p=c(apply(MHCLpevalN107.MC_K2000,1,mean),LSapproMHC_K2000[1:15]),
                         sd=c(sqrt(apply(MHCLpevalN107.MC_K2000,1,function(x) mean((x-mean(x))^2,na.rm = T))),rep(NA,15)),
                         method=factor(c(rep("MC",15),rep("Li-Siegmund",15)),levels = c("MC","Li-Siegmund")))
MHC2000Sp_dat=data.frame(q=rep(exp(seq(log(5),log(11),length.out = 15))[1:12],2),
                         p=-log10(c(apply(MHCSpevalN107_K2000,1,mean),LSapproMHC_K2000[16:27])),
                         sd=c(sqrt(apply(MHCSpevalN107_K2000,1,function(x) mse(x))),rep(NA,12)),
                         method=factor(c(rep("IS",12),rep("Li-Siegmund",12)),levels = c("IS","Li-Siegmund")))

MHC2000Sp=ggplot(MHC2000Sp_dat)+
  geom_point(aes(x=q,y=p,color=method,shape=method),alpha=0.7,size=5)+
  geom_errorbar(aes(x=q,ymin=p-sd,ymax=p+sd,color=method))+
  scale_colour_manual("",
                      breaks = c("IS",  "Li-Siegmund"),
                      values = c("blue","#D2691E"),
                      labels=c(expression("IS"["CE"]),"Li-Siegmund"))+
  scale_shape_manual("",
                     breaks = c("IS", "Li-Siegmund"),
                     values = c(19,17),
                     labels=c(expression("IS"["CE"]), "Li-Siegmund"))+
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"),
         linetype=guide_legend(title="Method"))+
  ylab("-log10 Estimated P-value") +
  xlab("statistic") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "none",
        legend.text = element_text(face="bold",size = 8),
        legend.title = element_text(face="bold"),
        plot.title = element_text(hjust = 0,size = 30, face = "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold",size=20),
        axis.text = element_text(face="bold",size=20))+
  ggtitle("")

MHC2000Lp=ggplot(MHC2000Lp_dat)+
  geom_point(aes(x=q,y=p,color=method,shape=method),alpha=0.7,size=5)+
  geom_errorbar(aes(x=q,ymin=p-sd,ymax=p+sd,color=method))+
  scale_colour_manual("",
                      breaks = c( "MC", "Li-Siegmund"),
                      values = c("red","#D2691E"),
                      labels=c("Monte Carlo","Li-Siegmund"))+
  scale_shape_manual("",
                     breaks = c("MC", "Li-Siegmund"),
                     values = c(18,17),
                     labels=c( "Monte Carlo","Li-Siegmund"))+
  guides(col=guide_legend(title="Method"),
         shape=guide_legend(title="Method"),
         linetype=guide_legend(title="Method"))+
  ylab("Estimated P-value") +
  xlab("statistic") +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.position = "none",
        legend.text = element_text(face="bold",size = 8),
        legend.title = element_text(face="bold"),
        plot.title = element_text(hjust = 0,size = 30, face = "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold",size=20),
        axis.text = element_text(face="bold",size=20))+
  ggtitle("K=2000")
MHC2000=ggarrange(MHC2000Lp,MHC2000Sp,nrow = 1)


plot_combined<-grid.arrange(ggarrange(FHC1000,FHC2000,nrow = 1),
                            ggarrange(MHC300,MHC500,nrow = 1),
                            ggarrange(MHC1000,MHC2000,nrow = 1),
                            layout_matrix = cbind(c(1,2,3)))
ggsave(plot_combined,filename = "FigureA5.png", width = 60, height = 60, units = "cm")


##### Generate Supplementary Figure A6 ######
prep.plots=function(orig.libs,truth=XHK50.lp[4]) {
  #n=length(orig.libs)
  #orig.libs=orig.libs %>% discard(is.null)
  libs=sapply(orig.libs,function(x) unlist(x))
  mean.libs=apply(log10(libs),1,mean)
  bias=log10(libs)-log10(truth)
  mean.bias=apply(bias,1,mean)
  var.bias=apply(bias,1,sd)
  #mean.libs=apply(log10(libs),1,function(x) sum(x)/n)
  mse.libs=apply(libs,1,function(x) mean((log10(x)-mean(log10(x)))^2))
  #mse.libs=apply(libs,1,function(x) sum((log10(x)-sum(log10(x))/n)^2)/n)
  cv.libs=apply(libs,1,function(x) sd(x)/mean(x))
  return(list(mean.libs,mse.libs,cv.libs,mean.bias,var.bias))
}
XHK50.lp=readRDS("XHK50.lp.rds")
XHK50.sp=readRDS("XHK50.sp.rds")
ISK50=list.load("ISK50.RData") #varying sample size q=25
ISK50=prep.plots(ISK50)
ISK50S=list.load("ISK50S.RData") #varying sample size q=25
ISK50S=prep.plots(ISK50S,truth=XHK50.sp[1])
MCK50=list.load("MCK50.RData")
MCK50=prep.plots(MCK50)
MCK50S=list.load("MCK50S.RData")
MCK50S=prep.plots(MCK50S,truth=XHK50.sp[1])
XHK500.lp=unlist(readRDS("XHK500.lp.rds"))
XHK500.sp=readRDS("XH500.sp.rds")
ISK500=list.load("ISK500.RData")
ISK500=prep.plots(ISK500,truth=XHK500.lp[4])
ISK500S=list.load("ISK500S.RData")
ISK500S=prep.plots(ISK500S,truth=XHK500.sp[1])
MCK500=list.load("MCK500.RData")
MCK500=prep.plots(MCK500,truth=XHK500.lp[4])
MCK500S=list.load("MCK500S.RData")
MCK500S=prep.plots(MCK500S,truth=XHK500.sp[1])
XHK5000.lp=readRDS("XHK5000.lp.rds")
ISK5000=list.load("ISK5000.RData") #varying sample size q=25
ISK5000=prep.plots(ISK5000,truth=XHK5000.lp[4])
MCK5000=list.load("MCK5000.RData")
MCK5000=prep.plots(MCK5000,truth=XHK5000.lp[4])#varying sample size q=25
ISK5000S=list.load("ISK5000S.RData")
ISK5000S=prep.plots(ISK5000S,truth=XHK5000.sp[1])
MCK5000S=list.load("MCK5000S.RData")
MCK5000S=prep.plots(MCK5000S,truth=XHK5000.sp[1])

XH.time=list.load("XH.time.RData")
XH.time=log(sapply(XH.time, function(x) as.numeric(x,units="secs")))
ISK50.time=list.load("ISK50.time.RData")
ISK50.time=log(sapply(ISK50.time, function(x) as.numeric(x,units="secs")))
MCK50.time=list.load("MCK50.time.RData")
MCK50.time=log(sapply(MCK50.time, function(x) as.numeric(x,units="secs")))
ISK500.time=list.load("ISK500.time.RData")
ISK500.time=log(sapply(ISK500.time, function(x) as.numeric(x,units="secs")))
MCK500.time=list.load("MCK500.time.RData")
MCK500.time=log(sapply(MCK500.time, function(x) as.numeric(x,units="secs")))
ISK5000.time=list.load("ISK5000.time.RData")
ISK5000.time=log(sapply(ISK5000.time, function(x) as.numeric(x,units="secs")))
MCK5000.time=list.load("MCK5000.time.RData")
MCK5000.time=log(sapply(MCK5000.time, function(x) as.numeric(x,units="secs")))
XHS.time=list.load("XHS.time.RData")
XHS.time=log(sapply(XHS.time, function(x) as.numeric(x,units="secs")))
ISK50S.time=list.load("ISK50S.time.RData")
ISK50S.time=log(sapply(ISK50S.time, function(x) as.numeric(x,units="secs")))
MCK50S.time=list.load("MCK50S.time.RData")
MCK50S.time=log(sapply(MCK50S.time, function(x) as.numeric(x,units="secs")))
ISK500S.time=list.load("ISK500S.time.RData")
ISK500S.time=log(sapply(ISK500S.time, function(x) as.numeric(x,units="secs")))
MCK500S.time=list.load("MCK500S.time.RData")
MCK500S.time=log(sapply(MCK500S.time, function(x) as.numeric(x,units="secs")))
ISK5000S.time=list.load("ISK5000S.time.RData")
ISK5000S.time=log(sapply(ISK5000S.time, function(x) as.numeric(x,units="secs")))
MCK5000S.time=list.load("MCK5000S.time.RData")
MCK5000S.time=log(sapply(MCK5000S.time, function(x) as.numeric(x,units="secs")))

data1<- data.frame(value=c(-log10(rep(XHK50.lp[4],4)),-ISK50[[1]][1:4],-MCK50[[1]][1:4],
                           -log10(rep(XHK500.lp[4],4)),-ISK500[[1]][1:4],-MCK500[[1]][1:4],
                           -log10(rep(XHK5000.lp[4],4)),-ISK5000[[1]][1:4],-MCK5000[[1]][1:4]),
                   mse=sqrt(c(rep(NA,4),ISK50[[2]][1:4],MCK50[[2]][1:4],
                              rep(NA,4),ISK500[[2]][1:4],MCK500[[2]][1:4],
                              rep(NA,4),ISK5000[[2]][1:4],MCK5000[[2]][1:4])),
                   bias=c(rep(0,4),ISK50[[4]][1:4],MCK50[[4]][1:4],
                          rep(0,4),ISK500[[4]][1:4],MCK500[[4]][1:4],
                          rep(0,4),ISK5000[[4]][1:4],MCK5000[[4]][1:4]),
                   sd=c(rep(NA,4),ISK50[[5]][1:4],MCK50[[5]][1:4],
                        rep(NA,4),ISK500[[5]][1:4],MCK500[[5]][1:4],
                        rep(NA,4),ISK5000[[5]][1:4],MCK5000[[5]][1:4]),
                   time=c(rep(XH.time[1],4),ISK50.time,MCK50.time,
                          rep(XH.time[2],4),ISK500.time,MCK500.time,
                          rep(XH.time[3],4),ISK5000.time,MCK5000.time),
                   samplesize=rep(c(3:6),9),
                   #samplesize=factor(rep(c("1e3","1e4","1e5","1e6"),9),levels =c("1e3","1e4","1e5","1e6")),
                   #p=rep(rep(c("small","large"),each=5),6),
                   #K=rep(c(50,500,5000),each=12),
                   K=factor(rep(c(50,500,5000),each=12),levels =c(50,500,5000)),
                   Method=factor(rep(rep(c("Barnett-Lin","IS","Monte Carlo"),each=4),3),levels = c("Barnett-Lin","IS","Monte Carlo")))

data2<- data.frame(value=c(-log10(rep(XHK50.sp[4],4)),-ISK50S[[1]][1:4],-MCK50S[[1]][1:4],
                           -log10(rep(XHK500.sp[4],4)),-ISK500S[[1]][1:4],-MCK500S[[1]][1:4],
                           -log10(rep(XHK5000.sp[4],4)),-ISK5000S[[1]][1:4],-MCK5000S[[1]][1:4]),
                   mse=sqrt(c(rep(NA,4),ISK50S[[2]][1:4],MCK50S[[2]][1:4],
                              rep(NA,4),ISK500S[[2]][1:4],MCK500S[[2]][1:4],
                              rep(NA,4),ISK5000S[[2]][1:4],MCK5000S[[2]][1:4])),
                   bias=c(rep(0,4),ISK50S[[4]][1:4],MCK50S[[4]][1:4],
                          rep(0,4),ISK500S[[4]][1:4],MCK500S[[4]][1:4],
                          rep(0,4),ISK5000S[[4]][1:4],MCK5000S[[4]][1:4]),
                   sd=c(rep(NA,4),ISK50S[[5]][1:4],MCK50S[[5]][1:4],
                        rep(NA,4),ISK500S[[5]][1:4],MCK500S[[5]][1:4],
                        rep(NA,4),ISK5000S[[5]][1:4],MCK5000S[[5]][1:4]),
                   time=c(rep(XHS.time[1],4),ISK50S.time,MCK50S.time,
                          rep(XHS.time[2],4),ISK500S.time,MCK500S.time,
                          rep(XHS.time[3],4),ISK5000S.time,MCK5000S.time),
                   samplesize=rep(c(3:6),9),
                   K=factor(rep(c(50,500,5000),each=12),levels =c(50,500,5000)),
                   Method=factor(rep(rep(c("Barnett-Lin","IS","Monte Carlo"),each=4),3),levels = c("Barnett-Lin","IS","Monte Carlo")))

data<-data1[which(data1$K==50),]
p1<-ggplot(data)+aes(x=samplesize,y=bias,color=Method,shape=Method)+geom_point(size=3.5, position = position_dodge(0.3))+
  geom_errorbar(aes(ymin=bias-sd, ymax=bias+sd), size=0.6,width =0.5,lwd=0.3,position=position_dodge(0.3))+
  theme_bw() +
  theme(axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=20,face="bold"),
        strip.text.x = element_text(size = 20),legend.position = 'none',
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
  scale_x_continuous(breaks=c(3:6))+xlab("log10 sampling size")+ylab("Bias at log10 scale")

data<-data1[which(data1$K==500),]
p2<-ggplot(data)+aes(x=samplesize,y=bias,color=Method,shape=Method)+geom_point(size=3.5, position = position_dodge(0.3))+
  geom_errorbar(aes(ymin=bias-sd, ymax=bias+sd), size=0.6,width =0.5,lwd=0.3,position=position_dodge(0.3))+
  theme_bw() +
  theme(axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=20,face="bold"),
        strip.text.x = element_text(size = 20),legend.position = 'none',
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
  scale_x_continuous(breaks=c(3:6))+xlab("log10 sampling size")+ylab("Bias at log10 scale")


data<-data1[which(data1$K==5000),]
p3<-ggplot(data)+aes(x=samplesize,y=bias,color=Method,shape=Method)+geom_point(size=3.5, position = position_dodge(0.3))+
  geom_errorbar(aes(ymin=bias-sd, ymax=bias+sd), size=0.6,width =0.5,lwd=0.3,position=position_dodge(0.3))+
  theme_bw() +
  theme(axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=20,face="bold"),
        strip.text.x = element_text(size = 20),legend.position = 'none',
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
  scale_x_continuous(breaks=c(3:6))+xlab("log10 sampling size")+ylab("Bias at log10 scale")

p4<-ggplot(data1)+
  geom_point(aes(x=K,y=time,color=Method,shape=Method),size=3.5)+
  geom_line(aes(x=K,y=time,color=Method,linetype=as.character(samplesize-2),group=interaction(Method,samplesize)))+
  theme_bw() +
  theme(axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=20,face="bold"),
        strip.text.x = element_text(size = 20),legend.position = 'none',
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
  scale_x_discrete(breaks=c(50,500,5000))+xlab("Dimension K")+ylab("log time in second")

data<-data2[which(data1$K==50),]
p5<-ggplot(data)+aes(x=samplesize,y=bias,color=Method,shape=Method)+geom_point(size=3.5, position = position_dodge(0.3))+
  geom_errorbar(aes(ymin=bias-sd, ymax=bias+sd), size=0.6,width =0.5,lwd=0.3,position=position_dodge(0.3))+
  theme_bw() +
  theme(axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=20,face="bold"),
        strip.text.x = element_text(size = 20),legend.position = 'none',
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
  scale_x_continuous(breaks=c(3:6))+xlab("log10 sampling size")+ylab("Bias at log10 scale")

data<-data2[which(data1$K==500),]
p6<-ggplot(data)+aes(x=samplesize,y=bias,color=Method,shape=Method)+geom_point(size=3.5, position = position_dodge(0.3))+
  geom_errorbar(aes(ymin=bias-sd, ymax=bias+sd), size=0.6,width =0.5,lwd=0.3,position=position_dodge(0.3))+
  theme_bw() +
  theme(axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=20,face="bold"),
        strip.text.x = element_text(size = 20),legend.position = 'none',
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
  scale_x_continuous(breaks=c(3:6))+xlab("log10 sampling size")+ylab("Bias at log10 scale")

data<-data2[which(data1$K==5000),]
p7<-ggplot(data)+aes(x=samplesize,y=bias,color=Method,shape=Method)+geom_point(size=3.5, position = position_dodge(0.3))+
  geom_errorbar(aes(ymin=bias-sd, ymax=bias+sd), size=0.6,width =0.5,lwd=0.3,position=position_dodge(0.3))+
  theme_bw() +
  theme(axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=20,face="bold"),
        strip.text.x = element_text(size = 20),legend.position = 'none',
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
  scale_x_continuous(breaks=c(3:6))+xlab("log10 sampling size")+ylab("Bias at log10 scale")

p8<-ggplot(data2)+
  geom_point(aes(x=K,y=time,color=Method,shape=Method),size=3.5)+
  geom_line(aes(x=K,y=time,color=Method,linetype=as.character(samplesize),group=interaction(Method,samplesize)))+
  theme_bw() +
  theme(axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=20,face="bold"),
        strip.text.x = element_text(size = 20),legend.position = 'none',
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
  scale_x_discrete(breaks=c(50,500,5000))+xlab("Dimension K")+ylab("log time in second")

library(gridExtra)
library(ggpubr)

plot1<-ggarrange(p1,p2,p3,p4, ncol = 4)
plot2<-ggarrange(p5,p6,p7,p8,ncol = 4)

plot_combined<-grid.arrange(plot1,plot2,
                            layout_matrix = cbind(c(1,2)))
ggsave(plot_combined,filename = "revisionFigure1.png", width = 60, height = 40, units = "cm")



##### Generate Supplementary Figure A7 ######
library(gridExtra)
library(ggpubr)
MCp=readRDS("US.MCp.rds")
XHp=readRDS("US.XHp.rds")
ISp=readRDS("US.ISp.rds")
MSTp=readRDS("US.MSTp.rds")
HCpvals=readRDS("US.HCpvals.rds")
US.data=data.frame(bias=c(log10(MCp)-log10(XHp),log10(ISp)-log10(XHp),log10(MSTp)-log10(XHp),log10(HCpvals)-log10(XHp)),
                   Method=factor(c(rep("Monte Carlo",length(XHp)),rep("IS",length(XHp)),rep("MST",length(XHp)),rep("UFI",length(XHp))),levels = c("MST","IS","Monte Carlo","UFI")))

US.plot=ggplot(US.data)+
  geom_boxplot(aes(x=Method,y=bias,fill=Method))+
  labs(y="Bias at log10 scale", x= "")+
  geom_hline(aes(yintercept=0),color="black",linetype="dashed")+
  theme_bw()+
  theme(legend.text.align = 0,
        legend.position ="bottom",
        legend.text=element_text(size=20,face = "bold"),
        legend.background = element_rect(fill='transparent'),
        legend.title = element_text(face="bold",size=20),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold",size=20),
        axis.text = element_text(face="bold",size=20))
ggsave(US.plot,filename = "revisionFigure2.png", width = 40, height = 40, units = "cm")

