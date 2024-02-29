
##laptop
#setwd("C:/Users/richa/Dropbox/Harvard/Research/Gary polarization manuscript work/Australia_2016_behavioral data")
setwd("C:/Users/Astaroph/Dropbox/Harvard/Research/Gary polarization manuscript work/Australia_2016_behavioral data")
area=read.table("R_ready_jalmenus_australia_polarized_trials_v2.csv",header=T, sep=",")
str(area)

library(nlme)
library(lme4)
library(car)
library(effects)
library(lsmeans)
library(pbkrtest)
str(area)

plot(Pol_diff~Treatment,data=area)
fit=aov(Pol_diff~Treatment,data=area)
summary(fit)

boxplot(Pol_diff~Treatment,data=area,notch=TRUE)
str(area)

blue=subset(area,Treatment=="3M+Blue R374")
yellow=subset(area,Treatment=="3M+Yellow R312")
red=subset(area,Treatment=="3M+Red G280")
gary=subset(area,Treatment=="HNPB")
red
t.test(gary$Pol_diff)
t.test(red$Pol_diff)
t.test(yellow$Pol_diff)
t.test(blue$Pol_diff)

library(ggplot2)

library(ggplot2)

violin_1=ggplot(
  data = area
  , mapping = aes(
    y =Pol_diff
    , x =factor(Treatment)
    , colour=factor(Treatment)
  )
)+
  labs(
    x = 'Color Treatment'
    ,title= 'Polarized approaches - Depolarized approaches'
    , y = 'Pol-Depol'
  )+
  geom_jitter(width=.1,
              mapping = aes(
                colour =factor(Treatment)
              )
              , size=2.5
  )+
  geom_violin(alpha=.4, draw_quantiles = c(0.25, 0.5, 0.75),
              mapping=aes(fill=factor(Treatment),colour=factor(Treatment)))+
  geom_hline(linetype='dashed',yintercept=0)


violin_1
#Changing the colors of the treatments
violin_2=violin_1+scale_color_manual(values=c("#00B2BA","#DD190D","#F8D500", "#7D7D7D"))+
  theme_bw()+scale_fill_manual(values=c("#00B2BA","#DD190D","#F8D500", "#7D7D7D"))
violin_2
##Resizing and moving the text in the graph
violin_3=violin_2+theme(legend.title = element_text(colour="Black", size=16, face="bold"))+
  theme(legend.text = element_text(colour="Black", size = 14, face = "bold"))+
  theme(legend.title = element_blank())+
  theme(axis.text=element_text(colour="Black", size=16, face="plain"))+
  theme(axis.title=element_text(colour="Black", size=16, face="bold"))+
  theme(plot.title=element_text(colour="Black", size=16, face="bold"))+
  theme(axis.title.y=element_text(vjust=1))+
  theme(strip.text=element_text(colour="Black", size = 16, face = "plain"))
violin_3
library(grid)
violin_4=violin_3+theme(panel.spacing = unit(0.5, "lines"),
                        legend.position="none")
violin_4
#remove gridlines and top/right borders and facets
violin_5=violin_4+ theme(axis.line = element_line(colour = "black"),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.border=element_blank())+
  theme(axis.line = element_line(color = 'black'))+
  #theme(strip.text=element_blank())+
  theme(strip.background=element_blank())
violin_5
  

