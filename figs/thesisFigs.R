setwd("/Users/yongsheng/Documents/PhD/TeX/PPP2R5C CYS PhD/figs")
require(grid) # for unit function in ggplot2
library(Cairo) # to print special character like † in ggplot2, and save to pdf
library(reshape2)
library(Hmisc)
library(gridExtra)

# define thesis theme based on theme_bw
theme_bw
# x axis text is italic
theme_thesis_it <- function(base_size=12, base_family="") {
    theme_bw(base_size = base_size, base_family = base_family) %+replace% 
        theme(axis.text=element_text(size=10),
              axis.text.x = element_text(angle = 30, vjust=0.9,hjust=0.8),
              axis.title=element_text(size=10),
              legend.position="right",
              legend.justification=c(3,0),
              legend.text = element_text(size = 7),
              legend.background = element_rect(fill="transparent", colour="transparent"),
              legend.key = element_blank(),
              legend.key.size = unit(0.3, "cm"),
              strip.background=element_rect(fill="transparent", size=0, colour="white"),
              axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank())
}
# x axis text is normal
theme_thesis <- function(base_size=12, base_family="") {
    theme_bw(base_size = base_size, base_family = base_family) %+replace% 
        theme(axis.text=element_text(size=10),
              axis.title=element_text(size=10),
              legend.position="right",
              legend.justification=c(3,0),
              legend.text = element_text(size = 7),
              legend.background = element_rect(fill="transparent", colour="transparent"),
              legend.key = element_blank(),
              legend.key.size = unit(0.3, "cm"),
              strip.background=element_rect(fill="transparent", size=0, colour="white"),
              axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank())
}


# read in overview data for in vivo miR knockdown
data.all <- read.xlsx("~/Data/20140121-0312 AAV full experiment/140121-0312-AAV-full-exp-overview with 31 mice.xlsx", sheetName="short overview", rowIndex=7:38, colIndex=1:33,header=TRUE)
data.all$trt <- factor(data.all$trt, levels=c("Random", "Fasting", "Refed"))
data.all$aav <- factor(data.all$aav, levels=c("miRNC", "miR12"))

str(data.all)

# Figure 2-1
###################################################
# fig2-1 import liver qPCR data and plot
###################################################

liverPPP2R5C <- read.xlsx("~/Data/20101115 qPCR analysis of mouse liver with PPP2R5C/liver.xls", sheetName="liver", colIndex=c(1:3, 4,6,8,10,12), rowIndex=147:171,header=TRUE)
liverPPP2R5C$mouse <- factor(liverPPP2R5C$mouse,levels=c("WT","db/db"))
liverPPP2R5C$trt <- factor(liverPPP2R5C$trt, levels=c("Random", "Fasting", "Refed"))
library(reshape2)
liverPPP2R5C.m <- melt(liverPPP2R5C[,-8], id=1:3)

# test for raw qPCR data liver
str(liverPPP2R5C.m)
liverPPP2R5C.m$test <- apply(liverPPP2R5C.m[,1:2], 1, paste, collapse="-")
liverPPP2R5C.m <- liverPPP2R5C.m[liverPPP2R5C.m$value < 10,]
liverPPP2R5C.m$test <- factor(liverPPP2R5C.m$test)
liver.list <- split(liverPPP2R5C.m, liverPPP2R5C.m$variable)
liver.t <- lapply(liver.list, function(x) pairwise.t.test(x$value, x$test, p.adjust="BH"))
liver.t

liver.aov.list <- lapply(liver.list, function(x) TukeyHSD(aov(value ~ mouse + trt, data=x)))
liver.aov.list
summary(liverPPP2R5C.m)
liverPPP2R5C.m$variable <- factor(liverPPP2R5C.m$variable, levels=c("variant1.4", "variant2", "variant3", "variant4"), labels=c("Variant 1", "Variant 2", "Variant 3", "Variant 4"))
df = data.frame(x=c(1.167,2, 1, 1.167,2,1.5,1.5,1.5), y=c(1.55,1.6,2.1,1.5,1.8,2.3,2.3,2.2), variable=c("Variant 1","Variant 1","Variant 4","Variant 4","Variant 4","Variant 1","Variant 2", "Variant 4"), trt=c("Refed","Fasting","Fasting", "Refed","Fasting","Fasting", "Fasting", "Fasting"),label=c("*","***","***", "*","*", "†††","†††","†"))
p <- ggplot(liverPPP2R5C.m, aes(x=mouse, y=value, fill=trt))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5), width=0.2) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=trt)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(x="", y="PPP2R5C mRNA Levels\n(TBP normalized)") + theme_bw() + theme(axis.text=element_text(size=10),axis.title=element_text(size=10), legend.position="right", legend.justification=c(3,0), legend.text = element_text(size = 7),legend.background = element_rect(fill="transparent"), legend.key = element_blank(), legend.key.size = unit(0.3, "cm"), strip.background=element_rect(fill="transparent", size=0, colour="white")) + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + scale_y_continuous(expand = c(0, 0), limits=c(0,2.5)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + facet_wrap(~ variable) + geom_text(data=df, mapping=aes(label=label,x,y,trt),size=4)


# this shorter version also works the same with theme defined as theme_thesis, and replace chuncky theme code in the following.
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5), width=0.2) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=trt)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(x="", y="PPP2R5C mRNA Levels\n(TBP normalized)") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,2.5)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + facet_wrap(~ variable) + geom_text(data=df, mapping=aes(label=label,x,y,trt),size=4)
ggsave("fig2-1 liver ppp2r5c.pdf",width=16,height=12,units="cm", device=cairo_pdf)



# Figure 2-2
###################################################
# fig 2-2 import adipose qPCR data and plot
###################################################

adiposePPP2R5C <- read.xlsx("~/Data/20101221 qPCR of mouse adipose tissue/20101222mouse adipose ppp2r5c .xls", sheetName="2006.12.20.csv", colIndex=c(1:3, 4,6,8,10,12), rowIndex=147:171,header=TRUE)
adiposePPP2R5C$mouse <- factor(adiposePPP2R5C$mouse,levels=c("WT","db/db"))
adiposePPP2R5C$trt <- factor(adiposePPP2R5C$trt, levels=c("Random", "Fasting", "Refed"))
library(reshape2)
adiposePPP2R5C.m <- melt(adiposePPP2R5C[,-8], id=1:3)

# test for raw qPCR data liver
str(adiposePPP2R5C.m)
summary(adiposePPP2R5C.m)
adiposePPP2R5C.m$test <- apply(adiposePPP2R5C.m[,1:2], 1, paste, collapse="-")
adipose.list <- split(adiposePPP2R5C.m, adiposePPP2R5C.m$variable)
adipose.t <- lapply(adipose.list, function(x) pairwise.t.test(x$value, x$test, p.adjust="BH"))
adipose.t

adipose.aov.list <- lapply(adipose.list, function(x) TukeyHSD(aov(value ~ mouse + trt, data=x)))
adipose.aov.list

adiposePPP2R5C.m$variable <- factor(adiposePPP2R5C.m$variable, levels=c("Variant1.4", "Variant2", "Variant3", "Variant4"), labels=c("Variant 1", "Variant 2", "Variant 3", "Variant 4"))
df = data.frame(x=c(1.5), y=c(9), variable=c("Variant 4"), trt=c( "Fasting"),label=c("†"))
p <- ggplot(adiposePPP2R5C.m, aes(x=mouse, y=value, fill=trt))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5), width=0.2) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=trt)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(x="", y="PPP2R5C mRNA Levels\n(TBP normalized)") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,10.5)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + facet_wrap(~ variable) + geom_text(data=df, mapping=aes(label=label,x,y,trt),size=4)
ggsave("fig2-2 adipose ppp2r5c.pdf",width=16,height=12,units="cm", device=cairo_pdf)


###################################################
# fig 2-3 import muscle qPCR data and plot
###################################################
musclePPP2R5C <- read.xlsx("~/Data/20101122 qPCR of mouse muscle RT from Katrin/muscle ppp2r5c.xls", sheetName="2006.12.20.csv", colIndex=c(1:3, 4,6,8,10,12), rowIndex=147:163,header=TRUE)
musclePPP2R5C$mouse <- factor(musclePPP2R5C$mouse,levels=c("WT","db/db"))
musclePPP2R5C$trt <- factor(musclePPP2R5C$trt, levels=c("Fasting", "Refed"))

musclePPP2R5C.m <- melt(musclePPP2R5C[,-8], id=1:3)

# test for raw qPCR data liver
summary(musclePPP2R5C.m)
musclePPP2R5C.m$test <- apply(musclePPP2R5C.m[,1:2], 1, paste, collapse="-")
muscle.list <- split(musclePPP2R5C.m, musclePPP2R5C.m$variable)
muscle.t <- lapply(muscle.list, function(x) pairwise.t.test(x$value, x$test, p.adjust="BH"))
muscle.t

muscle.aov.list <- lapply(muscle.list, function(x) TukeyHSD(aov(value ~ mouse + trt, data=x)))
muscle.aov.list

musclePPP2R5C.m$variable <- factor(musclePPP2R5C.m$variable, levels=c("variant1.4", "variant2", "variant3", "variant4"), labels=c("Variant 1", "Variant 2", "Variant 3", "Variant 4"))
df = data.frame(x=c(1.5,1.125,1.125,2.125,1.125), y=c(4,4.1,3,4.1, 4), variable=c("Variant 3", "Variant 4", "Variant 3","Variant 3", "Variant 2"), trt=c( "Fasting", "Refed","Refed","Refed", "Refed"),label=c("†††", "*", "***","***","*"))
p <- ggplot(musclePPP2R5C.m, aes(x=mouse, y=value, fill=trt))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5), width=0.2) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=trt)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(x="", y="PPP2R5C mRNA Levels\n(TBP normalized)") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,8)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + facet_wrap(~ variable) + geom_text(data=df, mapping=aes(label=label,x,y,trt),size=4)
ggsave("fig2-3 muscle ppp2r5c.pdf",width=16,height=12,units="cm", device=cairo_pdf)



###################################################
# fig 2-4 import liver hepatosteatosis qPCR data and plot
###################################################
liver_34sample <- read.xlsx("~/Data/20120427 qPCR with 34 sample from mauricio/20120430 qPCR with 34 sample from mauricio TBP norm.xls", sheetName="2006.12.20.csv", colIndex=c(1:3, 5), rowIndex=161:195,header=TRUE)

liver_34sample$fill <- factor(rep(c("Control","Treatment","Control","Treatment","Control","Treatment","Control","Treatment"), c(3,4,5,6,2,2,6,6)), levels=c("Control","Treatment"))
liver_34sample$trt <- factor(liver_34sample$trt, levels=c("Control Diet with PBS", "MCD Diet with PBS", "Control Diet with miRNA", "MCD Diet with miRNA", "STZ Liver Control", "STZ Liver", "Dex Liver Control", "Dex Liver"))

liver_34sample$test <- factor(rep(c("Dex","miRNA","PBS","STZ"), c(7,11,4,12)), levels=c("PBS","miRNA", "STZ", "Dex"))


liver34.list <- split(liver_34sample, liver_34sample$model)
liver34.list2 <- split(liver_34sample, liver_34sample$test)
liver34.t <- lapply(liver34.list, function(x) pairwise.t.test(x$PPP2R5C, x$trt, p.adjust="BH"))
liver34.t2 <- lapply(liver34.list2, function(x) pairwise.t.test(x$PPP2R5C, x$trt, p.adjust="BH"))


liver34.t
liver34.t2 # PBS 0.008 p-value
mcd <- liver34.list[[1]]
mcd$diet <- rep(c("control", "mcd", "control", "mcd"), c(5,6,2,2))
mcd$aav <- rep(c("miRNA", "PBS"), c(11,4))

mcd.aov <- TukeyHSD(aov(PPP2R5C ~ aav + diet, data=mcd))
mcd.aov

df = data.frame(x=2, y=0.65, model="Hepatosteatosis",trt="MCD Diet with PBS", fill="Treatment",label="**")

p <- ggplot(liver_34sample, aes(x=trt, y=PPP2R5C, fill=fill))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5), width=0.2) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=trt)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(x="", y="PPP2R5C mRNA Levels\n(TBP normalized)") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,1.5)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + facet_grid(.~ model, scales="free_x",space="free_x")+ geom_text(data=df, mapping=aes(label=label,x,y,trt),size=4) 
ggsave("fig2-4 liver 34sample ppp2r5c.pdf",width=16,height=10,units="cm", device=cairo_pdf)


############################
# fig 2-5 kd efficiency of shRNAs
###############################

luc.kd <- read.xlsx("~/Data/20110216 dose effect of shR miR on RL-Var2 with 0.02 to 2 of shR:RL/Dose effect of shR miR.xlsx", sheetName="Sheet1", colIndex=12:17, rowIndex=14:28,header=TRUE)

luc.kd$species <- factor(rep(c("shRNA", "miRNA","shRNA", "miRNA"), c(2,2,4,6)), levels=c("shRNA", "miRNA"))
luc.kd$trt <- factor(rep(c("Control", "Knockdown"), c(4,10)), levels=c("Control", "Knockdown"))
luc.kd$fill <- factor(rep(c(1,3,4,2,3,4),c(4,2,2,2,2,2)))
luckd.m <- melt(luc.kd, id=c(1,7,8,9))
luckd.m$kd <- factor(luckd.m$kd, levels=c("shRNC", "shR2", "shR3", "miRNC", "miR1", "miR2", "miR3"))
luckd.m$variable <- factor(luckd.m$variable, levels=c("r0.02","r0.1","r0.5","r1","r2"), labels=c("0.02", "0.1", "0.5", "1", "2"))

p <- ggplot(luckd.m, aes(x=variable, y=value, fill=fill,group=kd))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5), width=0.2) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=kd)) + scale_fill_grey(name="", start=0.8, end=0, labels=c("shRNC/miRNC", "miR1", "shR2/miR2","shR3/miR3")) + labs(x="shRNA or miRNA to Renilla Luc Ratio", y="Relative RL/FL") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + facet_grid(.~ species, scales="free_x",space="free_x") 
ggsave("fig2-5 shR miR kd luc.pdf",width=16,height=8,units="cm", device=cairo_pdf)



############################
# fig 2-6 kd efficiency of miRNAs
###############################
v2kd <- read.xlsx("~/Data/20120325 testing miR1-12 in 293 /kd efficiency.xlsx", sheetName="Sheet1", colIndex=c(1,6:8), rowIndex=5:18,header=TRUE)

v2kd$mir <- factor(v2kd$mir, levels=c("miRNC", "miR1", "miR2", "miR3", "miR4", "miR5", "miR6", "miR7", "miR8", "miR9", "miR10", "miR11","miR12"))

v2kd.m <- melt(v2kd, id=1)
v2kd.m$variable <- factor(v2kd.m$variable, levels=c("ratio.0.5.1", "ratio.1.1", "ratio.2.1"), labels=c("0.5", "1", "2"))

p <- ggplot(v2kd.m, aes(x=mir, y=value, fill=variable,group=variable))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5), width=0.2) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=variable)) + scale_fill_grey(name="", start=0.8, end=0, labels=c("miR/V2 = 0.5", "miR/V2 = 1", "miR/V2 = 2")) + labs(x="miRNAs Knockdown on Variant 2 of PPP2R5C", y="Relative anti-HA (miR1-12 / miRNC)") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5))) 
ggsave("fig2-6b miR kd v2.pdf",width=16,height=8,units="cm", device=cairo_pdf)

###################
#fig 2-9 ad-shR3 kd in hepa with different MOI
##################
hepakd <- read.xlsx("~/Data/20110530 qPCR analysis of Hepa 1c1 infected with AD-NC and AD-PPP2R5C/20110530 AD infected Hepa 1c1 ppp2r5c knockdown.xls", sheetName="2006.12.20.csv", colIndex=c(1:4, 6,8), rowIndex=141:159,header=TRUE)
str(hepakd.m)
hepakd.m <- melt(hepakd, id=c(1:3))
hepakd.m$variable <- factor(hepakd.m$variable, levels=c("variant2","variant3","variant4"), labels=c("Variant 2", "Variant 3","Variant 4"))
hepakd.m$shR <- factor(hepakd.m$shR, levels=c("shRNC", "shR3"))
hepakd.m$moi <- factor(hepakd.m$moi, levels=c(10,100,200))
p <- ggplot(hepakd.m, aes(x=moi, y=value, fill=shR))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5), width=0.2) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=shR)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(x="MOI of Adenovirus Infection on Hepa 1-6", y="PPP2R5C mRNA Levels\n(TBP normalized)") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,1.3), breaks=c(0,0.5,1)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + facet_wrap(~ variable)
ggsave("fig2-9 hepa ad kd.pdf",width=16,height=8,units="cm", device=cairo_pdf)



###################
#fig 2-10 glucose and lactate production upon kd
##################

glucose <- read.xlsx("~/Data/20110719 glucose level in hepa 1-6 after ad-ppp2r5c/20110719 glucose.xls", sheetName="Sheet0", colIndex=c(1,5), rowIndex=58:70,header=TRUE)


lactate <- read.xlsx("~/Data/20110725 lactate production of Hepa 1-6 with ad-ppp2r5c/lactate/20110725 lactate A1.xls", sheetName="Sheet0", colIndex=c(4), rowIndex=73:85,header=TRUE)

glulac <- cbind(glucose, lactate)
glulac$shR <- factor(rep(c("shRNC", "shR3", "shRNC", "shR3"),c(3,3,3,3)), levels=c("shRNC", "shR3"))
glulac$dos <- factor(rep(c(100,200), c(6,6)))
glulac.m <- melt(glulac, id=c(2,4,5))
glulac.m$variable <- factor(glulac.m$variable, labels=c("Glucose Consumption", "Lactate Production"))

glulac.list <- split(glulac.m, glulac.m$variable)
glulac.t <- lapply(glulac.list, function(x) pairwise.t.test(x$value, x$trt, p.adjust="BH"))
glulac.t

glulac.aov <- lapply(glulac.list, function(x) TukeyHSD(aov(value ~ dos + shR, data=x)))
glulac.aov

df = data.frame(x=c(1.125,2.125), y=c(1.45,1.45), shR=c("shR3", "shR3"),label=c("**", "**"))
p <- ggplot(glulac.list[[1]], aes(x=dos, y=value, fill=shR,group=shR))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5), width=0.2) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=shR)) + scale_fill_grey(name="", start=0.8, end=0, labels=c("shRNC", "shR3")) + labs(x="Dosage of AV-shR (MOI)", y="Relative Glucose Consumption") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,1.8)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + geom_text(data=df, mapping=aes(label=label,x,y,shR),size=4)
ggsave("fig2-10a glucose kd.pdf",width=8,height=6,units="cm", device=cairo_pdf)


df = data.frame(x=c(1.125,2.125), y=c(1.45,1.47), shR=c("shR3", "shR3"),label=c("**", "**"))
p <- ggplot(glulac.list[[2]], aes(x=dos, y=value, fill=shR,group=shR))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5), width=0.2) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=shR)) + scale_fill_grey(name="", start=0.8, end=0, labels=c("shRNC", "shR3")) + labs(x="Dosage of AV-shR (MOI)", y="Relative Lactate Production") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,1.8)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + geom_text(data=df, mapping=aes(label=label,x,y,shR),size=4)
ggsave("fig2-10b lactate kd.pdf",width=8,height=6,units="cm", device=cairo_pdf)


#########################################
# fig2-11 seahorse oligomycin 1uM time point 10,11,14 **,*, **
#########################################
seahorse <- read.xlsx("~/Documents/Lab Meetings/Draft/FigData.xlsx", sheetName="Sheet1", rowIndex=133:147, colIndex=1:22,header=TRUE)
str(seahorse)

seahorse.m <- melt(seahorse, id=22)
colnames(seahorse.m) <- c("time","well","rate")

seahorse.m$aav <- factor(rep(c("shRNC", "shR3"), c(12*14, 9*14)), levels=c("shRNC", "shR3"), labels=c("Control KD", "PPP2R5C KD"))
str(seahorse.m)


df = data.frame(x=c(62.1167,68.8167,89.15, 29, 58, 78), y=c(0.68,0.58,0,0.55,0.80,0.65), aav=c("shR3", "shR3", "shRNC", "shR3", "shR3", "shR3"), label=c("**","*", "**", "Glucose","Oligomycin", "2-Deoxy-\nGlucose"))
p <- ggplot(seahorse.m, aes(x=time, y=rate,group=aav))

# plot with t.test significance for mean gtt
set.seed(11)
p + stat_summary(fun.data = "mean_cl_boot", geom="errorbar", width=2) + stat_summary(fun.data = "mean_cl_boot", geom="line", aes(linetype=factor(aav))) + labs(x="Time (min)", y="ECAR (mpH/min)")  + scale_linetype_discrete(name="") + theme_thesis() + scale_y_continuous(expand = c(0, 0), limits=c(-0.3,1)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + geom_vline(xintercept=c(21.417,48.483,68.81), linetype="longdash", colour="grey80") + geom_text(data=df, mapping=aes(label=label,x,y,aav),size=4)

ggsave("fig2-11 seahorse.pdf", width=16,height=8,units="cm")


########################################
# fig2-12 2-13 facs filter and density
########################################

library(flowCore)
library(flowStats)
library(flowViz)
library(lattice)

load("~/Data/20130507 glucose uptake hepa16/20130507 2NBDG.RData")

pdf("./fig2-12 facs filter.pdf", width=8, height= 4)
xyplot(`SSC-A` ~ `FSC-W` | cell + trt , Data(wf[["base view"]])[c(1:2,5:8)],
       filter=filter(Data(wf[["base view"]])[c(1:2,5:8)],
                     norm2Filter("SSC-A", "FSC-W", filterId = "liveCell")),
       xlim=c(0,262144),ylim=c(0,262144)) # 262144 the limits of ssc and fsc
dev.off()

# grouped display of each experiment.
# refline could only add one line
pdf("./fig2-13 density 2NBDG.pdf", width=4, height=4)
densityplot(reorder(celltrt, index)~ `2NBDG`, 
            Data(wf[["2NBDG+"]])[c(1:2,5:8)],
            main = "Glucose Uptake Profiles",
            refline = mfi[c(29,1),4],
            subset = !(pd$name %in% c("Hepa16 Unstained")),
            xlim=c(6,10),
            xlab = "asinh Transformation of 2NBDG Intensity"
)
dev.off()




########################################
# fig2-14 facs measurement in starved inducible shRNA3 6 ***,**
#########################################
facs <- read.xlsx("~/Documents/Lab Meetings/Draft/FigData.xlsx", sheetName="Sheet2", rowIndex=65:71, colIndex=6:9,header=TRUE)
facs$cell <- factor(facs$cell, levels=c("Hepa 1-6", "shR3", "shR6"), labels=c("Hepa 1-6", "shR3-4", "shR6-5"))
facs$trt <- gl(2,1,2,labels=c("0x", "1x"))

df = data.frame(x=c(2.175,3.175), y=c(1.6,1.67), trt=c("1x", "1x"), label=c("***","**"))

p <- ggplot(facs, aes(x=cell, y=exp, fill=trt))

p + geom_bar(stat="identity", position="dodge", width=0.7) + geom_errorbar(aes(ymax=exp + se, ymin=exp - se), position=position_dodge(width=0.7), width=0.2) + scale_fill_grey(name="", start=0.8, end=0.2) + labs(x="", y="Relative 2NBDG\nUptake") + theme_thesis_it()  + scale_y_continuous(expand = c(0, 0), limits=c(0,2)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + geom_text(data=df, mapping=aes(label=label,x,y,trt),size=4)

ggsave("fig2-14 facs glucose uptake.pdf",width=16,height=8,units="cm")


########################################
# fig2-15 facs measure in primary hepatocyte ad-shR3 kd * for non-starve
########################################
hepa.1.facs <- read.xlsx("~/Data/20130929 glucose uptake primary hepa kd/glucose uptake primary hepa.xlsx", sheetName="glucose uptake primary hepa.csv", rowIndex=1:13, colIndex=4:6,header=TRUE)
hepa.1.facs$trt <- factor(hepa.1.facs$trt, levels=c("Non-Starved", "Starved"))
hepa.1.facs$aav <- factor(hepa.1.facs$aav, levels=c("PBS", "NC", "shR3"), labels=c("PBS","shRNC", "shR3"))

hepa.1.list <- split(hepa.1.facs[hepa.1.facs$aav!="PBS",], hepa.1.facs[hepa.1.facs$aav!="PBS",]$trt)
lapply(hepa.1.list, function(x) pairwise.t.test(x$nbdg, x$aav, p.adjust="none"))

df <- data.frame(x=3, y=1400, trt="Non-Starved", aav="shR3", label="*")
p <- ggplot(hepa.1.facs, aes(x=aav, y=nbdg, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5), width=0.2) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5, group=aav)) + scale_fill_grey(name="", start=0.8, end=0, labels=c("PBS", "shRNC", "shR3")) + labs(x="PPP2R5C KD in Primary Hepatocytes", y="2NBDG uptake (MFI)") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + facet_grid(.~trt) + geom_text(data=df, mapping=aes(label=label,x,y,trt),size=4)

ggsave("fig2-15 1st hepa glucose uptake facs.pdf",width=16,height=8,units="cm", device=cairo_pdf)


########################################
# fig2-16 TAG in 1st hepa 48 and 72 hrs * ** *** ns
########################################

temp <- read.xlsx("~/Data/20120221 TAG H48 24+48 K24+48/20120221 TAG H 48 24+48 K24+48.xls", sheetName="TAG", rowIndex=99:103, colIndex=2:13,header=TRUE)
temp <- as.data.frame(t(temp))
colnames(temp) <- c("trt","exp1","exp2","exp3")
temp$sample <- rep(c("K48", "K72", "H48", "H72"), each=3)

str(hepa.1.tag)
hepa.1.tag <- melt(temp, id=c(1,5))
hepa.1.tag <- hepa.1.tag[complete.cases(hepa.1.tag),]
hepa.1.tag$value <- as.numeric(hepa.1.tag$value) / 92.06 * 1000000
hepa.1.tag$sample <- factor(hepa.1.tag$sample, levels=c("K48","K72","H48","H72"), labels=c("Klingmueller 2 Day","Klingmueller 3 Day","Herzig 2 Day","Herzig 3 Day"))
hepa.1.tag$trt <- factor(hepa.1.tag$trt, levels=c("PBS","NC","shR3"), labels=c("PBS","shRNC","shR3"))

lapply(split(hepa.1.tag, hepa.1.tag$sample), function(x) pairwise.t.test(x$value, x$trt, p.adjust="BH"))

df <- data.frame(x=c(1.167,2.167,3.167), y=c(200, 148, 193), trt=c("shR3","shR3","shR3"), label=c("*", "**","***"))
p <- ggplot(hepa.1.tag, aes(x=sample, y=value, fill=trt))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5), width=0.2) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5, group=trt)) + scale_fill_grey(name="", start=0.8, end=0, labels=c("PBS", "shRNC", "shR3")) + labs(y="Relative Triglyceride\n(nmol/mg protein)", x="Triglyceride Storage upon PPP2R5C KD") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,250)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + geom_text(data=df, mapping=aes(label=label,x,y,trt),size=4)

ggsave("fig2-16 1st hepa tag.pdf",width=16,height=8,units="cm", device=cairo_pdf)

########################################
# fig2-17 ATP glycogen in 1st hepa ns
########################################

temp <- read.xlsx("~/Data/20110908 TAG Glycogen ATP level in primary hepa with ppp2r5c kd/20110908 glycogen.xls", sheetName="Sheet0", rowIndex=43:44, colIndex=2:7,header=TRUE)
hepa.1.glycogen <- as.data.frame(t(temp))
colnames(hepa.1.glycogen) <- c("glycogen")
hepa.1.glycogen$trt <- rep(c("PBS", "shRNC","shR3"), each=2)

str(hepa.1.glycogen)
hepa.1.glycogen$glycogen <- hepa.1.glycogen$glycogen * 1000 # ug/mg protein

hepa.1.glycogen$trt <- factor(hepa.1.glycogen$trt, levels=c("PBS","shRNC","shR3"))

with(hepa.1.glycogen, pairwise.t.test(glycogen, trt, p.adjust="BH"))

p <- ggplot(hepa.1.glycogen, aes(x=trt, y=glycogen, fill=trt))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5), width=0.2) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5)) + scale_fill_grey(name="", start=0.8, end=0, labels=c("PBS", "shRNC", "shR3")) + labs(y=expression("Relative Glycogen (" * mu * "g/mg protein)"), x="Adenovirus Infections") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5)))

ggsave("fig2-17a 1st hepa glycogen.pdf",width=8,height=8,units="cm", device=cairo_pdf)


temp <- read.xlsx("~/Data/20110908 TAG Glycogen ATP level in primary hepa with ppp2r5c kd/20110909 ATP.xls", sheetName="Sheet0", rowIndex=58:59, colIndex=2:7,header=TRUE)
hepa.1.atp <- as.data.frame(t(temp))
colnames(hepa.1.atp) <- c("atp")
hepa.1.atp$trt <- rep(c("PBS", "shRNC","shR3"), each=2)

str(hepa.1.atp)

hepa.1.atp$trt <- factor(hepa.1.atp$trt, levels=c("PBS","shRNC","shR3"))

with(hepa.1.atp, pairwise.t.test(atp, trt, p.adjust="BH"))

p <- ggplot(hepa.1.atp, aes(x=trt, y=atp, fill=trt))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5), width=0.2) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5)) + scale_fill_grey(name="", start=0.8, end=0, labels=c("PBS", "shRNC", "shR3")) + labs(y=expression("Relative ATP (nmol/mg protein)"), x="Adenovirus Infections") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5)))

ggsave("fig2-17b 1st hepa atp.pdf",width=8,height=8,units="cm", device=cairo_pdf)

########################################
# fig2-18 TAG in hepa with rapamycin 20nM * **
########################################

temp <- read.xlsx("~/Data/20120904 TAG measurement in hepa 1-6 with rapamycin treatment/20120831 TAG.xls", sheetName="Sheet1", rowIndex=83, colIndex=2:13,header=FALSE)
hepa.tag <- as.data.frame(t(temp))
colnames(hepa.tag) <- c("tag")
hepa.tag$trt <- factor(rep(c("PBS", "shRNC","shR3","PBS", "shRNC","shR3"), each=2), levels=c("PBS", "shRNC","shR3"))
hepa.tag$rapa <- factor(rep(c("DMSO", "Rapamycin"), each=6), level=c("DMSO", "Rapamycin"))

str(hepa.tag)


lapply(split(hepa.tag, hepa.tag$rapa), function(x) pairwise.t.test(x$tag, x$trt, p.adjust="BH"))

df <- data.frame(x=c(3,3), y=c(32,38.5), rapa=c("DMSO", "Rapamycin"),trt=c("shR3","shR3"), label=c("*", "**"))
p <- ggplot(hepa.tag, aes(x=trt, y=tag, fill=trt))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5), width=0.2) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5)) + scale_fill_grey(name="", start=0.8, end=0, labels=c("PBS", "shRNC", "shR3")) + labs(y=expression("Relative Triglyceride (nmol/mg protein)"), x="Adenovirus Infections") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0),limits=c(0,45)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + facet_grid(.~ rapa) + geom_text(data=df, mapping=aes(label=label,x,y,rapa),size=4)

ggsave("fig2-18 hepa tag rapamycin.pdf",width=16,height=8,units="cm")


########################################
# fig2-19 FA uptake and intracellular FA in hepa. ns
########################################

temp <- read.xlsx("~/Data/20120222 TAG glucose lactate in hepa 1-6 for profiling/20120222 TAG glucose lactate hepa 1-6.xls", sheetName="FFA", rowIndex=100:103, colIndex=2:4,header=TRUE)
hepa.ffa <- as.data.frame(t(temp))
colnames(hepa.ffa) <- c("exp1", "exp2", "exp3")
hepa.ffa$trt <- factor(c("PBS", "shRNC","shR3"), levels=c("PBS", "shRNC","shR3"))

hepa.ffa.m <- melt(hepa.ffa, id=4)

with(hepa.ffa.m, pairwise.t.test(value, trt, p.adjust="BH"))

p <- ggplot(hepa.ffa.m, aes(x=trt, y=value, fill=trt))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5), width=0.2) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5)) + scale_fill_grey(name="", start=0.8, end=0, labels=c("PBS", "shRNC", "shR3")) + labs(y=expression("Relative FFA (nmol/mg protein)"), x="Adenovirus Infections") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5)))

ggsave("fig2-19a hepa ffa uptake.pdf",width=8,height=8,units="cm", device=cairo_pdf)



temp <- read.xlsx("~/Data/20120305 TAG and FA in hepa 1-6 with KD and OE/20120305 hepa 16 OE free fatty acid.xls", sheetName="FA", rowIndex=70:73, colIndex=2:4,header=TRUE)
hepa.fa.in <- as.data.frame(t(temp))
colnames(hepa.fa.in) <- c("exp1", "exp2", "exp3")
hepa.fa.in$trt <- factor(c("PBS", "shRNC","shR3"), levels=c("PBS", "shRNC","shR3"))

hepa.fa.in.m <- melt(hepa.fa.in, id=4)

with(hepa.fa.in.m, pairwise.t.test(value, trt, p.adjust="BH"))

df <- data.frame(x=c(3), y=c(2.2), trt=c("shR3"), label=c("*"))
p <- ggplot(hepa.fa.in.m, aes(x=trt, y=value, fill=trt))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5), width=0.2) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5)) + scale_fill_grey(name="", start=0.8, end=0, labels=c("PBS", "shRNC", "shR3")) + labs(y=expression("Relative FFA (nmol/mg protein)"), x="Adenovirus Infections") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + geom_text(data=df, mapping=aes(label=label,x,y,trt),size=4)

ggsave("fig2-19b hepa ffa intracellular.pdf",width=8,height=8,units="cm", device=cairo_pdf)


###################################################
# fig2-20a import liver aav kd qPCR pilot data and plot
###################################################
pilot <- read.xlsx("~/Documents/Lab Meetings/Draft/FigData.xlsx", sheetName="Sheet2", rowIndex=43:45, colIndex=8:10,header=TRUE)

pilot$trt <- factor(pilot$trt, levels=c("miRNC", "miR12"))
pilot

p <- ggplot(pilot, aes(x=trt, y=exp))
p + geom_bar(stat="identity", position="dodge", width=0.47) + geom_errorbar(aes(ymax=exp + se, ymin=exp - se), position=position_dodge(width=0.47), width=0.2) + scale_fill_grey(name="", start=0.8, end=0.2) + labs(x="", y="PPP2R5C mRNA Levels\n(TBP normalized)") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,1.2)) + guides(colour = guide_legend(override.aes = list(size=0.5))) 
ggsave("fig2-20a pilot qPCR.pdf",width=8,height=6,units="cm")


########################################
# fig2-21 serum ALT
########################################

p <- ggplot(data.all, aes(x=trt, y=serum_ALT, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.47)) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=aav)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(y= "Serum ALT (U/L)", x="")  + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,50)) + guides(colour = guide_legend(override.aes = list(size=0.5)))  

ggsave("fig2-21 serum ALT.pdf", width=16,height=8,units="cm", device=cairo_pdf)

########################################
# fig2-22 body weight change
########################################
weight.melt <- melt(data.all[,1:16], id.vars=1:8, variable.name="Day", value.name="Weight")
weight.melt$Day <- rep(c(1,8,15,22,29,36,43,50), each=31)
p <- qplot(Day, Weight, data=weight.melt, facets=trt ~ aav, color=factor(gtt))
p + geom_line(aes(x=Day, y=Weight, group=mouse)) + scale_colour_discrete(name="GTT", labels=c("No GTT", "GTT")) + labs(title="Body weight change", y="Total Body Weight (g)")
p <- qplot(Day, Weight, data=weight.melt, facets=~aav)
p + geom_boxplot(aes(x=Day, y=Weight), alpha=0) 
p <- ggplot(weight.melt, aes(x=Day, y=Weight, group=aav, fill=aav))
# p + stat_summary(fun.data = "mean_cl_boot", geom="errorbar", width=2) + stat_summary(fun.data = "mean_cl_boot", geom="line", aes(linetype=aav)) + labs(x="Time Post AAV Injection (Day)", y="Body Weight (g)") + scale_linetype_discrete(name="") + theme_thesis() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + facet_grid(.~trt)
# ggsave("fig2-22 body weight.pdf", width=16,height=8,units="cm", device=cairo_pdf)

p4 <- p + stat_summary(fun.data = "mean_cl_boot", geom="errorbar", width=2) + stat_summary(fun.data = "mean_cl_boot", geom="line", aes(linetype=aav)) + labs(x="Time Post AAV Injection (Day)", y="Body Weight (g)") + scale_linetype_discrete(name="") + theme_thesis() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + theme(legend.position=c(1.2,0.8)) 
p4
ggsave("fig2-22a body weight.pdf", width=8,height=8,units="cm", device=cairo_pdf)

########################################
# fig2-23 echoMRI fat 
########################################
echo <- read.csv("~/Data/20140121-0312 AAV full experiment/echoMRI/echo.csv")
echo$aav <- factor(echo$aav, levels=c("miRNC", "miR12"))
echo$trt <- factor(echo$trt, levels=c("Random", "Fasted", "Refed"), labels=c("Random", "Fasting", "Refed"))
str(echo)
ind <- echo$Time=="0 Week" | echo$Time=="2 Week" | echo$Time=="4 Week" | echo$Time=="7 Week"
echo.7w <- echo[ind,]
echo.7wm <- melt(echo.7w, id=c(1:10,13)) 
echo.7wm$variable <- factor(echo.7wm$variable, levels=c("NormalFat", "NormalLean"), labels=c("Fat Mass", "Lean Mass"))
p <- ggplot(echo.7w, aes(x=Time, y=NormalFat, group=aav, fill=aav))
# p + stat_summary(fun.data = "mean_cl_boot", geom="errorbar", multi=1, width=0.2) + stat_summary(fun.data = "mean_cl_boot", geom="line", aes(linetype=factor(aav))) + scale_linetype_discrete(name="") + labs(y= "Normalized Fat Mass", x="Time after Virus Injection") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + facet_grid(.~ trt)
# 
# ggsave("fig2-23 echoMRI fat.pdf", width=16,height=8,units="cm", device=cairo_pdf)

p5 <- p + stat_summary(fun.data = "mean_cl_boot", geom="errorbar", multi=1, width=0.2) + stat_summary(fun.data = "mean_cl_boot", geom="line", aes(linetype=factor(aav))) + scale_linetype_discrete(name="") + labs(y= "Normalized Fat Mass", x="Time after Virus Injection") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + theme(legend.position=c(1.7,0.8)) 
p5
ggsave("fig2-22b echoMRI fat.pdf", width=8,height=8,units="cm", device=cairo_pdf)
p <- arrangeGrob(p4,p5,ncol=2)
ggsave("fig2-22 body weight and echoMRI fat.pdf", plot=p, width=16,height=8,units="cm", device=cairo_pdf)

########################################
# fig2-24 echoMRI lean 
########################################
p <- ggplot(echo.7w, aes(x=Time, y=NormalLean, group=aav, fill=aav))
p + stat_summary(fun.data = "mean_cl_boot", geom="errorbar", multi=1, width=0.2) + stat_summary(fun.data = "mean_cl_boot", geom="line", aes(linetype=factor(aav))) + scale_linetype_discrete(name="") + labs(y= "Normalized Lean Mass", x="Time after Virus Injection") + theme_thesis_it() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + facet_grid(.~ trt)

ggsave("fig2-24 echoMRI lean.pdf", width=16,height=8,units="cm", device=cairo_pdf)



########################################
# fig2-25 abdominal fat weight 
########################################

p <- ggplot(data.all, aes(x=trt, y=Abd_WAT_weight, group=aav, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5)) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=aav)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(y= "Abdominal WAT Weight (g)", x="")  + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,0.6)) + guides(colour = guide_legend(override.aes = list(size=0.5)))  

ggsave("fig2-25 abd wat weight.pdf", width=16,height=8,units="cm", device=cairo_pdf)

str(data.all)

##########################################
# fig2-26 blood glucose
#########################################
p <- ggplot(data.all, aes(x=trt, y=glucose, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5)) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=aav)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(y= "Blood Glucose (mg/dL)", x="")  + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,260)) + guides(colour = guide_legend(override.aes = list(size=0.5)))  

ggsave("fig2-26 blood glucose.pdf", width=16,height=8,units="cm" )

#############################
# Figure 2-27 liver weight * ** *
#############################
df <- data.frame(x=c(1.125,2.125,3.125), y=c(1.6,1.55,1.75), aav=c("miR12", "miR12", "miR12"), label=c("*","**","*"))
# plot liver weight
p <- ggplot(data.all, aes(x=trt, y=liver_weight, group=aav, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5)) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=aav)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(y= "Liver Weight (g)", x="")  + theme_thesis_it() + scale_y_continuous(expand = c(0, 0),limits=c(0,2.2)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + geom_text(data=df, mapping=aes(label=label,x,y,aav),size=4)

ggsave("fig2-27 liver weight.pdf", width=16,height=8,units="cm" )

##########################################
# fig2-28 serum insulin sensitivity increase
#########################################
df <- data.frame(x=c(1.125,3.125), y=c(0.6,1.25), aav=c("miR12", "miR12"), label=c("*","*"))

p <- ggplot(data.all, aes(x=trt, y=insulin_serum, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5)) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=aav)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(y= "Serum Insulin (ng/mL)", x="")  + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,2)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + geom_text(data=df, mapping=aes(label=label,x,y,aav),size=4)

ggsave("fig2-28 serum insulin.pdf", width=16,height=8,units="cm" )


##########################################
# fig2-29 homa index increase
#########################################
df <- data.frame(x=c(1.125), y=c(0.026), aav=c("miR12"), label=c("*"))

p <- ggplot(data.all, aes(x=trt, y=homa_index, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5)) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=aav)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(y= "HOMA Index", x="")  + theme_thesis_it() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + geom_text(data=df, mapping=aes(label=label,x,y,aav),size=4)

ggsave("fig2-29 homa index.pdf", width=16,height=8,units="cm" )

##########################################
# fig2-30 liver glycogen
#########################################
df <- data.frame(x=c(2.125, 3.125), y=c(40,55.5), aav=c("miR12","miR12"), label=c("**","*"))

p <- ggplot(data.all, aes(x=trt, y=glycogen_liver, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5)) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=aav)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(y= expression("Liver Glycogen (" * mu * "g/mg liver)"), x="")  + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,60)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + geom_text(data=df, mapping=aes(label=label,x,y,aav),size=4)

ggsave("fig2-30 liver glycogen.pdf", width=16,height=8,units="cm" )



##########################################
# fig2-31 liver TG
#########################################
df <- data.frame(x=c(1.125, 2.125), y=c(7.8,4.8), aav=c("miR12","miR12"), label=c("†","**"))

p <- ggplot(data.all, aes(x=trt, y=liver_TAG, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5)) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=aav)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(y= "Liver Triglyceride\n(nmol/mg liver)", x="")  + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,15)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + geom_text(data=df, mapping=aes(label=label,x,y,aav),size=4)

ggsave("fig2-31 liver tag.pdf", width=16,height=8,units="cm",device=cairo_pdf )

##########################################
# fig2-32 serum TG
#########################################
df <- data.frame(x=c(2.125, 3.125), y=c(2.2,2), aav=c("miR12","miR12"), label=c("**","***"))

p <- ggplot(data.all, aes(x=trt, y=serum_TAG, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5)) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=aav)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(y= "Serum Triglyceride\n(mmol/L)", x="")  + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,3)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + geom_text(data=df, mapping=aes(label=label,x,y,aav),size=4)

ggsave("fig2-32 serum tag.pdf", width=16,height=8,units="cm" )


##########################################
# fig2-33 liver NEFA
#########################################

p <- ggplot(data.all, aes(x=trt, y=liver_NEFA, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5)) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=aav)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(y= "Liver NEFA (nmol/mg)", x="")  + theme_thesis_it() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5)))

p <- ggplot(data.all, aes(x=trt, y=liver_NEFA, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5)) + stat_summary(fun.data = "mean_sdl", geom="bar",alpha=0.5, position="dodge", aes(width=0.5,group=aav)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(y= "Liver NEFA (nmol/mg)", x="")  + theme_thesis_it() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + geom_dotplot(binaxis = "y", stackdir = "center", binwidth=0.02, dotsize=3, position = position_dodge(width=0.5)) 

ggsave("fig2-33 liver NEFA.pdf", width=16,height=8,units="cm" )

##########################################
# fig2-34 serum NEFA
#########################################
df <- data.frame(x=c(3.125), y=c(0.65), aav=c("miR12"), label=c("*"))

p <- ggplot(data.all, aes(x=trt, y=serum_NEFA, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5)) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=aav)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(y= "Serum NEFA (mmol/L)", x="")  + theme_thesis_it() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + geom_text(data=df, mapping=aes(label=label,x,y,aav),size=4)

ggsave("fig2-34 serum NEFA.pdf", width=16,height=8,units="cm" )

##########################################
# fig2-35 liver cholesterol * **
#########################################
df <- data.frame(x=c(2.125, 3.125), y=c(1.9,2.15), aav=c("miR12", "miR12"), label=c("*", "**"))


p <- ggplot(data.all, aes(x=trt, y=liver_cholesterol, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5)) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=aav)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(y= "Liver Cholesterol\n(nmol/mg liver)", x="")  + theme_thesis_it() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + geom_text(data=df, mapping=aes(label=label,x,y,aav),size=4)

ggsave("fig2-35 liver chol.pdf", width=16,height=8,units="cm" )

##########################################
# fig2-36 serum cholesterol
#########################################
p <- ggplot(data.all, aes(x=trt, y=serum_cholesterol, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5)) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=aav)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(y= "Serum Cholesterol (mmol/L)", x="")  + theme_thesis_it() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5))) 

ggsave("fig2-36 serum chol.pdf", width=16,height=8,units="cm" )

##########################################
# fig2-37 serum TKB
#########################################
p <- ggplot(data.all, aes(x=trt, y=serum_TKB, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5)) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=aav)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(y= expression("Serum Total Ketone Body (" * mu * "mol/L)"), x="")  + theme_thesis_it() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5))) 

ggsave("fig2-37 serum TKB.pdf", width=16,height=8,units="cm" )


##########################################
# fig2-38 serum HB
#########################################
p <- ggplot(data.all, aes(x=trt, y=serum_HB, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5)) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=aav)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(y= expression("Serum Hydroxybutyrate (" * mu * "mol/L)"), x="")  + theme_thesis_it() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5))) 

ggsave("fig2-38 serum HB.pdf", width=16,height=8,units="cm" )

##########################################
# fig2-39 serum KB ratio
#########################################
p <- ggplot(data.all, aes(x=trt, y=KB_ratio, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.5)) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=aav)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(y= expression("HB/TKB Ratio"), x="")  + theme_thesis_it() + scale_y_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=0.5))) 

ggsave("fig2-39 serum KB ratio.pdf", width=16,height=8,units="cm" )


######################################
# fig2-40-42 serum FPLC profile
######################################
serum.fplc <- read.csv("~/Data/20140813-14 FPLC separation serum AAV/fplc_TgAndChol.csv")
fplc <- read.csv("~/Data/20140813-14 FPLC separation serum AAV/fplc.csv")
serum.fplc$trt <- factor(serum.fplc$trt, levels=c("Random", "Fasting", "Refed"))
serum.fplc$aav <- factor(serum.fplc$aav, levels=c("miRNC", "miR12"))
# labels=c("Control KD", "PPP2R5C KD")

fplc$trt <- factor(fplc$trt, levels=c("Random", "Fasting", "Refed"))
fplc$aav <- factor(fplc$aav, levels=c("miRNC", "miR12"))


# plot serum tag fplc
df <- data.frame(label = rep(c("VLDL","IDL","LDL", "HDL"),3), x = rep(c(22,29,36,41),3), y = c(0.072,0.022,0.016,0.022,0.115,0.022,0.019,0.025,0.14,0.032,0.025,0.028),trt=rep(c("Random","Fasting","Refed"),each=4),aav=rep("miR12",12),sample=rep(c("F","B","D"),each=4))
p <- ggplot(serum.fplc, aes(x=frac, y=tag, group=sample))
p + geom_line(aes(colour=aav)) + geom_point(aes(shape=aav), size=0.6)  + labs(x="Fractions", y="Triglycerid (mmol/L)")  + scale_shape_discrete(name="")  + scale_colour_manual(name="", values=c("grey60", "black"))+ annotate("text", size =2, ) + theme_thesis() + scale_y_continuous(expand = c(0, 0), limits=c(-0.05,0.17), breaks=c(0,0.05,0.1, 0.15)) + scale_x_continuous(expand = c(0, 0), limits=c(13,45), breaks=c(15,20,25,30,35,40,45))+ guides(colour = guide_legend(override.aes = list(size=0.5)))  + facet_grid(.~ trt) + geom_text(data=df, mapping=aes(label=label,x,y,trt,aav,sample),size=2)

ggsave("fig2-40 fplc_tag.pdf", width=16, height=8, units="cm")

# plot serum cholesterol fplc
df <- data.frame(label = rep(c("VLDL","IDL","LDL", "HDL"),3), x = rep(c(22,29,36,41),3), y = c(0.022,0.022,0.12,0.03,0.022,0.022,0.078,0.02,0.02,0.02,0.09,0.02),trt=rep(c("Random","Fasting","Refed"),each=4),aav=rep("miR12",12),sample=rep(c("F","B","D"),each=4))
p <- ggplot(serum.fplc, aes(x=frac, y=chol, group=sample))
p + geom_line(aes(colour=aav)) + geom_point(aes(shape=aav), size=0.6)  + labs(x="Fractions", y="Cholesterol (mmol/L)")  + scale_shape_discrete(name="")  + scale_colour_manual(name="", values=c("grey60", "black")) + theme_thesis() + scale_y_continuous(expand = c(0, 0), limits=c(-0.01,0.13), breaks=c(0,0.05, 0.1)) + scale_x_continuous(expand = c(0, 0), limits=c(13,45), breaks=c(15,20,25,30,35,40,45)) + guides(colour = guide_legend(override.aes = list(size=0.5)))  + facet_grid(.~trt) + geom_text(data=df, mapping=aes(label=label,x,y,trt,aav,sample),size=2)

ggsave("fig2-41 fplc_chol.pdf", width=16, height=8, units="cm")

# plot serum fplc profile
str(fplc)
p <- ggplot(fplc, aes(x=ml, y=X.mAU, group=aav))
p + geom_line(aes(colour=aav))  + labs(x="Volume (mL)", y="Absorbance 280 nm (mAU)")  + scale_colour_manual(name="", values=c("grey60", "black"))+ annotate("text", size =2,  label = c("VLDL","IDL","LDL", "HDL"), x = c(7.7, 12.5,15.2, 17), y = c(500,550, 800, 2200)) + theme_thesis() + scale_y_continuous(expand = c(0, 0), limits=c(-100,2500), breaks=c(0,1000, 2000)) + scale_x_continuous(expand = c(0, 0), limits=c(6.5,20), breaks=c(10,15,20)) + guides(colour = guide_legend(override.aes = list(size=0.5)))  + facet_grid(.~trt)

ggsave("fig2-42 fplc_profile.pdf", width=16, height=8, units="cm")



##########################################
# fig2-43 GTT
#########################################
# gtt assay data processing
gttres <- read.xlsx("~/Data/20140121-0312 AAV full experiment/GTT/20140225-26GTT YC.xlsx", sheetName="all data", rowIndex=1:25,colIndex=1:11)
gttres$AAV <- factor(gttres$AAV, levels=c("miRNC", "miR12"))
gttres$Mouse <- factor(gttres$Mouse)
gttres$gtt_time <- factor(gttres$gtt_time, labels=c("2014-02-25","2014-02-26"))
library(reshape2)
gttres.m <- melt(gttres, id=1:5)
colnames(gttres.m) <- c("aav","weight","label","gtt_time","operator","time","glucose")
gttres.m$time <- rep(c(0,20,60,90,120,150),each=24)
gttres.m$idx <- rep(factor(gttres$X20 > 410, labels=c("GTT_20 min < 410 mg/dL", "GTT_20 min > 410 mg/dL")), 6)

# bootstrap CI for all means from Hmisc, set same seed as plot
# use upper miRNC for label significance stars
library(Hmisc)
list.time <- split(gttres.m, f=gttres.m$time)
set.seed(11)
mean.cl.boot.list <- lapply(list.time, function(x) smean.cl.boot(x[x$aav == "miRNC",]$glucose))
mirnc.upper <- sapply(mean.cl.boot.list, "[", 3)


p <- ggplot(gttres.m, aes(x=time, y=glucose, group=factor(aav)))

# plot with t.test significance for mean gtt
df <- data.frame(x=c(60,90,120,150), y=c(392,309,262,215), label=c("*", "**"))


set.seed(11)
p + stat_summary(fun.data = "mean_cl_boot", geom="errorbar", width=5) + stat_summary(fun.data = "mean_cl_boot", geom="line", aes(linetype=factor(aav))) + labs(x="Time after Injection (min)", y="Blood Glucose (mg/dL)")  + scale_linetype_discrete(name="") +  annotate("text", x=c(0,20,60,90,120,150), y=c(mirnc.upper + 10), label=c(" "," ","**","*","*","*")) + theme_thesis() + scale_y_continuous(expand = c(0, 0), limits=c(0,500)) + guides(colour = guide_legend(override.aes = list(size=0.5)))  

ggsave("fig2-43 gtt.pdf", width=16,height=8,units="cm")


##########################################
# fig2-44 auc for gtt
##########################################
# AUC analysis 
library(pracma)
list <- split(gttres.m, f=gttres.m$label)
auc <- sapply(list, function(x) trapz(x$time, x$glucose))
names(auc)
auc2 <- data.frame(auc=auc[gttres$Mouse],aav=gttres$AAV)
with(auc2, t.test(auc[aav=="miRNC"], auc[aav=="miR12"])) # pvalue 0.008216

p1 <- ggplot(data=auc2, aes(y=auc, x=aav))
p1 + geom_boxplot(alpha=0.2) + geom_point()  + labs(y=expression("AUC for GTT"), x="") + annotate("text", x=c(2), y=c(45000), label=c("**")) + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(21000,54000)) + guides(colour = guide_legend(override.aes = list(size=0.5))) 
ggsave("fig2-44 AUC for gtt.pdf", width=8,height=8,units="cm" )

##########################################
# fig2-45 serum insulin for gtt
##########################################
serumIns.melt <- read.csv("~/Data/20140430 serum Insulin GTT aav big/serumIns.csv")
temp <- split(serumIns.melt, serumIns.melt$Time)
lapply(temp, function(x) t.test(repAverage ~ aav, data=x))
serumIns.melt$aav <- factor(serumIns.melt$aav, levels=c("miRNC", "miR12")) 
# * ns ** 

# repeated measure ANOVA analysis for GTT serum insulin
str(serumIns.melt)
serumIns.melt$timef <- as.factor(serumIns.melt$Time)
serumIns.melt$labelf <- as.factor(serumIns.melt$mouse)

a <- aov(repAverage ~ timef*aav + Error(labelf),data=serumIns.melt)
summary(a) # p-value 0.00611

summary(serumIns.melt$repAverage)
p3 <- ggplot(serumIns.melt, aes(x=Time, y=repAverage, group=aav))
p3 + stat_summary(fun.data = "mean_cl_boot", geom="errorbar", width=5) + stat_summary(fun.data = "mean_cl_boot", geom="line", aes(linetype=factor(aav))) +scale_linetype_discrete(name="")+ labs(x="Time after Injection (min)", y="O.D.450") + annotate("text", x=c(0,60), y=c(0.01235, 0.013), label=c("*","**")) + theme_thesis() + scale_y_continuous(expand = c(0, 0),limits=c(0,0.017)) + guides(colour = guide_legend(override.aes = list(size=0.5))) 
ggsave("fig2-45 serum ins for gtt.pdf", width=12,height=8,units="cm" )

#########################################
# Figure 2-53 qpcr of hif1a target in 1st hepa
# t.test ***, *, **, *,
#########################################
hif.target <- read.xlsx("~/Documents/Lab Meetings/Draft/FigData.xlsx", sheetName="Sheet2", rowIndex=160:168, colIndex=1:4,header=TRUE)
hif.target$gene <- factor(hif.target$gene, levels=c("NDRG","LDHa","HK2","PKM2"), labels=c("NDRG","LDHa","HK2","PKM2"))
hif.target$trt <- factor(hif.target$trt, levels=c("shRNC", "shR3"))

p <- ggplot(hif.target, aes(x=gene, y=exp, fill=trt))

p + geom_bar(stat="identity", position="dodge", width=0.7) + geom_errorbar(aes(ymax=exp + se, ymin=exp - se), position=position_dodge(width=0.7), width=0.2) + scale_fill_grey(name="", start=0.8, end=0.2) + labs(x="", y="HIF1a Targets mRNA\n Levels (TBP normalized)") + annotate("text", x=c(1.175,2.175,3.175,4.175), y=c(3,1.9,2.8, 1.7), label=c("***", "*","**","*")) + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,4.5)) + guides(colour = guide_legend(override.aes = list(size=0.5))) 
ggsave("fig2-53 hif target 1st hepa.pdf",width=12,height=8,units="cm")

#########################################
# Figure 2-54 qpcr of srebp target in 1st hepa
# t.test *, *, ***, **, ***,***
#########################################
srebp.target <- read.xlsx("~/Documents/Lab Meetings/Draft/FigData.xlsx", sheetName="Sheet2", rowIndex=104:116, colIndex=1:4,header=TRUE)
srebp.target$gene <- factor(srebp.target$gene, levels=c("HMGCS1","LDLR","SLC25A1","DGAT2","ACLY", "GPAT1"), labels=c("HMGCS1","LDLR","SLC25A1","DGAT2","ACLY", "GPAT1"))
srebp.target$trt <- factor(srebp.target$trt, levels=c("shRNC", "shR3"))

p <- ggplot(srebp.target, aes(x=gene, y=exp, fill=trt))

p + geom_bar(stat="identity", position="dodge", width=0.7) + geom_errorbar(aes(ymax=exp + se, ymin=exp - se), position=position_dodge(width=0.7), width=0.2) + scale_fill_grey(name="", start=0.8, end=0.2) + labs(x="", y="SREBP Targets mRNA \n Levels (TBP normalized)") + annotate("text", x=c(1.175,2.175,3.175,4.175,5.175,6.175), y=c(3.2,2.3,1.65,1.7,2.7,1.9), label=c("*", "*","***","**","***","**")) + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,4.5)) + guides(colour = guide_legend(override.aes = list(size=0.5))) 
ggsave("fig2-54 srebp 1st hepa.pdf",width=12,height=8,units="cm")

# t.test for qpcr data

srebp.t <- read.xlsx("~/Documents/Lab Meetings/Draft/FigData.xlsx", sheetName="Sheet4", rowIndex=1:9, colIndex=1:8,header=TRUE)
str(srebp.t)
sapply(srebp.t[,3:8], function(x) pairwise.t.test(x, factor(srebp.t$trt), p.adj = "none"))



#########################################
# Figure 2-55 qpcr of srebp target in Liver
# t.test *, *, ***, **, ***,***
#########################################

srebp.liver1 <- read.csv("~/Data/20140825 qPCR of Hif and SREBP target in aav big liver/srebp.csv")
srebp.liver2 <- read.csv("~/Data/20141018 qPCR pck1 g6pc srebp in aav big liver/qpcr results.csv")
srebp.liver2$aav <- factor(srebp.liver2$aav, levels=c("miRNC", "miR12"),labels=c("Control KD", "PPP2R5C KD"))

srebp.liver2[srebp.liver2$variable=="SREBP",]
srebp.liver <- rbind(srebp.liver1[,-1], srebp.liver2[srebp.liver2$variable=="SREBP",-1])

srebp.liver$variable <- factor(srebp.liver$variable, levels=c("SLC25A1", "DGAT2", "ACLY", "GPAT1","HILPDA", "SREBP"))
srebp.liver$trt <- factor(srebp.liver$trt, levels=c("Random","Fasting","Refed"))
srebp.liver$aav <- factor(srebp.liver$aav, levels=c("Control KD", "PPP2R5C KD"), labels=c("miRNC","miR12"))

df = data.frame(x=c(2.175,3.175,2.175,3.175,2.175,3.175,2.175,3.175,1.175,3.175), y=c(2.25,4,2.4,4,3,8,1.6,1.4,1.25,1.8), variable=c("SREBP","SREBP","ACLY","ACLY","GPAT1", "GPAT1","SLC25A1","SLC25A1","DGAT2","DGAT2"), trt=c("Fasting", "Refed","Fasting", "Refed", "Fasting", "Refed","Fasting", "Refed","Random","Refed"),aav=rep("miR12", 10),label=c("***","**","**","***", "***","**","**","***","***","**"))
p <- ggplot(srebp.liver[srebp.liver$variable!="HILPDA",], aes(x=trt, y=value, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.7)) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.7,group=aav)) + scale_fill_grey(name="", start=0.8, end=0.2) + labs(x="", y="SREBP Targets mRNA \n Levels (TBP normalized)") + theme_thesis_it() + guides(colour = guide_legend(override.aes = list(size=0.5))) + facet_wrap( ~ variable, ncol=2, scales="free_y", as.table=FALSE) + geom_text(data=df, mapping=aes(label=label,x,y,variable,trt),size=4) + theme(legend.position=c(1.3,0.9))                                                                                                                                                                                                                                                                                                                                                                                                                                                              
ggsave("fig2-55 liver srebp.pdf",width=12,height=13,units="cm")






###################################################
# fig2-57 pilot TAG liver **
###################################################
pilot.tg <- read.xlsx("/Users/yongsheng/Documents/Lab Meetings/Draft/FigData.xlsx", sheetName="Sheet2", rowIndex=174:184, colIndex=1:5,header=TRUE)

pilot.tg$aav <- factor(pilot.tg$aav, levels=c("miRNC-H", "miR12-H", "miRNC-L", "miR12-L"), labels=c("Control KD H", "PPP2R5C KD H", "Control KD L", "PPP2R5C KD L"))
pilot.tg$virus <- factor(pilot.tg$virus, levels=c("miRNC", "miR12"))

pilot.tg

p <- ggplot(pilot.tg, aes(x=virus, y=tag))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.47)) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5)) + scale_fill_grey(name="", start=0.6, end=0.2) + labs(y= "Liver Triglyceride\n(nmol/mg liver)", x="")  + annotate("text", x=c(2), y=c(7.8), label=c("**")) + theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,10)) + guides(colour = guide_legend(override.aes = list(size=0.5)))  

ggsave("fig2-57 liver pilot tag.pdf", width=8,height=8,units="cm" )

with(pilot.tg, pairwise.t.test(tag, factor(virus)))

#########################################
# fig 2-58 peak apex and area for FPLC serum VLDL
#########################################
vldl <- read.csv("~/Data/20140813-14 FPLC separation serum AAV/vldl_peak.csv")
vldl.m <- melt(vldl[,-1], id=c(3:5))
vldl.m$trt <- factor(vldl.m$trt, levels=c("Random", "Fasting", "Refed"))
vldl.m$aav <- factor(vldl.m$aav, levels=c("miRNC", "miR12"))
vldl.m$variable <- factor(vldl.m$variable, levels=c("apexIntensity", "peakArea"), labels=c("Peak Apex","Peak Area"))
vldl.m$method <- factor(vldl.m$method, levels=c("FPLC","TG"), labels=c("FPLC Profile", "Triglyceride"))

p <- ggplot(vldl.m, aes(x=trt, y=value, fill=aav))
p + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.5,group=aav)) + scale_fill_grey(name="", start=0.6, end=0.2) + facet_wrap(method ~ variable, scale="free_y") + labs(x="", y="Relative Intensity") +  theme_thesis_it() + guides(colour = guide_legend(override.aes = list(size=0.5)))

ggsave("fig2-58 vldl_peak.pdf", width=16,height=8,units="cm")


###################################################
# fig2-60--64 human liver ppp2r5c in T2D ***
###################################################
ppp2r5c.human <- read.xlsx("~/Data/20141119 human liver T2D ppp2r5c qPCR data from leipzig//PPP2R5C data Herzig.xls", sheetName="Tabelle1", colIndex=1:30, rowIndex=1:67, header=TRUE, colClasses=c("numeric","character", "numeric", rep("character",3),rep("numeric",24)))
str(ppp2r5c.human)
summary(ppp2r5c.human)
ppp2r5c.human$Diabetes <- factor(ppp2r5c.human$Diabetes, levels=c("no","T2D"), labels=c("Control", "T2D"))
colnames(ppp2r5c.human) <- c("PatNo", "DOB", "PPP2R5C", "Diabetes", "Gender", "Group", "VATarea","SATarea","BMI","WHR", "VATdiameter","SATdiameter","BodyFat","FPG","FPI","2hrOGTT","ClampGIR","HbA1c","Chol1","Chol2","HDL.Chol1", "HDL.Chol2", "LDL.Chol1", "LDL.Chol2","TG1", "TG2", "FFA","Leptin","Adip","IL6")
pairwise.t.test(ppp2r5c.human$PPP2R5C, ppp2r5c.human$Diabetes, p.adj=c("none"))

p1 <- ggplot(ppp2r5c.human, aes(y=PPP2R5C, x=Diabetes))

p1 <- p1 + geom_boxplot(size=0.3, outlier.size=1) + labs(y= "Human Liver PPP2R5C \n(18S rRNA Normalized)", x="")  + theme_thesis() + scale_y_continuous(expand = c(0, 0), limits=c(0,30)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + annotate("text", size = 3,  label = c("***","p=0.0003"), x = c(2, 1.5), y = c(28,15))
p1
ggsave("fig2-60a human liver ppp2r5c.pdf",width=8,height=8,units="cm")

range(ppp2r5c.human$ClampGIR)
p2 <- ggplot(ppp2r5c.human, aes(y=rank(PPP2R5C), x=ClampGIR))
p2 <- p2 + geom_point(aes(colour=Diabetes),size=0.7)+ labs(x=expression("Clamp GIR ("*mu*"mol/kg/min)"), y="Liver PPP2R5C Rank\n(18S rRNA normalized)") + geom_smooth(method="lm") + scale_colour_grey(name="", start=0.8, end=0.2) + annotate("text", size = 2, label = c("r=-0.34,\np-val=0.0047"), x = c(105), y = c(65)) + theme_thesis() + theme(legend.position=c(1.8,0.8))+ scale_y_continuous(expand = c(0, 0),limits=c(0,80)) + scale_x_continuous(expand = c(0, 0),limits=c(0,120)) + guides(colour = guide_legend(override.aes = list(size=1.5)))
p2
ggsave("fig2-60b human ppp2r5c vs ClampGIR.pdf", width=8,height=8,units="cm")
p <- arrangeGrob(p1,p2,ncol=2)
ggsave("fig2-60 human ppp2r5c vs t2d GIR.pdf",plot=p, width=16, height=8, units="cm")

lapply(split(ppp2r5c.human, ppp2r5c.human$Diabetes),function(x) TukeyHSD(aov(PPP2R5C ~ Group, data=x))) # sig in ctr: VIS-Lean,p 0.0008385;VIS-SC,p 0.0030885
ppp2r5c.human$Group <- factor(ppp2r5c.human$Group, levels=c("Lean","SC","VIS"), labels=c("Lean", "SC Obesity","VIS Obesity"))

df = data.frame(x=c(1,1,3,3,3,2,2),y=c(2.5,10,10,5.5,7.5,7.5,3.5),Diabetes=rep(c("Control","Control"), c(4,3)))
df2 = data.frame(x=c(2,2.5), y=c(11,8.5), Diabetes=c("Control","Control"),label=c("†††","††"))
p <- ggplot(ppp2r5c.human, aes(y=PPP2R5C, x=Group))
p + geom_boxplot() + facet_grid(.~Diabetes)+ labs(x="",y="Human Liver PPP2R5C\n(18S rRNA normalized)")+ theme_thesis_it() + scale_y_continuous(expand = c(0, 0), limits=c(0,30)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + geom_path(data=df, mapping=aes(x,y,Diabetes)) + geom_text(data=df2, mapping=aes(label=label,x,y,Diabetes),size=4) 
ggsave("fig2-61 human ppp2r5c vs ob group.pdf",width=16, height=8, units="cm", device=cairo_pdf)

list.T2D <- split(ppp2r5c.human,ppp2r5c.human$Diabetes)
list <- list()
list.cor.total <- list()
# pearson and spearman method for correlation and test
for (i in 1:2) {
    temp <- list.T2D[[i]]
    temp2 <- list()
    for (j in 5:30) {
        r <- cor(temp[,3],as.numeric(temp[,j]),use="pairwise.complete.obs")
        test <- cor.test(temp[,3], as.numeric(temp[,j]),use="pairwise.complete.obs")
        pval <- test$p.value
        r.spearman <- cor(temp[,3],as.numeric(temp[,j]),use="pairwise.complete.obs",method="spearman")
        test.spearman <- cor.test(temp[,3], as.numeric(temp[,j]),use="pairwise.complete.obs",method="spearman")
        pval.spearman <- test.spearman$p.value
        r.total <- cor(ppp2r5c.human[,3], as.numeric(ppp2r5c.human[,j]), use="pairwise.complete.obs")
        test.total <- cor.test(ppp2r5c.human[,3], as.numeric(ppp2r5c.human[,j]),use="pairwise.complete.obs")
        pval.total <- test.total$p.value
        r.total.spearman <- cor(ppp2r5c.human[,3], as.numeric(ppp2r5c.human[,j]), use="pairwise.complete.obs",method="spearman")
        test.total.spearman <- cor.test(ppp2r5c.human[,3], as.numeric(ppp2r5c.human[,j]),use="pairwise.complete.obs",method="spearman")
        pval.total.spearman <- test.total.spearman$p.value
        temp2[[j]] <- data.frame(coef=r, pvalue=pval,coef.spear=r.spearman, pvalue.spear=pval.spearman)
        list.cor.total[[j]] <- data.frame(coef=r.total, pvalue=pval.total,coef.spear=r.total.spearman, pvalue.spear=pval.total.spearman)
    }
    list[[i]] <- temp2
}

list
list.cor.total
names(list.T2D)
corResults <- do.call(rbind, lapply(list, function(x) do.call(rbind, x)))
corResTotal <- rbind(corResults, do.call(rbind, list.cor.total))

corResTotal$Diabetes <- rep(c(names(list.T2D),"Total"), each=26)
corResTotal$Variables <- colnames(ppp2r5c.human)[5:30]
nrRes <- corResTotal[c(c(1:15,17,19,21,23:26), c(1:15,17,19,21,23:26) + 26, c(1:15,17,19,21,23:26) + 52),]

# plot pearson method result
sigIdx <- nrRes$pvalue[45:66] < 0.05
corResNew <- nrRes[rep(sigIdx,3),]
corResNew$Variables <- factor(corResNew$Variables, levels=c("Group","VATarea","HbA1c","TG1","2hrOGTT","FPG","ClampGIR","Chol1","LDL.Chol1","Leptin","IL6"), labels=c("Group","VATarea","HbA1c","Triglyceride","2hrs OGTT","FPG","ClampGIR","Cholesterol","LDL Cholesterol","Leptin","IL6"))

sigStar <- corResNew[corResNew$pvalue < 0.05,c(1,2,5,6)]

sigStar$star <- ifelse(sigStar$pvalue < 0.05 & sigStar$pvalue > 0.01, "*", ifelse(sigStar$pvalue < 0.01 & sigStar$pvalue > 0.001, "**", ifelse(sigStar$pvalue < 0.001,"***","***")))

sigStar$y <- ifelse(sigStar$coef >= 0, sigStar$coef + 0.02, sigStar$coef - 0.04)
sigStar$x <- ifelse(sigStar$Diabetes == "Control", as.numeric(sigStar$Variables) - 0.2, ifelse(sigStar$Diabetes == "T2D", as.numeric(sigStar$Variables), ifelse(sigStar$Diabetes == "Total", as.numeric(sigStar$Variables) + 0.2, as.numeric(sigStar$Variables) + 0.2)))
sigStar
p <- ggplot(corResNew, aes(x=Variables, y=coef, fill=Diabetes))
p + stat_summary(fun.data = "mean_sdl", geom="errorbar", width=0.2, mult=1, position = position_dodge(width=0.6)) + stat_summary(fun.data = "mean_sdl", geom="bar",position="dodge", aes(width=0.6, group=Diabetes)) + scale_fill_grey(name="", start=0.8, end=0.2) + labs(y= "PPP2R5C Correlation Coefficients", x="")  + theme_thesis_it() + guides(colour = guide_legend(override.aes = list(size=0.5))) + annotate("text", size = 3, label = sigStar$star, x = sigStar$x, y = sigStar$y)
ggsave("fig2-62 corr change.pdf", width=16, height=8, units="cm")


p2 <- ggplot(ppp2r5c.human, aes(y=rank(PPP2R5C), x=FPG))
p2 <- p2 + geom_point(aes(colour=Diabetes),size=0.7)+ labs(x=expression("Fasting Plasma Glucose (mmol/L)"), y="Liver PPP2R5C Rank\n(18S rRNA normalized)") + geom_smooth(method="lm") + scale_colour_grey(name="", start=0.8, end=0.2) + annotate("text", size = 2, label = c("r=0.36,\np-val=0.006"), x = c(5), y = c(70)) + theme_thesis() + theme(legend.position=c(1.8,0.8))+ scale_y_continuous(expand = c(0, 0),limits=c(0,80)) + guides(colour = guide_legend(override.aes = list(size=1.5)))
p2
ggsave("fig2-63a human ppp2r5c vs fpg.pdf", width=8,height=8,units="cm")
range(ppp2r5c.human$IL6,na.rm=TRUE)
table(ppp2r5c.human$Diabetes)
p3 <- ggplot(ppp2r5c.human, aes(y=rank(PPP2R5C), x=IL6))
p3 <- p3 + geom_point(aes(colour=Diabetes),size=0.7)+ labs(x=expression("Serum IL6 (pmol/L)"), y="Liver PPP2R5C Rank\n(18S rRNA normalized)") + geom_smooth(method="lm") + scale_colour_grey(name="", start=0.8, end=0.2) + annotate("text", size = 2, label = c("r=0.42,\np-val=0.0006"), x = c(2), y = c(70)) + theme_thesis() + theme(legend.position=c(1.8,0.8))+ scale_y_continuous(expand = c(0, 0),limits=c(0,80)) + scale_x_continuous(expand = c(0, 0)) + guides(colour = guide_legend(override.aes = list(size=1.5)))
p3
ggsave("fig2-63b human ppp2r5c vs IL6.pdf", width=8,height=8,units="cm")
p <- arrangeGrob(p2,p3,ncol=2)
ggsave("fig2-63 human ppp2r5c vs fpg IL6.pdf",plot=p, width=16, height=8, units="cm")


#########################################
# fig2-64 human adipose ppp2r5c 
#########################################
ppp2r5c.adipose <- read.xlsx("~/Data/20141203 human adipose PPP2R5C/PPP2R5C_jvo.xlsx", sheetName="Sheet1", colIndex=1:23, rowIndex=1:48, header=TRUE, colClasses=rep("numeric",23))

ppp2r5c.adipose$dm <- factor(ppp2r5c.adipose$dm, labels=c("Control", "T2D"))
pairwise.t.test(ppp2r5c.adipose$PPP2R5C.RQsc, ppp2r5c.adipose$dm) # 0.022
table(ppp2r5c.adipose$dm) # 36 control, 11 T2D
pairwise.t.test(ppp2r5c.adipose$PPP2R5C.RQv, ppp2r5c.adipose$dm) # 0.79

plot.t2d <- ppp2r5c.adipose[,c(5,22,23)]
plot.t2d.m <- melt(plot.t2d,id=1)
plot.t2d.m$variable <- factor(plot.t2d.m$variable, levels=c("PPP2R5C.RQsc", "PPP2R5C.RQv"), labels=c("SAT", "VAT"))

df = data.frame(x=c(0.775,0.775,1.775,1.775),y=c(2.2,3,3,2.5))
df2 = data.frame(x=c(1.275), y=c(2.8), label=c("p=0.02"))
p <- ggplot(plot.t2d.m, aes(y=value, x=dm))
p + geom_boxplot(aes(fill=variable), position = position_dodge(width = 0.9), size=0.3, outlier.size=1.5) + labs(x="", y="PPP2R5C mRNA Levels\n(PPIA normalized)") + scale_fill_grey(name="", start=1, end=0.8) + theme_thesis() + scale_y_continuous(expand = c(0, 0),limits=c(0,3)) + guides(colour = guide_legend(override.aes = list(size=0.5))) + geom_path(data=df, mapping=aes(x,y)) + geom_text(data=df2, mapping=aes(label=label,x,y),size=4) 

ggsave("fig2-64 human adipose ppp2r5c.pdf",width=12,height=8,units="cm")












# Figure 4
#########################################
# fig 4b GO enrich from bioID
#########################################
go.function <- read.xlsx("./FigData.xlsx", sheetName="Sheet2", rowIndex=147:154, colIndex=2:5,header=FALSE)
colnames(go.function) <- c("GO Terms", "gene", "perGene", "perPath")

p <- ggplot(go.function, aes(x=factor(1),y=perGene, fill=factor(go.function[,1]))) + geom_bar(stat="identity", width = 1) 

p + coord_polar(theta="y") + xlab('') + ylab('') + labs(fill="go") + theme_bw() + theme(axis.ticks = element_blank(), axis.text = element_blank(), legend.position="right", legend.text = element_text(size = 5),legend.background = element_rect(fill="transparent"), legend.key.size = unit(0.3, "cm"), legend.title=element_text(size=5)) + theme(axis.line = element_line(colour = "NA"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + scale_y_continuous(expand = c(0, 0))  +  scale_fill_discrete(name="GO Terms")

ggsave("fig4b.pdf",width=12,height=5.5,units="cm")

go.path <- read.xlsx("./FigData.xlsx", sheetName="Sheet2", rowIndex=122:144, colIndex=2:5,header=FALSE)
colnames(go.path) <- c("go", "gene", "perGene", "perPath")

p <- ggplot(go.path, aes(x=factor(1),y=perGene, fill=factor(go))) + geom_bar(width = 1) 

p + coord_polar(theta="y") + xlab('') + ylab('') + labs(fill="go") + theme_bw() + theme(axis.ticks = element_blank(), axis.text = element_blank(), legend.position="right", legend.text = element_text(size = 5),legend.background = element_rect(fill="transparent"), legend.key.size = unit(0.3, "cm"), legend.title=element_text(size=5)) + theme(axis.line = element_line(colour = "NA"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + scale_y_continuous(expand = c(0, 0))  +  scale_fill_discrete(name="GO Terms") 

ggsave("fig4b GO path.pdf",width=20,height=12,units="cm")








