**Step 1: Prepare the Environment in R**

These commands will load various libraries into the active R environment. If you receive an error that one of these libraries isn't recognized, simply use the command install.packages()`

```
library(RCurl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(cowplot)
```

**Step 2: Acquire Data to Recreate Figure 1 from Github**

For this one we are going to build four graphs and then combine them together at the end. So first, we need to load data from two different sources. Then we'll pull the data from the second data set
into three different subsets for each figure B through D.

```
Fig4A_data<-read.csv(text=getURL("https://raw.githubusercontent.com/surtlab/data_for_figures/master/Final_Putida_pMP_comps_Results.csv"))
Fig4BCD_data<-read.csv(text=getURL("https://raw.githubusercontent.com/surtlab/data_for_figures/master/mSphere_megaplasmid_2020/pMPPutida_cross-enviroment-relative_fitness.csv"))
Fig4B_data<-subset(Fig4BCD_data,Experiment=="Nal")
Fig4C_data<-subset(Fig4BCD_data,Experiment=="Cipro")
Fig4D_data<-subset(Fig4BCD_data,Experiment=="Temp")

```
**Step 3: Create the Plot Data Order**

This step does a couple of different things at once. First, note the `%>%` function which acts as a pipe between commands. Then I'm using 'mutate' to reorder the data on the x-axis. I do this for both data sets for consistency, and call the data for each graph into a different variable, `p` or `q`

```
p<- Fig4A_data %>% mutate(Strain_Pair=recode(Strain_Pair,'C305-1604'="DBL1604",'D305-1620'="DBL1620")) 
q<- Fig4B_data %>% mutate(Treatment=recode(Treatment,'Nal 0 ng/uL'="0",'Nal 20 ng/uL'="20",'Nal 40 ng/uL'="40"))
r<- Fig4C_data %>% mutate(Treatment=recode(Treatment,'Cip 0 ng/uL'="0",'Cip 0.5 ng/uL'="0.5"))
s<- Fig4D_data %>% mutate(Treatment=recode(Treatment,'27oC'="27",'37oC'="37"))

```

**Step 4: Add Boxplot and Data Points to the Graph**

I am a huge fan of showing each specific data point as well as a summary of the data using a boxplot. To do this, I create a new variables for the graphs called
`p_plots` and `q_plots`, and use `geom_boxplot` and `geom_point` to populate the data in the graph. 
. 
Have set the theme as `theme_bw` with a bunch of modifications to the axis/legend fonts and sizes. These are all done through the `theme()` command using commands like `axis.text.x(element_text=()`. I get rid of grid lines using `theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()`. The axis titles and text can all be configured using commdands such as `axis.text.y=element(size=FONTSIZE, family=FONT_FAMILY)` as well as `axis.title.x=element(size=FONTSIZE, family=FONT_FAMILY`. 

For three of the graphs I use `stat_compare_means` to put significance values in the graph itself

*also important to note that if you are copy/pasting the lines below that the '+' cannot be on a new line or the code will not work*

```
p_plot<-p %>% ggplot(aes(x=Strain_Pair, y=Comp_Fitness))+geom_boxplot(width=0.4, outlier.shape=NA)+geom_point(alpha=0.5,fill="grey",pch=21,position=position_jitter(width=0.11))+theme_bw(base_size=16)+scale_fill_manual(values=c("grey","white"))+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(), strip.text=element_text(size=10,family = "Helvetica"),axis.text.x=element_text(size=13,family = "Helvetica"),legend.position="none",axis.title.x=element_text(size=15,family="Helvetica"))+ylab("Relative Fitness")+stat_compare_means(method="t.test",label="p.signif",label.y=0,label.x=0.75,size=5.,symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns")))+xlab("Fitness Relative to DBL305")

q_plot<-q %>% ggplot(aes(x=Treatment, y=Relative.Fitness))+geom_boxplot(width=0.4, outlier.shape=NA)+geom_point(alpha=0.5,fill="grey",pch=21,position=position_jitter(width=0.11))+theme_bw(base_size=16)+scale_fill_manual(values=c("grey","white"))+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(), strip.text=element_text(size=10,family = "Helvetica"),axis.text.x=element_text(size=12,family = "Helvetica"),axis.title.x=element_text(size=15, family ="Helvetica"),axis.title.y=element_text(size=16,family="Helvetica"),legend.position="none")+ylab("Relative Fitness")+xlab(expression("Naldixic Acid Concentration " * mu*"g/mL"))

r_plot<-r %>% ggplot(aes(x=Treatment, y=Relative.Fitness))+geom_boxplot(width=0.4, outlier.shape=NA)+geom_point(alpha=0.5,fill="grey",pch=21,position=position_jitter(width=0.11))+theme_bw(base_size=16)+scale_fill_manual(values=c("grey","white"))+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(), strip.text=element_text(size=10,family = "Helvetica"),axis.text.x=element_text(size=12,family = "Helvetica"),axis.title.x=element_text(size=16, family ="Helvetica"),axis.title.y=element_text(size=15,family="Helvetica"),legend.position="none")+ylab("Relative Fitness")+xlab(expression("Ciprofloxaxin Concentration " * mu*"g/mL"))+stat_compare_means(method="t.test",label="p.signif",label.y=0,label.x=0.75,size=5.,symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns")))

s_plot<-s %>% ggplot(aes(x=Treatment, y=Relative.Fitness))+geom_boxplot(width=0.4, outlier.shape=NA)+geom_point(alpha=0.5,fill="grey",pch=21,position=position_jitter(width=0.11))+theme_bw(base_size=16)+scale_fill_manual(values=c("grey","white"))+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(), strip.text=element_text(size=10,family = "Helvetica"),axis.text.x=element_text(size=12,family = "Helvetica"),axis.title.x=element_text(size=16, family ="Helvetica"),axis.title.y=element_text(size=15,family="Helvetica"),legend.position="none")+ylab("Relative Fitness")+xlab(expression("Temperature ( "^"o"*"C)"))+stat_compare_means(method="t.test",label="p.signif",label.y=0,label.x=0.75,size=5.,symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns")))
```

**Step 5: Combine Graphs Together into One Grid***

Now I call all the graphs together using `plot_grid`. I can set the label, and use `hjust` and `vjust` to put the label in the correct position

```
pqrs_plot<-plot_grid(p_plot,q_plot,r_plot,s_plot,labels=c("A)","B)","C)","D)"),hjust=0,vjust=1.2,label_size=20)
```

**Step 6: Export Graph to a Figure File**

The last step here will be to export the graph to a readable figure file using the ggsave command. In this case, I will export as both `.png` files and `.svg` files on my desktop and called `Fig1.png` or `Fig1.svg`.

```
ggsave(file="~/Desktop/Fig4.png",plot=pqrs_plot,width=8,height=6)
ggsave(file="~/Desktop/Fig4.svg",plot=pqrs_plot,width=8,height=6)
```


**Step 6: Statistics**

OK, one figure at a time for the statistics. First up is Figure 4A, which compares wild type white/blue megaplasmid- strains as well as white megaplasmid- and blue megaplasmid+ strain in the P. putida background.

Set up the ANOVA framework:

```
Comp_Fitness_Putida<-aov(Comp_Fitness~Strain_Pair+Strain_Pair:Assay, data=Fig4A_data)
```
Check the summary stats of the ANOVA.

```
summary(Comp_Fitness_Putida)
```
Which yields the following stats:
```                         Df Sum Sq Mean Sq F value Pr(>F)    
Tailocin            4 1592.2   398.1   298.6 <2e-16 ***
Tailocin:Replicate  5  557.8   111.6    83.7 <2e-16 ***
Residuals          47   62.6     1.3                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

So there are significant effects for both the treatment (Tailocin) as well as a block assay effect. Now I'll set up a Post-hoc Tukey's HSD test to step through the Nal results
```
TukeyHSD(Tailocin_Xize)
```
which yields (there are values for the assay block effect too, but I'm not going to include those below
```
 Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Area ~ Tailocin + Tailocin:Replicate, data = Fig1B_data)

$Tailocin
                         diff         lwr         upr     p adj
Ptt-A-11-A          9.6198168   8.2234469  11.0161868 0.0000000
Ptt-B-11-A         -0.9512992  -2.3182686   0.4156702 0.2943747
PttPEG-A-11-A      10.3152374   8.9188675  11.7116074 0.0000000
PttPEG-B-11-A      -1.4865542  -2.8535235  -0.1195848 0.0268581
Ptt-B-Ptt-A       -10.5711160 -11.9380854  -9.2041466 0.0000000
PttPEG-A-Ptt-A      0.6954206  -0.7009493   2.0917906 0.6227797
PttPEG-B-Ptt-A    -11.1063710 -12.4733404  -9.7394016 0.0000000
PttPEG-A-Ptt-B     11.2665366   9.8995673  12.6335060 0.0000000
PttPEG-B-Ptt-B     -0.5352550  -1.8721774   0.8016674 0.7869030
PttPEG-B-PttPEG-A -11.8017916 -13.1687610 -10.4348222 0.0000000

```
Which means that the size of the overlay of the Ptt tailocin (and PEG prepped version of this tailocin) against sensitivity class A strains signficantly differs by post hoc test from the tailocin from strain 011 against sensitivity class A strains!

There does not appear to be a significant difference of PEG treatment on the size of the overlay, and the size of the overlay for strain 011 tailocin against sensitivity class A is no statistically different than the size of the Ptt tailocin overlay against sensitivity class B. 

