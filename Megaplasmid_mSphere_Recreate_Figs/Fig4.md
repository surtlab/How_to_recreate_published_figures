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
```            Df  Sum Sq Mean Sq F value   Pr(>F)    
Strain_Pair        1 0.21936 0.21936   600.3  < 2e-16 ***
Strain_Pair:Assay  2 0.08519 0.04260   116.6 9.47e-12 ***
Residuals         20 0.00731 0.00037                     
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
Now I test to see whether there is any effect of the blue/white pair of control strains by testing whether this ratio is different than 1. I pull this particular strain pair out of the overall data set.
```
wt_comp<-subset(Fig4A_data, Strain_Pair=="C305-1604",select=c(Comp_Fitness))
```
And test using a ttest against the expectation of "1".
```
t.test(wt_comp,mu=1,var.equal=FALSE)
```
Which yields the answer that, no, there doesn't seem to be an effect of the marker (blue) gene
```
	One Sample t-test
data:  wt_comp
t = 1.5792, df = 11, p-value = 0.1426
alternative hypothesis: true mean is not equal to 1
95 percent confidence interval:
 0.995927 1.024764
sample estimates:
mean of x 
 1.010345 
 ```
Now I simply do the same ANOVA for each of the comparisons in Figures 4B, 4C, and 4D.
```
Nal_comp<-subset(Fig4BCD_data, Experiment=="Nal",select=c(Relative.Fitness,Treatment,Replicate))
Cip_comp<-subset(Fig4BCD_data, Experiment=="Cipro",select=c(Relative.Fitness,Treatment,Replicate))
Temp_comp<-subset(Fig4BCD_data, Experiment=="Temp",select=c(Relative.Fitness,Treatment,Replicate))
```
And set up the ANOVA framework for each experiment.
```
Nal_Fitness<-aov(Relative.Fitness~Treatment+Treatment:Replicate, data=Nal_comp)
Cip_Fitness<-aov(Relative.Fitness~Treatment+Treatment:Replicate, data=Cip_comp)
Temp_Fitness<-aov(Relative.Fitness~Treatment+Treatment:Replicate, data=Temp_comp)
```
and now take the ANOVA summary of each in step, first up the Naldixic acid comparison.
```
summary(Nal_Fitness)
```
which yields
```
                  Df Sum Sq Mean Sq F value   Pr(>F)    
Treatment            2 0.2897 0.14484  264.82 1.54e-13 ***
Treatment:Replicate  3 0.1477 0.04924   90.04 1.23e-10 ***
Residuals           17 0.0093 0.00055                     
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
So there are significant effects for both the treatment (Nal level) as well as a block assay effect. Now I'll set up a Post-hoc Tukey's HSD test to step through the Nal results
```
TukeyHSD(Nal_Fitness)
```
which yields (there are values for the assay block effect too, but I'm not going to include those below
```
  Tukey multiple comparisons of means
    95% family-wise confidence level
Fit: aov(formula = Relative.Fitness ~ Treatment + Treatment:Replicate, data = Nal_comp)
$Treatment
                               diff       lwr       upr p adj
Nal 20 ng/uL-Nal 0 ng/uL  0.1413158 0.1113179 0.1713137     0
Nal 40 ng/uL-Nal 0 ng/uL  0.2782291 0.2471783 0.3092798     0
Nal 40 ng/uL-Nal 20 ng/uL 0.1369133 0.1058626 0.1679641     0
```
Which means that all of the comparisons are signficantly different by post hoc test!

Now, I can simply do the same kind of ANOVA for each of the remaining two datasets (although I don't have to really perform a Tukey's HSD on these because there are only two levels in each case so the significance terms will pop out of the summary stats as a significant `Treatment` effect.
```
summary(Cip_Fitness)
```
Gives a significant treatment result.
```
                    Df Sum Sq Mean Sq F value   Pr(>F)    
Treatment            1 1.8124  1.8124  178.89 8.64e-11 ***
Treatment:Replicate  4 0.0482  0.0121    1.19    0.349    
Residuals           18 0.1824  0.0101                     
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
And the same for the temperature effect.
```
summary(Temp_Fitness)
```
```
                  Df Sum Sq Mean Sq F value   Pr(>F)    
Treatment            1 0.3228  0.3228 319.380 6.66e-13 ***
Treatment:Replicate  4 0.0396  0.0099   9.799 0.000217 ***
Residuals           18 0.0182  0.0010                     
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

