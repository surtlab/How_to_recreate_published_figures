**Step 1: Prepare the Environment in R**

These commands will load various libraries into the active R environment. If you receive an error that one of these libraries isn't recognized, simply use the command install.packages()`

```
library(RCurl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
```

**Step 2: Acquire Data to Recreate Figure 1 from Github**

For this one we are going to build two different graphs and then combine them together at the end. So first, we need to load data for the two different graphs:

```
Fig3A_data<-read.csv(text=getURL("https://raw.githubusercontent.com/surtlab/data_for_figures/master/Final_Putida_Nal_Results.csv"))
Fig3B_data<-read.csv(text=getURL("https://raw.githubusercontent.com/surtlab/data_for_figures/master/Final_Putida_Cip_Results.csv"))
```
**Step 3: Create the Plot Data Order**

This step does a couple of different things at once. First, note the `%>%` function which acts as a pipe between commands. Then I'm using 'mutate' to reorder the data on the x-axis. I do this for both data sets for consistency, and call the data for each graph into a different variable, `p` or `q`

```
p<-Fig3A<-Fig3A_data %>% mutate(Treatment=fct_relevel(Treatment,"wt","pMP"))%>% mutate(Treatment=recode(Treatment,'wt'="-pMP",'pMP'="+pMP"))%>% mutate(Strain=fct_relevel(Strain,"DBL305","DBL1604","DBL759","DBL1620"))
q<-Fig3B<-Fig3B_data %>% mutate(Treatment=fct_relevel(Treatment,"wt","pMP"))%>% mutate(Treatment=recode(Treatment,'wt'="-pMP",'pMP'="+pMP"))%>% mutate(Treatment=recode(Treatment,'wt'="-pMP",'pMP'="+pMP"))%>%mutate(Strain=fct_relevel(Strain,"DBL305","DBL1604","DBL759","DBL1620"))

```

**Step 4: Add Boxplot and Data Points to the Graph**

I am a huge fan of showing each specific data point as well as a summary of the data using a boxplot. To do this, I create a new variables for the graphs called
`p_plots` and `q_plots`, and use `geom_boxplot` and `geom_point` to populate the data in the graph. 
. 
Have set the theme as `theme_bw` with a bunch of modifications to the axis/legend fonts and sizes. These are all done through the `theme()` command using commands like `axis.text.x(element_text=()`. I get rid of grid lines using `theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()`. The lengend and axis titles can all be configured using commdands such as `axis.text.y=element(size=FONTSIZE, family=FONT_FAMILY)` as well as `legend.text=element_text()`.

I've also set the particular colors for each treatment using the `scale_fill_manual` command and by designating colors for each treatment as they appear after the reordering from above. Fill for this is set in the `aes` line above.

I also place the legend inside of the first figure using `legend.position=c(.2,.1)`

*also important to note that if you are copy/pasting the lines below that the '+' cannot be on a new line or the code will not work*

```
p_plots<-p%>%ggplot(aes(x=Strain, y=Normalized_Area,fill=Treatment))+geom_boxplot(width=0.4, outlier.shape=NA)
+geom_point(alpha=0.5,fill="grey",pch=21,position=position_jitter(width=0.11))
+theme_bw(base_size=14)+scale_fill_manual(values=c("white","grey"))
+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_text(size=12,family = "Helvetica"),axis.text.y=element_text(size=12,family="Helvetica"),legend.position=c(.2,.1),axis.title.x=element_blank(),axis.title.y=element_text(size=14,family="Helvetica"),legend.title=element_blank(),legend.text=element_text(size=12,family="Helvetica"))
+ylab("Halo Size Around Naldixic Acid Disc (Normalized)")
```
I do this for both sets of data
```
q_plots<-q%>%ggplot(aes(x=Strain, y=Normalized_Area,fill=Treatment))
+geom_boxplot(width=0.4, outlier.shape=NA)+geom_point(alpha=0.5,fill="grey",pch=21,position=position_jitter(width=0.11))
+scale_fill_manual(values=c("white","grey"))
+theme_bw()
+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_text(size=12,family = "Helvetica"),axis.text.y=element_text(size=12,family="Helvetica"),legend.position="none",axis.title.x=element_blank(),axis.title.y=element_text(size=14,family="Helvetica"),legend.title=element_blank(),legend.text=element_text(size=10,family="Helvetica"))
+ylab("Halo Size Around Ciprofloxacin Disc (Normalized)")
```

**Step 5: Combine Graphs Together into One Grid***

```
pq_plot<-plot_grid(p_plots,q_plots,labels=c("A)","B)"),label_size=20)
```

**Step 6: Export Graph to a Figure File**

The last step here will be to export the graph to a readable figure file using the ggsave command. In this case, I will export as both `.png` files and `.svg` files on my desktop and called `Fig1.png` or `Fig1.svg`.

```
ggsave(file="~/Desktop/Fig3.png",plot=pq_plot,width=8,height=6)
ggsave(file="~/Desktop/Fig3.svg",plot=pq_plot,width=8,height=6)
```


**Step 7: Statistics**

Going to show the statistics for each graph separately below. First up is performing a ttest on the comparisons from Figure 3A, the Naldixic acid disc overlays for different P. putida strain pairs. First thing I will do is pull the data from each strain into it's own set from the overall dataset:

```
DBL305_Nal<-subset(Fig3A_data, Strain=="DBL305",select=c(Normalized_Area))
DBL759_Nal<-subset(Fig3A_data, Strain=="DBL759",select=c(Normalized_Area))
DBL1604_Nal<-subset(Fig3A_data, Strain=="DBL1604",select=c(Normalized_Area))
DBL1620_Nal<-subset(Fig3A_data, Strain=="DBL1620",select=c(Normalized_Area))
```

Next, we can simply call the `t.test` function to compare pairs of strains. Don't necessarily need to do this in a full ANOVA context because the data within easy assay is normalized to that assay (which limits some of the assay to assay variance).

```
t.test(DBL305_Nal,DBL759_Nal, var.equal=FALSE)
```

Which gives the following results for strain pair DBL305 and DBL759:

```
	Welch Two Sample t-test

data:  DBL305_Nal and DBL759_Nal
t = 4.1626, df = 32.132, p-value = 0.00022
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.06903193 0.20129320
sample estimates:
mean of x mean of y 
1.0000198 0.8648572 
```

and then we do the same for strain pair DBL1604 and DBL1620:

```
t.test(DBL1604_Nal,DBL759_Nal,var.equal=FALSE)
```

giving the following results for strain pair DBL1604 and DBL1620:

```
	Welch Two Sample t-test

data:  DBL1604_Nal and DBL759_Nal
t = 5.6875, df = 37.037, p-value = 1.659e-06
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.08701157 0.18331103
sample estimates:
mean of x mean of y 
1.0000185 0.8648572 
```

Next, I repeat the exact same analysis for the Ciprofloxacin pairs from the `Fig3B_data` dataset.

```
DBL305_Cip<-subset(Fig3B_data, Strain=="DBL305",select=c(Normalized_Area))
DBL759_Cip<-subset(Fig3B_data, Strain=="DBL759",select=c(Normalized_Area))
DBL1604_Cip<-subset(Fig3B_data, Strain=="DBL1604",select=c(Normalized_Area))
DBL1620_Cip<-subset(Fig3B_data, Strain=="DBL1620",select=c(Normalized_Area))
```

```
t.test(DBL305_Cip,DBL759_Cip, var.equal=FALSE)
```

Which gives the following results for strain pair DBL305 and DBL759:

```
	Welch Two Sample t-test

data:  DBL305_Cip and DBL759_Cip
t = -3.6281, df = 35.75, p-value = 0.000884
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.14264362 -0.04033493
sample estimates:
mean of x mean of y 
 1.000010  1.091499 
```

and then we do the same for strain pair DBL1604 and DBL1620:

```
t.test(DBL1604_Cip,DBL1620_Cip, var.equal=FALSE)
```

giving the following results for strain pair DBL1604 and DBL1620:

```
		Welch Two Sample t-test

data:  DBL1604_Cip and DBL1620_Cip
t = -1.0747, df = 26.284, p-value = 0.2923
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.07118778  0.02228794
sample estimates:
mean of x mean of y 
 1.000009  1.024459 
```
