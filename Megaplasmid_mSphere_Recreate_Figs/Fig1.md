# A Basic Guide for Recreating Figure 1 from This Paper

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

```
Fig1_data<-read.csv(text=getURL("https://raw.githubusercontent.com/surtlab/data_for_figures/master/mSphere_megaplasmid_2020/Final_Stutzeri28a24_Nal_Results.csv"))
```

**Step 3: Create the Plot Background**

This step does a couple of different things at once. First, note the `%>%` function which acts as a pipe between commands. Then I'm using 'mutate' to reorder the data on the x-axis. Then I'm loading the data into ggplot. Them I'm adding the boxplot (adding this before points so that it is behind points. 
Also, want to make sure that outliers are taken away from the boxplot (or they'll show up as redundant once the points are added).'

*also important to note that if you are copy/pasting the lines below that the '+' cannot be on a new line or the code will not work*

```
q<-Fig1_data %>% mutate(Strain=fct_relevel(Strain,"DBL332","DBL386","DBL453","DBL408"))


```

**Step 4: Add Boxplot and Data Points to the Graph**

I am a huge fan of showing each specific data point as well as a summary of the data using a boxplot. To do this, I create a new variable for the graph called `q_plots`, and use `geom_boxplot` and `geom_point` to populate the data in the graph. I've set the color for the boxplot to vary based off of the tratement. 
Have set the theme as `theme_bw` with a bunch of modifications to the axis/legend fonts and sizes. These are all done through the `theme()` command using commands like `axis.text.x(element_text=()`. I get rid of grid lines using `theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()`. I've also set the particular colors for each treatment using the `scale_fill_manual` command and by designating colors for each treatment as they appear after the reordering from above. Fill for this is set in the `aes` line above.
Lastly, I label the Y axis title using `ylab()`.

```
q_plots<-q %>% ggplot(aes(x=Strain, y=Normalized_Area,fill=Treatment))
+geom_boxplot(width=.5,outlier.shape=NA)+geom_point()
+theme_bw()+ scale_fill_manual(values=c("grey","white"),labels=c("+Megaplasmid","-Megaplasmid"),guide=guide_legend(reverse=TRUE))
+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_text(size=12,family = "Helvetica"),axis.text.y=element_text(size=12,family="Helvetica"),legend.position="bottom",axis.title.x=element_blank(),axis.title.y=element_text(size=14,family="Helvetica"),legend.title=element_blank(),legend.text=element_text(size=10,family="Helvetica"))
+ylab("Halo Size Around Naldixic Acid Disc (Normalized)")```
```
**Step 5: Export Graph to a Figure File**

The last step here will be to export the graph to a readable figure file using the ggsave command. In this case, I will export as both `.png` files and `.svg` files on my desktop and called `Fig1.png` or `Fig1.svg`.

```
ggsave(file="~/Desktop/Fig1.png",plot=q_plots,width=8,height=6)
ggsave(file="~/Desktop/Fig1.svg",plot=q_plots,width=8,height=6)
```


**Step 6: Statistics

First thing to do is perform an ANOVA to see if the megaplasmid increases halo size in these strains

```
CompPMP28a24<-aov(Normalized_Area~Treatment+Strain+Treatment:Strain, data=Fig1_data)
summary(CompPMP28a24)
```

Which gives you the following result

```
	   Df Sum Sq Mean Sq F value   Pr(>F)    
Treatment    1 1.6352  1.6352 132.143 1.32e-13 ***
Strain       2 0.0033  0.0016   0.133    0.876    
Residuals   36 0.4455  0.0124                     
---

```
There is a statistically significant difference! The megaplasmid sensitizes P. stutzeri strain 28a24 (two versions of it) to Nalidixic acid according to Kirby-Bauer diffusion assays
