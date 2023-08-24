# A Basic Guide for Recreating Figure 1 from This Paper

**Step 1: Prepare the Environment in R**

These commands will load various libraries into the active R environment. If you receive an error that one of these libraries isn't recognized, simply use the command install.packages()`

```
library(RCurl)
library(forcats)
library(dplyr)
library(showtext)
library(ggplot2)
library(svglite)
```

**Step 2: Set Basic Formats For Figure**

This will load the theme for the graph as well as fonts for the axes. In particular the `font_add_google` function will This function will search the [Google Fonts Repository](https://fonts.google.com/) for a specified family name, download the proper font files, and then add them to sysfonts.

```
font_add_google("Poppins", "Poppins")
font_add_google("Roboto Mono", "Roboto Mono")
showtext_auto()
theme_set(theme_light(base_size=24, base_family="Poppins"))
```

**Step 3: Acquire Data to Recreate Figure 1 from Github**

The raw data underlying Figure 1 is openly available as a .csv file in Github. The RCurl function allows you to import this data from the webpage, however, you need to make sure that the address is correct and you are pulling in in the raw data.
We'll pull this raw data into a variable called `Fig1_data` through the use of the `read.csv` function.

```
Fig1B_data<-read.csv(text=getURL("https://raw.githubusercontent.com/surtlab/data_for_figures/master/LLB_Fig1B_data.csv"))
```


**Step 4: Create the Plot Background**

This step does a couple of different things at once. First, note the `%>%` function which acts as a pipe between commands. Next, note that I am reordering the factors on the x-axis to be shown in the order I want. You can alter the orders here in case you'd like to rearrange how the data is presented. Next I am labelling the axes and getting rid of the legend. All of these commands create a graph which I've called `q`.

*also important to note that if you are copy/pasting the lines below that the '+' cannot be on a new line or the code will not work*

```
q<-Fig1B_data%>%
mutate(Tailocin=fct_relevel(Tailocin,"Ptt-A","PttPEG-A","11-A","Ptt-B","PttPEG-B"))%>%
ggplot(aes(x=Tailocin,y=Area,color=Tailocin))
+theme(legend.position="none",axis.title=element_text(size=24),panel.grid=element_blank())+labs(x="Tailocin", y="Area(mm2)")
```

**Step 5: Add Boxplot and Data Points to the Graph**

I am a huge fan of showing each specific data point as well as a summary of the data using a boxplot. To do this, I create a new variable for the graph called `qbp_jitter_color`, and use `geom_boxplot` and `geom_jitter` to populate the data in the graph. I've set the color for the boxplot as a somewhat dark grey, and I've set the alpha value for geom_jitter (which determines the transparency of the data points). Within this command, I've also set the particular colors for each of the strain's data points through the `scale_color_manual` command and by designating colors for each strain as they appear after the reordering from above.

```
qbp_jitter_color<-q
+geom_boxplot(color = "gray50", outlier.alpha = 0)
+geom_jitter(aes(shape=Replicate),size = 2, alpha = 0.4, width = 0.2)
+scale_color_manual(values=c("blue","springgreen4","mediumpurple","grey38","red"))
```

**Step 6: Export Graph to a Figure File**

The last step here will be to export the graph to a readable figure file using the ggsave command. In this case, I will export as both `.png` files and `.svg` files on my desktop and called `Fig1.png` or `Fig1.svg`.

```
 ggsave(file="~/Desktop/Fig1B.png",plot=qbp_jitter_color,width=10,height=8)
 
 ggsave(file="~/Desktop/Fig1B.svg",plot=qbp_jitter_color,width=10,height=8)
```
**Step 6: Statistics**

OK, one figure at a time for the statistics. First up is Figure 4A, which compares wild type white/blue megaplasmid- strains as well as white megaplasmid- and blue megaplasmid+ strain in the P. putida background.

Set up the ANOVA framework:

```
Tailocin_size<-aov(Area~Tailocin+Tailocin:Replicate, data=Fig1B_data)
```
Check the summary stats of the ANOVA.

```
summary(Tailocin_size)
```
Which yields the following stats:
```            Df  Sum Sq Mean Sq F value   Pr(>F)    
Strain_Pair        1 0.21936 0.21936   600.3  < 2e-16 ***
Strain_Pair:Assay  2 0.08519 0.04260   116.6 9.47e-12 ***
Residuals         20 0.00731 0.00037                     
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
