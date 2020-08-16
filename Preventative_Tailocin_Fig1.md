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
Fig1_data<-read.csv(text=getURL("https://raw.githubusercontent.com/surtlab/data_for_figures/master/Cit7_final.csv"))
```

**Step 4: Create the Plot Background**

This step does a couple of different things at once. First, note the `%>%` function which acts as a pipe between commands. Next, note that I am reordering the factors on the x-axis to be shown in the order I want. You can alter the orders here in case you'd like to rearrange how the data is presented. Next I am labelling the axes and getting rid of the legend. All of these commands create a graph which I've called `q`.

*also important to note that if you are copy/pasting the lines below that the '+' cannot be on a new line or the code will not work*

```
q<-Fig1_data%>%
mutate(Strain=fct_relevel(Strain,"PsyCit7","PsyB728a","PaePA01","No Tailocin"))%>%
ggplot(aes(x=Strain,y=LogCFU,color=Strain))
+theme(legend.position="none",axis.title=element_text(size=24),panel.grid=element_blank())+labs(x="Strain", y="Log10 Colony Forming Units (CFU)")
```

**Step 5: Add Boxplot and Data Points to the Graph**

I am a huge fan of showing each specific data point as well as a summary of the data using a boxplot. To do this, I create a new variable for the graph called `qbp_jitter_color`, and use `geom_boxplot` and `geom_jitter` to populate the data in the graph. I've set the color for the boxplot as a somewhat dark grey, and I've set the alpha value for geom_jitter (which determines the transparency of the data points). Within this command, I've also set the particular colors for each of the strain's data points through the `scale_color_manual` command and by designating colors for each strain as they appear after the reordering from above.

```
qbp_jitter_color<-q
+geom_boxplot(color = "gray50", outlier.alpha = 0)
+geom_jitter(size = 2, alpha = 0.4, width = 0.2)
+scale_color_manual(values=c("blue","springgreen4","mediumpurple","grey38"))
```

**Step 6: Export Graph to a Figure File**

The last step here will be to export the graph to a readable figure file using the ggsave command. In this case, I will export as both `.png` files and `.svg` files on my desktop and called `Fig1.png` or `Fig1.svg`.

```
ggsave(file="~/Desktop/Fig1.png",plot=qbp_jitter_color,width=10,height=8)
ggsave(file="~/Desktop/Fig1.svg",plot=qbp_jitter_color,width=10,height=8)
```


**Step 7: Statistics

First thing to do is perform a non-parametric Kruskal-Wallis test to see if any strain is different than any other strain in terms of prevention

```
kruskal.test(LogCFU~Strain, data = Fig1_data)
```

Which gives you the following result

```
		Kruskal-Wallis rank sum test

data:  LogCFU by Strain
Kruskal-Wallis chi-squared = 34.812, df = 3, p-value =
1.335e-07
```
There is a statistically significant difference! That justifies a follow up test, a pairwise Wilcox test with correction for multiple testing. 

```
 pairwise.wilcox.test(Fig1_data$LogCFU, Fig1_data$Strain,
                 p.adjust.method = "BH")

```
Which gives you to following result

```
	Pairwise comparisons using Wilcoxon rank sum test with continuity correction 

data:  Fig1_data$LogCFU and Fig1_data$Strain 

         No Tailocin PaePA01 PsyB728a
PaePA01  0.53        -       -       
PsyB728a 0.18        0.69    -       
PsyCit7  5.2e-06     5.2e-06 5.2e-06 
```
Which means that PsyCit7 is grouped into A, and the other treatments are significantly different than group A and are therefore grouped into group B
