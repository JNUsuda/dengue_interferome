# http://www.sthda.com/english/wiki/ggplot2-box-plot-quick-start-guide-r-software-and-data-visualization
# http://sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/

library(ggplot2)
library(tidyr)
library(ggpubr)
library(svglite)

tabela1 = read.delim("../Metanalise/Metanalysis2_R/data/GSE51808_Systems/GSE51808_biomarker_boxplot.txt", 
                     header = T)
tabela1$Phases <- as.factor(tabela1$Phases)
tabela1$Severity <- as.factor(tabela1$Severity)


tabela2 = read.delim("../Metanalise/Metanalysis2_R/data/GSE51808_Systems/GSE51808_biomarker_boxplot2.txt", 
                     header = T)
tabela3 <- tidyr::pivot_longer(data = tabela2, 
                               cols = c("IFI27", "ISG15", "CYBRD1"),
                               names_to = "Gene", values_to = "Expression" )

tabela3$Classification <- factor(tabela3$Classification, levels = c("Acute_DF","Acute_DHF", "Convalescent", "Control"))
tabela3$Gene = factor(tabela3$Gene, levels = c("IFI27", "ISG15", "CYBRD1"))

#orderGene <- c("IFI27", "ISG15", "CYBRD1")
#tabela3 <- transform(tabela3, Gene=factor(Gene,levels=orderGene))

my_comparisons <- list( c("Acute_DF", "Convalescent"), c("Acute_DHF", "Convalescent"))

# Basic box plot
ggplot(tabela3, aes(x=Classification, y=Expression, color = Phases, fill = Phases)) + 
  geom_boxplot() + 
  #scale_x_discrete(limits=c("Acute_DF","Acute_DHF", "Convalescent", "Control")) + 
  geom_jitter(shape=16, position=position_jitter(0.15)) +
  scale_color_manual(values=c("#6a329f", "grey50", "#ff6f00")) +
  scale_fill_manual(values=c("#ab9dc5", "grey80", "#ffb77f")) +
  #ggtitle("IFI27") +
  labs(x = element_blank(), y = "Normalized Expression") +
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        panel.background = element_rect( fill = "white", color  = "grey"),
        panel.grid = element_blank(),
        strip.text.y = element_text(vjust = 1, size = 8), 
        strip.background = element_rect(colour="grey", fill="grey90")) +
  facet_wrap(vars(Gene), scales = "free" ) +
  #stat_compare_means(method = "anova")
  #stat_compare_means(method = "anova", label.y = 40) +   
  stat_compare_means(aes(group = Classification), label = "p.signif", 
                     method = "wilcox.test", comparisons = my_comparisons, 
                     bracket.size = 0.5, step.increase = 0.08, size = 5, 
                     hide.ns=T,  vjust = 0.5)
  
ggsave(filename = "markers_boxplot.svg", path = "./output/figures", width = 2800, height = 1600, 
       units = "px", dpi = 300, device = "svg", scale = 1)
