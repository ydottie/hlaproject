
rm(list = ls())

library(tidyverse)
library(data.table)

raw_data_memory <- fread("work_memory.csv")[,-1]
colnames(raw_data_memory)
colnames(raw_data_memory)<-c("arcasHLA","hlaforest","phlat","seq2HLA", "hlaminer (alignment)", "hlaminer (assembly)", "hlahd","hlavbseq","optitype")
colnames(raw_data_memory)

data_memory <- raw_data_memory %>% as_tibble %>% 
  gather(Tools, rawMemory) %>%
  mutate(Memory = rawMemory / 10^3)

mem <- ggplot(data_memory, aes(x=reorder(Tools,Memory,FUN=median), y=Memory, fill=Tools))+
  theme_bw()+
  geom_boxplot()+
  coord_trans(y="log10")+
  scale_fill_brewer(palette="Set3") +
  ggtitle("Memory used by the process")+
  labs(x="",y="Max Memory (MB)")+
  theme(axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.justification=c(1,0),legend.position=c(1, 0),legend.title=element_blank(),
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
mem
#####################################################################
#####################################################################
#####################################################################

raw_data_time <- fread("work_time.csv")[,-1]
colnames(raw_data_time)
colnames(raw_data_time)<-c("arcasHLA","hlaforest","phlat","seq2HLA", "hlaminer (alignment)", "hlaminer (assembly)", "hlahd","hlavbseq","optitype")
colnames(raw_data_time)

data_time <- raw_data_time %>% as_tibble %>% 
  gather(Tools, time) 

core <- ggplot(data_time, aes(x=reorder(Tools, time, FUN=median), y=time, fill=Tools))+
  theme_bw()+
  geom_boxplot()+
  coord_trans(y="log10")+
  scale_fill_brewer(palette="Paired") +
  ggtitle("CPU time used by the process")+
  labs(x="",y="CPU Time('s)")+
  theme(axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
core

cowplot::plot_grid(core, mem, align = "hv") +
  ggsave(file="work_memcore.png", width = 10, height = 6, units = "in")
