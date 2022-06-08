
rm(list = ls())
library(tidyverse)
library(data.table)
source("method_hla.R")
se <- function(x) sd(x) / sqrt(length(x))

files <- list.files("ERP000101_summary_hla/")
raw_res <-lapply(1:length(files), function(x){
  dat <- fread(paste0("ERP000101_summary_hla/", files[x])) %>% as_tibble %>% select(1:7)
  cat(files[x], ": ")
  data.frame(fourdigit=mean(cal_acc(dat,gsfile="gold_standard.csv")$fourdigit),
             twodigit=mean(cal_acc(dat,gsfile="gold_standard.csv",type = 2)$fourdigit),
             se_fourdigit=se(cal_acc(dat,gsfile="gold_standard.csv")$fourdigit),
             se_twodigit=se(cal_acc(dat,gsfile="gold_standard.csv",type = 2)$fourdigit))
})

res <- data.frame(names = sub("(.*).csv","\\1",files), rbindlist(raw_res))
res %>% arrange(desc(fourdigit))

#####################################################################
#####################################################################
#####################################################################

tools <- c("optitype", "seq2HLA", "hlaforest_JJ", "phlat",
           "arcasHLA_300_caucasian","arcasHLA_3100_caucasian", "arcasHLA_3200_caucasian", "arcasHLA_3300_caucasian", "arcasHLA_3400_caucasian",
           "hlahd_300", "hlahd_3100", "hlahd_3200", "hlahd_3300", "hlahd_3400",
           "hlaminer_alignment_300_default", "hlaminer_alignment_3100_default", "hlaminer_alignment_3200_default", "hlaminer_alignment_3300_default", "hlaminer_alignment_3400_default",
           "hlaminer_assembly_300_100_old", "hlaminer_assembly_3100_100_old", "hlaminer_assembly_3200_100_old", "hlaminer_assembly_3300_100_old", "hlaminer_assembly_3400_100_old",
           'hlavbseq_300', 'hlavbseq_3100', "hlavbseq_3200", "hlavbseq_3300", "hlavbseq_3400")

res <- res[match(tools,res$names),] %>% 
  mutate(version = c(rep("The reference packaged\nwith the software",4),
                     rep(c("3.0.0","3.10.0","3.20.0","3.30.0","3.40.0"), 5))) %>% 
  mutate(fac_tool=c("optitype", "seq2HLA", "hlaforest", "phlat", 
                    rep("arcasHLA",5), rep("hlahd",5),
                    rep("hlaminer (alignment)",5), rep("hlaminer (assembly)",5),
                    rep("hlavbseq",5))) %>% 
  mutate(fac_tool=factor(fac_tool, levels = c("optitype", "seq2HLA", "arcasHLA", "phlat", "hlaforest", "hlavbseq", "hlahd","hlaminer (assembly)", "hlaminer (alignment)")))

gg_res_mean <- res %>% 
  select(names, fac_tool, version, fourdigit, twodigit) %>% 
  gather("digit","acc",fourdigit:twodigit) %>% 
  mutate(digit=factor(digit,
                      level = c("twodigit", "fourdigit"),
                      labels = c("Two-digit resolution", "Four-digit resolution")))

gg_res_se <- res %>% 
  select(names, fac_tool, version, se_fourdigit, se_twodigit) %>% 
  gather("se_digit","se",se_fourdigit:se_twodigit) %>% 
  select(names, se_digit, se)

gg <- cbind(gg_res_mean,gg_res_se)
all.equal(gg[,1],gg[,6])
gg <- gg[,-6]


ggplot(gg, aes(x=fac_tool, y = acc, fill = version)) +
  geom_hline(yintercept = 0.5, color="grey60",linetype = "dashed")+
  geom_bar(stat = "identity", position=position_dodge(), color = "black") +
  geom_errorbar(aes(ymin=acc - se, ymax=acc+ se), width=0, position=position_dodge(.9)) +
  xlab("") + ylab("Accuracy") + labs(fill = "IPD-IMGT/HLA\nversion") +
  # scale_fill_brewer(palette="Set3")+
  scale_fill_manual(values = c("orange", "dark green", "purple","firebrick3","royalblue3", "grey60"))+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.1)) +
  theme(axis.text.x = element_text(size=12)) +
  facet_wrap(digit~., ncol = 2)+
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  ggsave("performanc.png",width = 9, height = 5, units = "in")


