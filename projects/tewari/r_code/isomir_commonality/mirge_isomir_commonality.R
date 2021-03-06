setwd(here::here())
library(tidyverse)
library(ggplot2)
library(cowplot)
library(edgeR)
library(matrixStats)
source("r_code/functions.R")

theme_set(theme_bw(base_size = 11))

meta_pilot = read_csv("meta_pilot.csv")

complete = read_tsv("tools/mirge/mirtop/expression_counts.tsv.gz")
keys = read_csv("tools/mirge/sample_fn_key.txt", col_names = c("fn", "inside")) %>%
    mutate(fn = gsub("_isomiRs.gff", "", fn))

dds = DGEList(complete[, 13:ncol(complete)])
dds = calcNormFactors(dds)
counts = cpm(dds, normalized.lib.sizes = TRUE)

dds$samples %>% rownames_to_column("sample") %>%
    inner_join(keys, by = c("sample" = "inside")) %>%
    filter(fn  %in%  meta_pilot[["fixed_name"]]) %>%
    left_join(meta_pilot, by = c("fn" = "fixed_name")) %>%
    ggplot(aes(x = replicate, y = lib.size)) +
    geom_bar(stat = "identity") +
    facet_grid(lab~lib_method_simple) +
    ggsave("figures/replicates/mirge_libsize.png", width = 7, height = 9)

new_average<- cbind(counts, average_val = rowMeans2(counts))
probs <- c(0.25, 0.5, 0.75)
pre_pilot = cbind(complete[, 1], new_average)

counts_avg <- data.matrix(new_average[, "average_val"])
cQuans<-colQuantiles(counts_avg, probs=probs)
minQ <- cQuans[[1]]
maxQ <- cQuans[[3]]
all_exprn = cbind(complete[, 1], counts)
low_exprn<-subset(pre_pilot, pre_pilot[,"average_val"] <= minQ)
medium_exprn<-subset(pre_pilot, pre_pilot[,"average_val"] >minQ & pre_pilot[,"average_val"] < maxQ)
high_exprn<-subset(pre_pilot, pre_pilot[,"average_val"] >= maxQ)

exprn = function(fn_exrp, fname)
{
  filename3 = paste("figures/replicates/mirge_counts_per_isomir_type_",fname,".jpg",sep="")
  filename2 = paste("figures/replicates/mirge_",fname,".jpg",sep="")  

  pilot = fn_exrp %>%
    gather(sample, value, -UID) %>%
    inner_join(keys, by = c("sample" = "inside")) %>%
    filter(fn  %in%  meta_pilot[["fixed_name"]]) %>%
    left_join(complete[,1:12]) %>%
    filter(value >= 1) %>%
    mutate(Variant = ifelse(is.na(Variant), "Reference", Variant)) %>%
    left_join(meta_pilot, by = c("fn" = "fixed_name")) %>%
    filter(abs(iso_5p)<4, abs(iso_3p)<4, abs(iso_add)<4 )
  
  lapply(1:3, function(x){
    filter(pilot, value >= x) %>%
        summarize_isomir %>%
        mutate(min_counts = x)
  }) %>% bind_rows() %>%
    plot_summarize_isomir +
    ggsave(filename2, width = 9, height = 9)


  pilot %>% expression_isomirs_by_lab_protocol_isomir %>%
    ggplot(aes(x=lab,y=counts,fill=as.factor(reps))) +
    geom_boxplot() + scale_y_log10() +
    facet_grid(lib_method_simple~isomir_type) +
    ggsave(filename3,width = 9, height = 9)
  rm(pilot)
}

exprn(all_exprn, "all_exprn")
exprn(high_exprn, "high_exprn")
exprn(medium_exprn, "medium_exprn")
exprn(low_exprn, "low_exprn")