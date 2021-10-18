load("../data/source/Svenson_DO_HFD_v10_03.26.2021.RData", svenson <- new.env()) # Gene expression
load("../data/source/do_proteomics_qtlviewer.Rdata", svenson) # Protein Expression
ls(svenson)

library(rstatix)
library(tidyverse)
library(bmediatR)
library(doFuture)
registerDoFuture()
plan(multisession)

#################################################################################
## Find diet and sex effects
#################################################################################

diet_phenotypes <- svenson$dataset.phenotype$data %>%
  as.data.frame() %>% rownames_to_column("mouse.id") %>%
  pivot_longer(!mouse.id, names_to = "phenotype", values_to = "value") %>%
  merge(svenson$dataset.phenotype$covar.matrix %>%
          as.data.frame %>% rownames_to_column("mouse.id") %>%
          select(mouse.id, sexM, diethf)) %>%
  group_by(phenotype) %>%
  na.omit() %>%
  anova_test(value ~ diethf, type = 1) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  add_significance()

diet_rna <- svenson$dataset.mrna$data$rz[,1:5] %>%
  as.data.frame() %>% rownames_to_column("mouse.id") %>%
  pivot_longer(!mouse.id, names_to = "gene.id", values_to = "value") %>%
  merge(svenson$dataset.phenotype$covar.matrix %>%
          as.data.frame %>% rownames_to_column("mouse.id") %>%
          select(mouse.id, sexM, diethf)) %>%
  plyr::ddply(~gene.id, function(df) df %>% na.omit() %>% anova_test(value ~ diethf, type = 1),
              .parallel = T) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  add_significance()

sex_phenotypes <- svenson$dataset.phenotype$data %>%
  as.data.frame() %>% rownames_to_column("mouse.id") %>%
  pivot_longer(!mouse.id, names_to = "phenotype", values_to = "value") %>%
  merge(svenson$dataset.phenotype$covar.matrix %>%
          as.data.frame %>% rownames_to_column("mouse.id") %>%
          select(mouse.id, sexM, diethf)) %>%
  group_by(phenotype) %>%
  na.omit() %>%
  anova_test(value ~ diethf, type = 1) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  add_significance()

sex_rna <- svenson$dataset.mrna$data$rz %>%
  as.data.frame() %>% rownames_to_column("mouse.id") %>%
  pivot_longer(!mouse.id, names_to = "gene.id", values_to = "value") %>%
  merge(svenson$dataset.phenotype$covar.matrix %>%
          as.data.frame %>% rownames_to_column("mouse.id") %>%
          select(mouse.id, sexM, diethf)) %>%
  plyr::ddply(~gene.id, function(df) df %>% na.omit() %>% anova_test(value ~ diethf, type = 1),
              .parallel = T) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  add_significance()

#################################################################################
## Mediate diet effects on phenotypes with trsnacript data
#################################################################################

phenotypes <- diet_phenotypes %>%
  as_tibble() %>%
  filter(p.adj.signif != "ns") %>%
  pull(phenotype)

dataset_mrna <- svenson$dataset.mrna$data$rz
transcript_mediation <- svenson$dataset.phenotype$data[,1:5, drop = F] %>%
  as.data.frame %>% rownames_to_column("mouse.id") %>%
  merge(svenson$dataset.phenotype$covar.matrix[,"diethf",drop = F] %>% as.data.frame() %>% rownames_to_column("mouse.id")) %>%
  pivot_longer(!c(mouse.id, diethf), names_to = "phenotype") %>%
  plyr::dlply(~phenotype, function(df) {
    y = matrix(df$value, ncol = 1, dimnames = list(df$mouse.id, unique(df$phenotype)))
    M = dataset_mrna
    X = matrix(df$diethf, ncol = 1, dimnames = list(df$mouse.id, "diethf"))
    mice = Reduce(intersect, list(rownames(y), rownames(M), rownames(X)))
    bmediatR(y[mice,,drop = F],
             M[mice,,drop = F],
             X[mice,,drop = F],
             options_X = list(scale = T, center = T, sum_to_zero = F)
    )
  },
  .parallel = T)
saveRDS(transcript_mediation, file = "~/Desktop/svenson_diet_mediation.rds")

transcript_mediation[["bmd1"]]$ln_post_odds %>%  as.data.frame() %>%
  rownames_to_column("gene.id") %>% merge(svenson$dataset.mrna$annot.mrna) %>%
  ggplot(aes(middle, complete)) + geom_point() +
  facet_wrap(~chr, nrow = 1, strip.position = "bottom") +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank())  +
  ggrepel::geom_text_repel(data = . %>% filter(complete > 0), aes(label = symbol))
