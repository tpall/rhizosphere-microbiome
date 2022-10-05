library("tidyverse")
library("ggvenn")
library("ggVennDiagram")
library("patchwork")
library("phyloseq")
library("here")
library("cowplot")
library("centroidr")
library("microViz")
library("brms")
library("tidybayes")
library("microbiome")
library("phyloseqGraphTest")
library("ggsci")
library("ggvenn")
library("modelr")

old <- theme_set(theme_cowplot())


#'
#' ## Import dataset
#'
#' Detatiled mothur tax.summary output.
#' 
#+
biom1 <- here("results/Galaxy427-[Make.biom_on_data_421_and_data_418__biom_files__0.03].biom1")
if (!file.exists(biom1)) {
  if (!dir.exists(dirname(biom1))) {
    dir.create(dirname(biom1))
  }
  url <- "https://galaxy.hpc.ut.ee/datasets/273afcfe5a96c3df/display?to_ext=data&hdca_id=c8699e897d691a3e&element_identifier=0.03"
  download.file(
    url = url, 
    destfile = biom1)
}

ps <- import_biom(biom1)
ps <- phyloseq_validate(ps, remove_undetected = TRUE, verbose = FALSE)
ps <- tax_fix(ps, unknowns = c("Actinobacteria_unclassified"))
ps <- subset_samples(ps, str_detect(SAMPLE, "Z\\d"))

#' Drop missing taxa
ps <- ps %>% 
  ps_filter()

#' Update sample_table
ps <- ps %>% 
  ps_mutate(
    sample_name = str_remove(SAMPLE, "UH-3046-") %>% 
      str_split("_") %>% 
      map(1) %>% 
      unlist(),
    zone = str_extract(sample_name, "^Z\\d"),
    rhizo = as.numeric(str_detect(SAMPLE, "plus")),
    zoneplus = paste0(zone, ifelse(as.logical(rhizo), "+", ""))
   )


#' Update tax_table
colnames(tax_table(ps)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

#' 
#' Have look at the summary of our dataset
#'
#+
microbiome::summarize_phyloseq(ps)

#'
#' OTU table
#+  
otu_tab <- microbiome::abundances(ps)

#'
#' Head of the OTU table
#+
otu_tab[1:5,1:5] 

#' 
#' Taxonomy table
#+
tax_tab <- phyloseq::tax_table(ps)

#' 
#' Head of the taxonomy table
#+ 
tax_tab[1:5,1:5]

#'
#' Drop singleton OTUs
#+
pss <- subset_taxa(ps, rowSums(otu_table(ps)) > 0) # rowSums(otu_table(ps)) != 1 # apply(otu_table(ps) / colSums(otu_table(ps)), 1, max) < 0.001
psr <- rarefy_even_depth(ps, rngseed = 11)

#'
#' Top10 tax overview.
#+
ps %>% 
  comp_barplot(tax_level = "Genus", sample_order = "default", n_taxa = 10, label = "sample_name") +
  coord_flip()

#'
#' Graph-based permutation test
#+
gt <- graph_perm_test(ps, "zoneplus", grouping = "sample_name", distance = "jaccard", type = "knn")
plotNet1 <- plot_test_network(gt) + 
  ggsci::scale_color_d3() +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))
plotPerm1 <- plot_permutations(gt)
plotNet1 + plotPerm1

#' 
#' Number of OTUs in *S. alba* cuttings with and without *S. europea* rhizosphere 
#' supplement.
#' 
#+
salba <- ps %>% 
  ps_filter(!str_detect(SAMPLE, "plus"))
nrow(otu_table(salba))
6697 #12133

salba_plus <- ps %>% 
  ps_filter(str_detect(SAMPLE, "plus"))
nrow(otu_table(salba_plus))
8693 #19731

nrow(otu_table(ps))
ntaxa(ps)
10137 #26611

#'
#' Total number of unique sequences
#'
#+
sum(otu_table(ps))
1439201 #1455675

#'
#' ## Rarefaction
#'
#+
psr <- rarefy_even_depth(ps, rngseed = 11)
salbar <- psr %>% 
  ps_filter(!str_detect(SAMPLE, "plus"))
ntaxa(salbar)
3231 #4209

salbar_plus <- psr %>% 
  ps_filter(str_detect(SAMPLE, "plus"))
ntaxa(salbar_plus)
4380 #5917

sum(otu_table(psr))
245448
ntaxa(psr)
8255 # number of unique otus

otus <- otu_table(psr)
colnames(otus) <- sample_data(psr)$sample_name
otus


#'
#' - Venn diagrams of the shared OTUs across the zones 
#' (including individual samples 1-6). 
#'It is  with numbers indicating the total number of OTU for each group
#'
#'
#' 1. All samples merged by Zone and Rhizo
#' 
#'
#+
groups <- as.factor(paste0(sample_data(ps)$"zone", ifelse(as.logical(sample_data(pss)$"rhizo"), "+", "")))
pssm <- ps %>% 
  rarefy_even_depth(rngseed = 11) %>% 
  merge_samples(group = groups, fun = sum)

#' Parsing sample data
#+
samples <- sample_data(pssm) %>% 
  as.data.frame()
samples$SAMPLE <- row.names(samples)
samples$zone <- str_extract(samples$SAMPLE, "^Z\\d")
samples$rhizo <- ifelse(samples$rhizo > 0, 1, 0)
samples$sample_name <- samples$SAMPLE
sample_data(pssm) <- samples

otu_table(pssm)[1:6, 1:10]

#' 
#' Subset a groups for OTU sets comparison.
#+ 
d <- t(otu_table(pssm)) > 0
d <- as_tibble(d@.Data, rownames = "otu")

#'
#' Venn diagrams.
#'
#'
#' Kui palju OTUsid on per tsoon
#' 
#' 
#' Rarefied data
#+
colSums(t(otu_table(pssm)) > 0)
# Z1  Z1+   Z2  Z2+   Z3  Z3+ 
# 2172 3690 1894 1654 1311 1893 

#' 
#' Non-rarefied data
#+
colSums(t(otu_table(merge_samples(ps, group = groups, fun = sum))) > 0)
# Z1  Z1+   Z2  Z2+   Z3  Z3+ 
# 6356 9423 4989 4488 3537 9414 

#' 
#' Common OTUs in unsupplemented soils
#'
#+
fill_col <- ggsci::pal_npg()(6)
text_size <-  3.5
set_name_size <- 5
show_perc <- FALSE
no_rhiz <- ggvenn(d, c("Z1", "Z2", "Z3"), fill_color = fill_col, text_size = text_size, set_name_size = set_name_size, show_percentage = show_perc, stroke_color = NA)

get_common_otus <- function(x) {
  x %>% 
    rowwise() %>% 
    apply(1, all) %>% 
    sum()
}

d[c("Z1", "Z2", "Z3")] %>% get_common_otus()
321 #323

#' 
#' Common OTUs in +rhizo supplemented soils
#'
#+
rhiz <- d %>% 
  rename(`Z1+rhiz` = `Z1+`, `Z2+rhiz` = `Z2+`, `Z3+rhiz` = `Z3+`) %>% 
  ggvenn(c("Z1+rhiz", "Z2+rhiz", "Z3+rhiz"), fill_color = fill_col, text_size = text_size, set_name_size = set_name_size, show_percentage = show_perc, stroke_color = NA)
d[c("Z1+", "Z2+", "Z3+")] %>% get_common_otus()
346 #341

a <- no_rhiz + rhiz
ggsave("plots/venn_panela.png", a)


#'
#' Proportion of OTUs common across three sites w/wo treatment
#'
#+
data <- tibble(
  ser = c("No", "Yes"),
  total = c(ntaxa(ps_filter(pssm, !str_detect(SAMPLE, "\\+"))), ntaxa(ps_filter(pssm, str_detect(SAMPLE, "\\+")))),
  common = c(d[c("Z1", "Z2", "Z3")] %>% get_common_otus(), d[c("Z1+", "Z2+", "Z3+")] %>% get_common_otus())
  )

get_prior(
  common | trials(total) ~ ser,
  data = data,
  family = binomial()
)

m1 <- brm(
  common | trials(total) ~ ser,
  data = data,
  family = binomial(),
  prior = prior(normal(0, 0.5), class = b),
  file = here("models/common | trials(total) ~ ser"),
  file_refit = "on_change"
)

pp_check(m1)
plot(conditional_effects(m1))

hypothesis(m1, "serYes = 0")
as.data.frame(m1) %>% 
  mean_hdci(
    no = inv_logit_scaled(b_Intercept), 
    yes = inv_logit_scaled(b_Intercept + b_serYes),
    es = no - yes, 
    h1 = as.numeric(es > 0),
    .width = 0.95) %>% 
  mutate_at(vars(matches("no|yes|es")), scales::percent, accuracy = 0.01)

#'
#' Is there any difference in number of common OTUs within each site w/wo 
#' treatment, what's the proportion of common OTUs
#+
d_renamed <- d %>% 
  rename(`Z1+rhiz` = `Z1+`, `Z2+rhiz` = `Z2+`, `Z3+rhiz` = `Z3+`)
z1 <- ggvenn(d_renamed, c("Z1", "Z1+rhiz"), fill_color = fill_col, text_size = text_size, set_name_size = set_name_size, show_percentage = show_perc, stroke_color = NA)
z2 <- ggvenn(d_renamed, c("Z2", "Z2+rhiz"), fill_color = fill_col, text_size = text_size, set_name_size = set_name_size, show_percentage = show_perc, stroke_color = NA)
z3 <- ggvenn(d_renamed, c("Z3", "Z3+rhiz"), fill_color = fill_col, text_size = text_size, set_name_size = set_name_size, show_percentage = show_perc, stroke_color = NA)

b <- z1 + z2 + z3
ggsave("plots/venn_panelb.png", b)

#'
#+
data <- tibble(
  zone = rep(str_c("Z", 1:3), each = 2),
  ser = rep(c("No", "Yes"), 3),
  common = rep(c(
    d[c("Z1", "Z1+")] %>% get_common_otus(), 
    d[c("Z2", "Z2+")] %>% get_common_otus(),
    d[c("Z3", "Z3+")] %>% get_common_otus()
  ), each = 2),
  total = c(
    sum(d["Z1"]),
    sum(d["Z1+"]),
    sum(d["Z2"]),
    sum(d["Z2+"]),
    sum(d["Z3"]),
    sum(d["Z3+"])
    )
)

get_prior(
  common | trials(total) ~ ser + zone + ser:zone,
  data = data,
  family = binomial()
)

m2 <- brm(
  common | trials(total) ~ ser + zone + ser:zone,
  data = data,
  family = binomial(),
  prior = prior(normal(0, 0.5), class = b),
  iter = 2600,
  file = here("models/common | trials(total) ~ ser + zone + ser:zone")
)

pp_check(m2)
mp <- plot(conditional_effects(m2), plot = FALSE)
pd <- position_dodge(0.3)
venn_df_model <- data %>% 
  add_linpred_draws(m2) %>% 
  ggplot(aes(x = zone, y = inv_logit_scaled(.linpred), color = ser)) +
  stat_pointinterval(position = pd) +
  scale_color_discrete("S. europaea\nrhizosphere") +
  labs(y = "Proportion of\ncommon OTUs") +
  theme(
    axis.title.x = element_blank(),
    legend.position = "bottom", 
    legend.title = element_text(vjust = 2)
    )

#'
#' OTUs common to all sites
#'
#+
d[-1] %>% get_common_otus()
198 #203

nrow(d)
5776 #8255

comm1 <- brm(
  common | trials(total) ~ 1,
  data = tibble(common = d[-1] %>% get_common_otus(), total = nrow(d)),
  family = binomial()
)

fixef(comm1) %>% inv_logit_scaled()


c <- d %>% 
  rowwise() %>% 
  mutate(
    `-rhiz` = any(Z1, Z2, Z3),
    `+rhiz` = any(`Z1+`, `Z2+`, `Z3+`)
  ) %>% 
  ggvenn(c("-rhiz", "+rhiz"), fill_color = fill_col, text_size = text_size, set_name_size = set_name_size, show_percentage = show_perc, stroke_color = NA)
ggsave(here("plots/venn_panelc.png"), c)

venns <- (no_rhiz | rhiz | c) / (z1 | z2 | z3) / ((plot_spacer() | venn_df_model | plot_spacer()) + plot_layout(widths = c(1, 6, 1))) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "plain"))
ggsave(here("plots/venn_diagrams.png"), width = 20, height = 18, units = "cm")
tiff(here("plots/venn_diagrams.tiff"), width = 20, height = 18, units = "cm", res = 300)
venns
dev.off()

#'
#' ## Diversity indices
#'
#+
rich <- estimate_richness(ps, measures = c("Observed", "Shannon", "InvSimpson")) %>% 
  as_tibble(rownames = "library")
alpha <- as_tibble(sample_data(ps)) %>% 
  mutate_at("zone", ~as.numeric(str_extract(.x, "\\d"))) %>% 
  bind_cols(rich) %>% 
  mutate(zone = str_c("Z", zone), rhizo = ifelse(as.logical(rhizo), "+rhiz", "-rhiz"))

libsize <- otu_table(ps) %>% 
  colSums()
alpha$log_libsize <- log(libsize)
mean(log(libsize))

#'
#' Observed nr of OTUs.
#'
#+
get_prior(
  Observed ~ log_libsize + rhizo + zone + rhizo:zone,
  data = alpha,
  family = negbinomial()
)

m3o <- brm(
  Observed ~ log_libsize + rhizo + zone + rhizo:zone,
  data = alpha,
  family = negbinomial(),
  prior = c(prior("normal(0, 1)", class = "b")),
  file = here("models/Observed ~ log_libsize + rhizo + zone + rhizo:zone"),
  file_refit = "on_change",
  iter = 2000
)

summary(m3o)
pp_check(m3o)

po <- plot(conditional_effects(m3o), plot = FALSE)
po$zone
po$rhizo

newdata <- alpha %>% 
  data_grid(rhizo, zone, log_libsize = 10) %>% 
  add_epred_draws(m3o)

oh <- newdata %>% 
  group_by(rhizo) %>% 
  mutate(id = row_number()) %>% 
  select(id, zone, rhizo, .draw, .epred) %>% 
  pivot_wider(names_from = rhizo, values_from = .epred) %>% 
  mutate(es = `+rhiz` - `-rhiz`, h1 = as.numeric(es > 0)) %>% 
  group_by(zone) %>% 
  mean_qi(`+rhiz`, `-rhiz`, es, h1)

stars <- oh %>% 
  filter(h1 >= 0.95) %>% 
  select(zone) %>% 
  inner_join(alpha) %>% 
  group_by(zone) %>% 
  filter(Observed == max(Observed)) %>% 
  select(zone, Observed) %>% 
  mutate(label = "*")

observed <- newdata %>% 
  ggplot() +
  stat_pointinterval(aes(x = factor(zone), y = .epred, color = rhizo), position = pd) +
  geom_jitter(data = alpha, aes(x = factor(zone), y = Observed, color = rhizo), position = position_jitterdodge(jitter.width = 0.1, dodge.width = pd$width)) +
  geom_text(data = stars, aes(zone, Observed, label = label), inherit.aes = FALSE, nudge_y = 100, size = 5) +
  scale_color_discrete("S. europaea\nrhizosphere", labels = c("No", "Yes")) +
  labs(y = "Observed") +
  theme(
    axis.title.x = element_blank()
    )
ggsave(here("plots/diversity_observed.png"))
hypothesis(m3o, "zoneZ3 < 0")
hypothesis(m3o, "zoneZ2 < 0")
hypothesis(m3o, "exp(zoneZ3) < exp(zoneZ2)")
newdata %>% 
  filter(rhizo == "-rhiz") %>% 
  group_by(zone) %>% 
  mutate(id = row_number()) %>% 
  select(id, zone, .draw, .epred) %>% 
  pivot_wider(names_from = zone, values_from = .epred) %>% 
  mutate(
    es12 = Z1 - Z2,
    es13 = Z1 - Z3,
    es23 = Z2 - Z3,
    h12 = as.numeric(es12 > 0),
    h13 = as.numeric(es13 > 0),
    h23 = as.numeric(es23 > 0)
  ) %>% 
  pivot_longer(cols = c("Z1", "Z2", "Z3", "es12", "es13", "es23", "h12", "h13", "h23")) %>% 
  group_by(name) %>% 
  mean_hdci(value) %>% 
  view()


# Observed number of taxa
newdata %>% 
  group_by(zone, rhizo) %>% 
  mean_qi(.epred)

#'
#' Shannon diversity index.
#'
#+
get_prior(
  Shannon ~ log_libsize + zone + rhizo + zone:rhizo,
  data = alpha,
  family = gaussian()
)

m3s <- brm(
  Shannon ~ log_libsize + zone + rhizo + zone:rhizo,
  data = alpha,
  family = student(),
  prior = prior("normal(0, 1)", class = "b"),
  file = here("models/Shannon ~ log_libsize + zone + rhizo + zone:rhizo"),
  file_refit = "on_change",
  sample_prior = "yes"
)

pp_check(m3s)
psh <- plot(conditional_effects(m3s), plot = FALSE)
psh$zone
psh$rhizo

summary(m3s)

newdata <- alpha %>% 
  data_grid(rhizo, zone, log_libsize = 10) %>% 
  add_epred_draws(m3s)

sh <- newdata %>% 
  group_by(rhizo) %>% 
  mutate(id = row_number()) %>% 
  select(id, zone, rhizo, .draw, .epred) %>% 
  pivot_wider(names_from = rhizo, values_from = .epred) %>% 
  mutate(es = `+ser` - `-ser`, h1 = as.numeric(es > 0)) %>% 
  group_by(zone) %>% 
  mean_qi(`+ser`, `-ser`, es, h1)

hypothesis(m3s, "zoneZ3 < 0", alpha = 0.05)
hypothesis(m3s, "zoneZ2 < 0", alpha = 0.05)
hypothesis(m3s, "zoneZ3 < zoneZ2", alpha = 0.05)

stars <- sh %>% 
  filter(h1 >= 0.95) %>% 
  select(zone) %>% 
  inner_join(alpha) %>% 
  group_by(zone) %>% 
  filter(Shannon == max(Shannon)) %>% 
  select(zone, Shannon) %>% 
  mutate(label = "*")

shannon <- newdata %>% 
  ggplot() +
  stat_pointinterval(aes(x = factor(zone), y = .epred, color = rhizo), position = pd) +
  geom_jitter(data = alpha, aes(x = factor(zone), y = Shannon, color = rhizo), position = position_jitterdodge(jitter.width = 0.1, dodge.width = pd$width)) +
  geom_text(data = stars, aes(zone, Shannon, label = label), inherit.aes = FALSE, nudge_y = 0.1, size = 5) +
  scale_color_discrete("S. europaea\nrhizosphere", labels = c("No", "Yes")) +
  labs(y = "Shannon index") +
  theme(
    axis.title.x = element_blank()
  )

ggsave("plots/diversity_shannon.png", shannon)

#'
#' Inverse Simpson diversity index.
#'
#+
get_prior(
  InvSimpson ~ log_libsize + zone + rhizo + zone:rhizo,
  data = alpha,
  family = gaussian()
)

m3is <- brm(
  InvSimpson ~ log_libsize + zone + rhizo + zone:rhizo,
  data = alpha,
  family = student(),
  prior = prior("student_t(3, 0, 1)", class = "b"),
  file = here("models/InvSimpson ~ log_libsize + zone + rhizo + zone:rhizo"),
  file_refit = "on_change"
)

pp_check(m3is)
pis <- plot(conditional_effects(m3is), plot = FALSE)
pis$zone

newdata <- alpha %>% 
  data_grid(rhizo, zone, log_libsize = 10) %>% 
  add_epred_draws(m3is)

newdata %>% 
  group_by(rhizo) %>% 
  mutate(id = row_number()) %>% 
  select(id, zone, rhizo, .draw, .epred) %>% 
  pivot_wider(names_from = rhizo, values_from = .epred) %>% 
  mutate(es = `+rhiz` - `-rhiz`, h1 = as.numeric(es > 0)) %>% 
  group_by(zone) %>% 
  mean_qi(`+rhiz`, `-rhiz`, es, h1)

invsimpson <- newdata %>% 
  ggplot() +
  stat_pointinterval(aes(x = factor(zone), y = .epred, color = rhizo), position = pd) +
  geom_jitter(data = alpha, aes(x = factor(zone), y = InvSimpson, color = rhizo), position = position_jitterdodge(jitter.width = 0.1, dodge.width = pd$width)) +
  scale_color_discrete("S. europaea\nrhizosphere", labels = c("No", "Yes")) +
  labs(y = "InvSimpson") +
  theme(
    axis.title.x = element_blank()
  )

ggsave("plots/diversity_invsimpson.png", shannon)
hypothesis(m3is, "zoneZ3 < 0", alpha = 0.05)
hypothesis(m3is, "zoneZ2 < 0", alpha = 0.05)
hypothesis(m3is, "zoneZ3 < zoneZ2", alpha = 0.05)

#'
#' Alpha diversity plot
#'
#+
diversity <- ((observed | shannon | invsimpson) + 
  plot_layout(guides = "collect") & theme(legend.position = "bottom")) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "plain"))

ggsave(here("plots/ritchness_plot.png"), width = 25, height = 16, units = "cm")
ggsave(here("plots/beta_diversity.png"), plot = betap, width = 20, height = 10, units = "cm", dpi = 300)
tiff(here(glue::glue("plots/alpha_diversity.tiff")), width = 20, height = 10, units = "cm", res = 300)
diversity
dev.off()

#' 
#' ## Ordination plot
#' 
#' - Nonmetric multidimensional plot displaying beta diversity by Bray-Curis 
#' dissimilarities (based on Hellinger-transformed abundances) for microbial 
#' communities
#' 
#' NMD plot shows gradient across zones. No clear separation of rhizo treated or untreated samples within zones.
#'
#+ fig.cap="NMD plot." 
set.seed(11)
ord <- ordinate(ps, "NMDS", "bray")
p2 <- plot_ordination(
  ps %>% ps_mutate(rhizo = ifelse(rhizo == 0, "no supp", "+Se rhiz")), 
  ord, 
  type = "samples", 
  color = "zone",
  shape = "rhizo",
  justDF = TRUE
) 
pnmds <- p2 %>% 
  ggplot(aes(NMDS1, NMDS2, color = zone)) +
  # stat_ellipse(aes(linetype = rhizo)) +
  stat_ellipse(linetype = 2) +
  geom_point(aes(shape = rhizo), size = 3) +
  coord_fixed() +
  theme(legend.title = element_blank())
ggsave(here("plots/nmds.png"), pnmds, width = 12, height = 8, units = "cm", dpi = 300)

#' 
#' ## 
#'
#' - Boxplots for each of the sites, the distances to the centroid
#'
#+
dist <- p2 %>% 
  group_by(zone, rhizo) %>% 
  nest() %>% 
  mutate(
    cd = map(data, cendist),
    cd = map(cd, matrix, ncol =  1)
  ) %>% 
  select(zone, cd) %>% 
  unnest(cd) %>% 
  mutate_at("cd", as.numeric)

dist_centroid <- p2 %>% 
  mutate(dist = dist$cd) 

m5 <- brm(
  bf(dist ~ zone + rhizo, sigma ~ zone + rhizo),
  data = dist_centroid,
  family = student(),
  prior = prior("normal(0, 1)", class = "b"),
  sample_prior = "yes",
  file = here("models/bf(dist ~ zone + rhizo, sigma ~ zone + rhizo)"),
  file_refit = "on_change"
)
pp_check(m5)
plot(conditional_effects(m5), points = TRUE, ask = FALSE)
hyp <- c("exp(sigma_Intercept) = 0",
         "exp(sigma_Intercept + sigma_zoneZ2) = 0",
         "exp(sigma_Intercept + sigma_zoneZ3) = 0"
)
hypothesis(m5, hyp)
hyp <- "exp(sigma_Intercept + sigma_zoneZ2) > exp(sigma_Intercept)"
(hyp <- hypothesis(m5, hyp))
hyp <- "exp(sigma_Intercept + sigma_zoneZ2 + sigma_rhizonosupp) < exp(sigma_Intercept + sigma_zoneZ2)"
(hyp <- hypothesis(m5, hyp))

pdc <- dist_centroid %>% 
  ggplot(aes(x = zone, y = dist, color = rhizo))  +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) +
  labs(x = "Zone", y  = "Dist. to centroid") +
  theme(legend.title = element_blank(), axis.title.x = element_blank()) 
ggsave(here("plots/distance_to_centroid_zones_rhizo.png"), pdc, width = 12, height = 8, units = "cm", dpi = 300)

betap <- pnmds + pdc + plot_annotation(tag_levels = "A")
ggsave(here("plots/beta_diversity.png"), plot = betap, width = 20, height = 10, units = "cm", dpi = 300)
tiff(here(glue::glue("plots/beta_diversity.tiff")), width = 20, height = 10, units = "cm", res = 300)
betap
dev.off()

#'
#'Halomonadaceae – unclassified
#'Idiomarina
#'Halomonas
#'Marinobacteria
#'
#+
n <- 20
rank <- "Genus"
top_taxa <- tax_top(ps, n = n, rank = rank)
highlight <- str_to_lower(top_taxa) %>% 
  str_subset("mari|halo") %>% 
  str_to_sentence()
htmp <- ps %>% 
  tax_transform("clr", rank = rank, zero_replace = "halfmin") %>%
  comp_heatmap(
    taxa = top_taxa, 
    taxon_renamer = function(x) str_replace(x, "_unclassified", " u."),
    tax_seriation = "Heatmap",
    row_names_gp = grid::gpar(col = ifelse(top_taxa %in% highlight, "green", "black")),
    sample_seriation = "Identity",
    sample_anno = sampleAnnotation(
      Zone = anno_sample("zone"),
      Rhizo = anno_sample("rhizo", fun = as.character),
      col = list(
        Zone = c("Z1" = "green", "Z2" = "yellow", "Z3" = "red"),
        Rhizo = c("0" = "gray85", "1" = "gray25"))
    ),
    show_row_dend = FALSE,
    # colors = heat_palette(palette = "BluYl", rev = TRUE),
    row_names_side = "left", sample_side = "bottom",
    use_raster = TRUE
  )
tiff(here(glue::glue("plots/heatmap_top{n}.tiff")), width = 20, height = 20, units = "cm", res = 300)
ComplexHeatmap::draw(
  object = htmp,
  annotation_legend_list = attr(htmp, "AnnoLegends"), merge_legends = TRUE
)
dev.off()

rank <- "Phylum"
top_taxa <- tax_top(ps, n = n, rank = rank)
subset_taxa(pss, Phylum %in% top_taxa) %>% 
  tax_table()
p_phy <- ps %>% 
  comp_barplot(tax_level = "Phylum", sample_order = "default", n_taxa = 10, label = "sample_name") +
  coord_flip()
ggsave(here("plots/comp_barplot_phy.png"), width = 20, height = 20, units = "cm", dpi = 300)


proteobacteria <- subset_taxa(ps, Phylum == "Proteobacteria")
p_pro <- proteobacteria %>% 
  comp_barplot(tax_level = "Genus", sample_order = "default", n_taxa = 10, label = "sample_name") +
  coord_flip()
ggsave(here("plots/comp_barplot_proteobacteria.png"), width = 20, height = 20, units = "cm", dpi = 300)

p_bar <- (p_phy + theme(plot.tag.position  = c(0, 1.015))) + 
  (p_pro + theme(plot.tag.position  = c(0, 1.015))) + 
  plot_annotation(tag_levels = "A")
ggsave(here("plots/comp_barplot.png"), width = 40, height = 20, units = "cm", dpi = 300)

