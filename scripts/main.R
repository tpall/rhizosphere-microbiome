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
ps <- phyloseq_validate(ps, remove_undetected = TRUE)
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
    rhizo = as.numeric(str_detect(SAMPLE, "plus"))
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
# pss <- ps
# otu_table(pss)[otu_table(ps) / colSums(otu_table(ps)) < 0.001] <- 0
# all.equal(otu_table(ps), otu_table(pss))
# pss <- ps_filter(pss)

pss <- subset_taxa(ps, rowSums(otu_table(ps)) > 1) # rowSums(otu_table(ps)) != 1 # apply(otu_table(ps) / colSums(otu_table(ps)), 1, max) < 0.001

#'
#' Top10 tax overview.
#+
pss %>% 
  comp_barplot(tax_level = "Genus", sample_order = "default", n_taxa = 10, label = "sample_name") +
  coord_flip()

#' 
#' Number of OTUs in *S. alba* cuttings with and without *S. europea* rhizosphere 
#' supplement.
#' 
#+
salba <- pss %>% 
  ps_filter(!str_detect(SAMPLE, "plus"))
nrow(otu_table(salba))
6697 #12133

salba_plus <- pss %>% 
  ps_filter(str_detect(SAMPLE, "plus"))
nrow(otu_table(salba_plus))
8693 #19731

nrow(otu_table(pss))
10137 #26611

#'
#' Total number of unique sequences
#'
#+
sum(otu_table(pss))
1439201

#'
#'
#'
#+
psr <- rarefy_even_depth(pss, rngseed = 11)
salbar <- psr %>% 
  ps_filter(!str_detect(SAMPLE, "plus"))
nrow(otu_table(salbar))
3231 #4209

salbar_plus <- psr %>% 
  ps_filter(str_detect(SAMPLE, "plus"))
ntaxa(salbar_plus)
4380 #5917

sum(otu_table(psr))
245448
ntaxa(psr)
8255 # number of unique otus

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
groups <- as.factor(paste0(sample_data(pss)$"zone", ifelse(sample_data(pss)$"rhizo" == 1, "+Se-rhiz", "")))
pssm <- merge_samples(psr, group = groups, fun = sum)
samples <- sample_data(pssm) %>% 
  as.data.frame()
samples$SAMPLE <- row.names(samples)
samples$zone <- str_extract(samples$SAMPLE, "^Z\\d")
samples$rhizo <- if_else(samples$rhizo > 0, 1, 0)
samples$sample_name <- samples$SAMPLE
sample_data(pssm) <- samples

otu_table(pssm)[1:6, 1:10]

#' 
#' Subset a groups for OTU sets comparison.
#+ 
samples$SAMPLE
z1 <- pssm %>% 
  subset_samples(SAMPLE == "Z1") %>% 
  filter_taxa(function(x) x > 0, TRUE) %>% 
  taxa_names()
z1rhiz <- pssm %>% 
  subset_samples(SAMPLE == "Z1+Se-rhiz") %>% 
  filter_taxa(function(x) x > 0, TRUE) %>% 
  taxa_names()
z2 <- pssm %>% 
  subset_samples(SAMPLE == "Z2") %>% 
  filter_taxa(function(x) x > 0, TRUE) %>% 
  taxa_names()
z2rhiz <- pssm %>% 
  subset_samples(SAMPLE == "Z2+Se-rhiz") %>% 
  filter_taxa(function(x) x > 0, TRUE) %>% 
  taxa_names()
z3 <- pssm %>% 
  subset_samples(SAMPLE == "Z3") %>% 
  filter_taxa(function(x) x > 0, TRUE) %>% 
  taxa_names()
z3rhiz <- pssm %>% 
  subset_samples(SAMPLE == "Z3+Se-rhiz") %>% 
  filter_taxa(function(x) x > 0, TRUE) %>% 
  taxa_names()
otu.sets <- list("Z1" = z1, "Z1+Se-rhiz" = z1rhiz, "Z2" = z2, "Z2+Se-rhiz" = z2rhiz, "Z3" = z3, "Z3+Se-rhiz" = z3rhiz)


#'
#' Venn diagrams.
#'
#'
#' Kui palju OTUsid on per tsoon
#+
map_dbl(otu.sets, length) %>% sort()

#' 
#' Common OTUs in unsupplemented soils
#'
#+
no_rhiz <- ggplot() + 
  geom_venn(otu.sets[c(1, 3, 5)]) +
  theme_void() +
  theme(legend.position = "none")  +
  coord_fixed(ratio = 1, clip = "off")
length(Reduce(intersect, otu.sets[c(1, 3, 5)]))
321 #323

#' 
#' Common OTUs in Se rhizo supplemented soils
#'
#+
rhiz <- ggplot() + 
  geom_venn(otu.sets[c(2, 4, 6)]) +
  theme_void() +
  theme(legend.position = "none") +
  coord_fixed(ratio = 1, clip = "off")
length(Reduce(intersect, otu.sets[c(2, 4, 6)]))
346 #341

a <- no_rhiz + rhiz
ggsave("plots/venn_panela.png", a)


#'
#' Proportion of OTUs common across three sites w/wo treatment
#'
#+
data <- tibble(
  ser = c("No", "Yes"),
  total = c(ntaxa(salbar), ntaxa(salbar_plus)),
  common = c(length(Reduce(intersect, otu.sets[c(1, 3, 5)])), length(Reduce(intersect, otu.sets[c(2, 4, 6)])))
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
conditional_effects(m1)

hypothesis(m1, "serYes = 0")
as.data.frame(m1) %>% 
  mean_qi(
    no = exp(b_Intercept), 
    yes = exp(b_Intercept + b_serYes),
    es = no - yes) %>% 
  mutate_at(vars(matches("no|yes|es")), scales::percent, accuracy = 0.1)

#'
#' Is there any difference in number of common OTUs within each site w/wo 
#' treatment, what's the proportion of common OTUs
#+
z1 <- ggplot() + 
  geom_venn(otu.sets[c(1, 2)]) +
  theme_void() +
  theme(legend.position = "none") +
  coord_fixed(ratio = 1/1.5, clip = "off")

z2 <- ggplot() + 
  geom_venn(otu.sets[c(3, 4)]) +
  theme_void() +
  theme(legend.position = "none") +
  coord_fixed(ratio = 1/1.5, clip = "off")

z3 <- ggplot() + 
  geom_venn(otu.sets[c(5, 6)]) +
  theme_void() +
  theme(legend.position = "none") +
  coord_fixed(ratio = 1/1.5, clip = "off")

b <- z1 + z2 + z3
ggsave("plots/venn_panelb.png", b)

#'
#+
data <- tibble(
  zone = rep(str_c("Z", 1:3), each = 2),
  ser = rep(c("No", "Yes"), 3),
  common = rep(c(
    length(Reduce(intersect, otu.sets[c(1, 2)])), 
    length(Reduce(intersect, otu.sets[c(3, 4)])),
    length(Reduce(intersect, otu.sets[c(5, 6)]))
  ), each = 2),
  total = c(
    length(otu.sets[[1]]),
    length(otu.sets[[2]]),
    length(otu.sets[[3]]),
    length(otu.sets[[4]]),
    length(otu.sets[[5]]),
    length(otu.sets[[6]])
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
data %>% 
  add_epred_draws(m2) %>% 
  mutate(.epred = .epred / total) %>% 
  ggplot(aes(x = zone, y = .epred, color = ser)) +
  stat_pointinterval() +
  stat_lineribbon(.width = 0.95, alpha = 0.3)


#'
#' Decrease in OTUs  across zones
#' 
get_prior(
  total ~ ser + zone + ser:zone,
  data = data,
  family = negbinomial()
)

m2.1 <- brm(
  total ~ ser + zone + ser:zone,
  data = data,
  family = poisson(),
  prior = prior(normal(0, 0.5), class = b),
  iter = 2600,
  file = here("models/total ~ ser + zone + ser:zone"),
  file_refit = "on_change"
)

plot(conditional_effects(m2.1), ask = FALSE)

m2.2 <- brm(
  total ~ 1 + (ser | zone),
  data = data,
  family = poisson(),
  # prior = prior(normal(0, 0.5), class = b),
  iter = 2000,
  file = here("models/total ~ 1 + (ser | zone)"),
  file_refit = "on_change"
)

pp_check(m2.2)
plot(conditional_effects(m2.2), ask = FALSE)

draws <- data %>% 
  add_epred_draws(m2.2)
draws %>% 
  ggplot(aes(zone, .epred, color = ser)) +
  stat_pointinterval()
draws %>% 
  mean_qi()


#'
#' OTUs common to all sites
#'
#+
all <- Reduce(intersect, otu.sets[c(1, 3, 5, 2, 4, 6)])
length(all)
198 #203

otu.sets %>% 
  unlist() %>% 
  unique() %>% 
  length()
5776 #8255


no_rhiz_all <- otu.sets[c(1, 3, 5)] %>% unlist() %>% unique()
rhiz_all <- otu.sets[c(2, 4, 6)] %>% unlist() %>% unique()
c <- ggplot() + 
  geom_venn(list("no supp" = no_rhiz_all, "+Se-rhiz" = rhiz_all)) +
  theme_void() +
  theme(legend.position = "none") +
  coord_fixed(ratio = 0.67,  clip = "off")
ggsave(here("plots/venn_panelc.png"), c)

design <- "
AABBCC
AABBCC
EEEEEE
EEEEEE
"
venns <- a + c + b + 
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A")
ggsave(here("plots/venn_diagrams.png"), width = 20, height = 12, units = "cm")

#'
#' ## Diversity indices
#'
#+
richness <- estimate_richness(pss, measures = c("Observed", "Shannon", "InvSimpson"))
samples <- sample_data(pss) %>% 
  as.matrix() %>% 
  as_tibble() %>% 
  mutate_at("SAMPLE", str_replace_all, "-", ".") %>% 
  mutate_at("zone", ~as.numeric(str_extract(.x, "\\d")))
alpha <- richness %>% 
  rownames_to_column(var="SAMPLE") %>% 
  as_tibble() %>% 
  mutate(Evenness = Shannon / log(Observed)) %>% 
  left_join(samples)

#'
#' Observed nr of OTUs.
#'
#+
get_prior(
  Observed ~ zone + rhizo + zone:rhizo,
  data = alpha,
  family = negbinomial()
)

m3o <- brm(
  Observed ~ zone + rhizo + zone:rhizo,
  data = alpha,
  family = negbinomial(),
  prior = prior("student_t(3, 0, 1)", class = "b"),
  file = here("models/Observed ~ zone + rhizo + zone:rhizo"),
  file_refit = "on_change",
  sample_prior = "yes"
)

m3o1 <- brm(
  Observed ~ 1 + (rhizo | zone),
  data = alpha,
  family = negbinomial(),
  # prior = prior("student_t(3, 0, 1)", class = "b"),
  file = here("models/Observed ~ 1 + (rhizo | zone)"),
  file_refit = "on_change",
  sample_prior = "yes"
)

pp_check(m3o1)
loo(m3o, m3o1)

pp_check(m3o)
po <- plot(conditional_effects(m3o), plot = FALSE)
summary(m3o)
as.data.frame(m3o)
hypothesis(m3o, "zone = 0", alpha = 0.05)
hypothesis(m3o, "zone + rhizo1 + zone:rhizo1 = 0", alpha = 0.05)
pp1o <- po$zone + scale_x_continuous("Zone", breaks = c(1, 2, 3))
pp2o <- po$`zone:rhizo` +   
  geom_text(data = tibble(zone = c(2, 3), Observed = c(1700, 1600), label = c("*", "*")), aes(zone, Observed, label = label), inherit.aes = FALSE) +
  scale_x_continuous("Zone", breaks = c(1, 2, 3)) +
  theme(axis.title.y = element_blank(),
        legend.position = "none")

pd <- position_dodge(0.8)
pp3o <- alpha %>% 
  ggplot(aes(as.factor(zone), Observed, color = rhizo)) +
  geom_boxplot(position = pd, outlier.color = NA, witdh = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) + 
  scale_x_discrete("Zone", breaks = c(1, 2, 3)) +
  theme(legend.position = "none")
observed <- pp3o / (pp1o + pp2o)
ggsave("plots/diversity_observed.png", observed)

newdata <- alpha %>% 
  add_epred_draws(m3o)
newdata %>% 
  group_by(rhizo) %>% 
  mutate(id = row_number()) %>% 
  select(id, zone, rhizo, .draw, .epred) %>% 
  pivot_wider(names_from = rhizo, values_from = .epred) %>% 
  group_by(zone) %>% 
  mutate(es = `1` - `0`, h1 = es > 0) %>% 
  mean_qi(`0`, `1`, es, h1)

# Observed number of taxa
newdata %>% 
  group_by(zone, rhizo) %>% 
  mean_qi(.epred)

#'
#' Shannon diversity index.
#'
#+
get_prior(
  Shannon ~ zone + rhizo + zone:rhizo,
  data = alpha,
  family = gaussian()
)

m3s <- brm(
  Shannon ~ zone + rhizo + zone:rhizo,
  data = alpha,
  family = student(),
  prior = prior("normal(0, 1)", class = "b"),
  file = here("models/Shannon ~ zone + rhizo + zone:rhizo"),
  file_refit = "on_change",
  sample_prior = "yes"
)

pp_check(m3s)
ps <- plot(conditional_effects(m3s), plot = FALSE)
summary(m3s)
hypothesis(m3s, "zone = 0", alpha = 0.05)
hypothesis(m3o, "zone + rhizo1 + zone:rhizo1 = 0", alpha = 0.05)
pp1s <- ps$zone + scale_x_continuous("Zone", breaks = c(1, 2, 3))
newdata <- alpha %>% 
  select(zone, rhizo) %>% 
  distinct() %>% 
  add_epred_draws(m3s)
newdata %>% 
  group_by(zone) %>% 
  select(zone, .draw, rhizo, .epred) %>% 
  pivot_wider(names_from = rhizo, values_from = .epred) %>% 
  mutate(es = `1` - `0`, h1 = es > 0) %>% 
  mean_qi(`0`, `1`, es, h1)

pp2s <- ps$`zone:rhizo` +   
  # geom_text(data = tibble(zone = c(2, 3), Shannon = c(4.5, 4.3), label = c("*", "*")), aes(zone, Shannon, label = label), inherit.aes = FALSE) +
  scale_x_continuous("Zone", breaks = c(1, 2, 3)) +
  theme(axis.title.y = element_blank(),
        legend.position = "none")

pp3s <- alpha %>% 
  ggplot(aes(as.factor(zone), Shannon, color = rhizo)) +
  geom_boxplot(position = pd, outlier.color = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) + 
  scale_x_discrete("Zone", breaks = c(1, 2, 3)) +
  theme(legend.position = "none")
shannon <- pp3s / (pp1s + pp2s) 
ggsave("plots/diversity_shannon.png", shannon)

#'
#' Inverse Simpson diversity index.
#'
#+
get_prior(
  InvSimpson ~ zone + rhizo + zone:rhizo,
  data = alpha,
  family = gaussian()
)

m3is <- brm(
  InvSimpson ~ zone + rhizo + zone:rhizo,
  data = alpha,
  family = student(),
  prior = prior("student_t(3, 0, 1)", class = "b"),
  file = here("models/InvSimpson ~ zone + rhizo + zone:rhizo"),
  file_refit = "on_change"
)

pp_check(m3is)
pis <- plot(conditional_effects(m3is), plot = FALSE)
summary(m3is)
hypothesis(m3is, "zone = 0", alpha = 0.05)
hypothesis(m3is, "rhizo1 = 0", alpha = 0.05)
hypothesis(m3is, "zone + rhizo1 < 0", alpha = 0.05)

pp1is <- pis$zone + scale_x_continuous("Zone", breaks = c(1, 2, 3))
newdata <- alpha %>% 
  select(zone, rhizo) %>% 
  distinct() %>% 
  add_epred_draws(m3is)
newdata %>% 
  group_by(zone) %>% 
  select(zone, .draw, rhizo, .epred) %>% 
  pivot_wider(names_from = rhizo, values_from = .epred) %>% 
  mutate(es = `1` - `0`, h1 = es > 0) %>% 
  mean_qi(`0`, `1`, es, h1)

pp2is <- pis$`zone:rhizo` +   
  scale_x_continuous("Zone", breaks = c(1, 2, 3)) +
  theme(axis.title.y = element_blank()) +
  scale_color_discrete("", labels = c("no supp", "+Se rhiz")) +
  scale_fill_discrete("", labels = c("no supp", "+Se rhiz"))

pp3is <- alpha %>% 
  ggplot(aes(as.factor(zone), InvSimpson, color = rhizo)) +
  geom_boxplot(position = pd, outlier.color = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) + 
  scale_x_discrete("Zone", breaks = c(1, 2, 3)) +
  scale_color_discrete("", labels = c("no supp", "+Se rhiz"))
invsimpson <- pp3is / (pp1is + pp2is) 
ggsave("plots/diversity_invsimpson.png", invsimpson)

#'
#' Evenness
#'
#+
get_prior(
  Evenness ~ zone + rhizo + zone:rhizo,
  data = alpha,
  family = gaussian()
)

m3e <- brm(
  Evenness ~ zone + rhizo + zone:rhizo,
  data = alpha,
  family = Beta(),
  prior = prior("student_t(3, 0, 1)", class = "b"),
  file = here("models/Evenness ~ zone + rhizo + zone:rhizo"),
  file_refit = "on_change"
)

pp_check(m3e)
pe <- plot(conditional_effects(m3e), plot = FALSE)
summary(m3e)
hypothesis(m3e, "zone = 0", alpha = 0.05)
hypothesis(m3e, "rhizo1 = 0", alpha = 0.05)
hypothesis(m3e, "zone + rhizo1 + zone:rhizo1 = 0", alpha = 0.05)

pp1e <- pe$zone + scale_x_continuous("Zone", breaks = c(1, 2, 3))
newdata <- alpha %>% 
  select(zone, rhizo) %>% 
  distinct() %>% 
  add_epred_draws(m3e)
newdata %>% 
  group_by(zone) %>% 
  select(zone, .draw, rhizo, .epred) %>% 
  pivot_wider(names_from = rhizo, values_from = .epred) %>% 
  mutate(es = `1` - `0`, h1 = es > 0) %>% 
  mean_qi(`0`, `1`, es, h1)

pp2e <- pe$`zone:rhizo` +   
  scale_x_continuous("Zone", breaks = c(1, 2, 3)) +
  theme(axis.title.y = element_blank()) +
  scale_color_discrete("", labels = c("no supp", "+Se rhiz")) +
  scale_fill_discrete("", labels = c("no supp", "+Se rhiz"))

pp3e <- alpha %>% 
  ggplot(aes(as.factor(zone), Evenness, color = rhizo)) +
  geom_boxplot(position = pd, outlier.color = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) + 
  scale_x_discrete("Zone", breaks = c(1, 2, 3)) +
  scale_color_discrete("", labels = c("no supp", "+Se rhiz"))
evenness <- pp3e / (pp1e + pp2e) 
ggsave("plots/diversity_evenness.png", evenness)

#'
#' Alpha diversity plot
#'
#+
diversity <- observed | shannon | invsimpson + 
  plot_layout(guides = "collect")

diversity + plot_annotation(tag_levels = "A")
ggsave(here("plots/ritchness_plot.png"), width = 25, height = 16, units = "cm")

evenness2 <- (pp3e + theme(legend.position = "none")) / (pp1e + (pp2e + theme(legend.position = "none"))) 
diversity4 <- observed | shannon | evenness2 | invsimpson + 
  plot_layout(guides = "collect")

diversity4 + plot_annotation(tag_levels = "A")
ggsave(here("plots/ritchness_plot4.png"), width = 25, height = 16, units = "cm")

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
ord <- ordinate(pss, "NMDS", "bray")
p2 <- plot_ordination(
  pss %>% ps_mutate(rhizo = ifelse(rhizo == 0, "no supp", "+Se rhiz")), 
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

pnmds + pdc + plot_annotation(tag_levels = "A")
ggsave(here("plots/beta_diversity.png"), width = 20, height = 10, units = "cm", dpi = 300)

#'
#'Halomonadaceae – unclassified
#'Idiomarina
#'Halomonas
#'Marinobacteria
#'
#+
n <- 20
rank <- "Genus"
top_taxa <- tax_top(psr, n = n, rank = rank)
highlight <- str_to_lower(top_taxa) %>% 
  str_subset("mari|halo") %>% 
  str_to_sentence()
htmp <- pss %>% 
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
png(here(glue::glue("plots/heatmap_top{n}.png")), width = 20, height = 20, units = "cm", res = 300)
ComplexHeatmap::draw(
  object = htmp,
  annotation_legend_list = attr(htmp, "AnnoLegends"), merge_legends = TRUE
)
dev.off()

rank <- "Phylum"
top_taxa <- tax_top(pss, n = n, rank = rank)
subset_taxa(pss, Phylum %in% top_taxa) %>% 
  tax_table()
p_phy <- pss %>% 
  comp_barplot(tax_level = "Phylum", sample_order = "default", n_taxa = 10, label = "sample_name") +
  coord_flip()
ggsave(here("plots/comp_barplot_phy.png"), width = 20, height = 20, units = "cm", dpi = 300)


proteobacteria <- subset_taxa(pss, Phylum == "Proteobacteria")
p_pro <- proteobacteria %>% 
  comp_barplot(tax_level = "Genus", sample_order = "default", n_taxa = 10, label = "sample_name") +
  coord_flip()
ggsave(here("plots/comp_barplot_proteobacteria.png"), width = 20, height = 20, units = "cm", dpi = 300)

p_bar <- (p_phy + theme(plot.tag.position  = c(0, 1.015))) + 
  (p_pro + theme(plot.tag.position  = c(0, 1.015))) + 
  plot_annotation(tag_levels = "A")
ggsave(here("plots/comp_barplot.png"), width = 40, height = 20, units = "cm", dpi = 300)

