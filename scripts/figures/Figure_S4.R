# =============================================================================
#  Figure S4 - RINENG-D-26-06081 (Revision R1)
#  Robustesse au retrait d'un temps (leave-one-time-out, LOTO)
#
#  (a) ACP : variance cumulee PC1+PC2 et similarite de Procrustes de la
#            configuration des variables, par rapport a l'analyse complete
#  (b) HCA : indice de Rand ajuste et accord d'appartenance, par rapport a la
#            partition de reference G = 3
#  (c) Stabilite par espece : part des sept analyses LOTO ou chaque espece
#      conserve son groupe de reference
#
#  Les panneaux (a) et (b) sont exprimes en pourcentage, ce qui permet de
#  superposer deux indicateurs par panneau sur une echelle unique, sans second
#  axe.
#
#  Sources : leave_one_time_out_pca_sensitivity.csv
#            leave_one_time_out_hca_sensitivity_G3.csv
#            leave_one_time_out_variable_hca_stability_G3.csv
# =============================================================================

source(file.path("scripts", "_helpers.R"))
assert_repository_root()
load_required_packages(c("ggplot2", "scales", "patchwork"))


# -----------------------------------------------------------------------------
# Chemins et parametres
# -----------------------------------------------------------------------------
TABDIR <- file.path("results", "tables")
OUTDIR <- file.path("results", "figures")

CSV_PCA <- file.path(TABDIR, "leave_one_time_out_pca_sensitivity.csv")
CSV_HCA <- file.path(TABDIR, "leave_one_time_out_hca_sensitivity_G3.csv")
CSV_STA <- file.path(TABDIR, "leave_one_time_out_variable_hca_stability_G3.csv")

W <- 11.0 ; H <- 9.4 ; UNITS <- "in" ; DPI <- 1200

# Times New Roman. Sous Windows la police est disponible d'office ; ailleurs,
# l'alias "serif" du peripherique graphique y renvoie.
if (.Platform$OS.type == "windows") {
  windowsFonts(Times = windowsFont("Times New Roman"))
  FAM <- "Times"
} else {
  FAM <- "serif"
}

for (f in c(CSV_PCA, CSV_HCA, CSV_STA))
  if (!file.exists(f))
    stop("Fichier introuvable : ", f,
         "\n  R exige / ou \\\\ dans les chemins, pas un antislash simple.")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Lecture
#    Les temps sont ordonnes numeriquement : sans cela t12 precederait t0.
# -----------------------------------------------------------------------------
pca <- read.csv(CSV_PCA, stringsAsFactors = FALSE)
hca <- read.csv(CSV_HCA, stringsAsFactors = FALSE)
sta <- read.csv(CSV_STA, stringsAsFactors = FALSE)

ordonner <- function(d) {
  d$t_num <- as.numeric(sub("^t", "", d$Omitted_time))
  d <- d[order(d$t_num), ]
  d$Omitted_time <- factor(d$Omitted_time, levels = d$Omitted_time)
  d
}
pca <- ordonner(pca); hca <- ordonner(hca)

cat("ACP : similarite de Procrustes (variables, 2D) minimale pour",
    as.character(pca$Omitted_time[which.min(pca$Variable_Procrustes_similarity_2D)]),
    sprintf("(%.3f)\n", min(pca$Variable_Procrustes_similarity_2D)))
cat("HCA : indice de Rand ajuste minimal pour",
    as.character(hca$Omitted_time[which.min(hca$Adjusted_Rand_Index)]),
    sprintf("(%.3f)\n", min(hca$Adjusted_Rand_Index)))
cat("Especes parfaitement stables :",
    sum(sta$Stability_percent == 100), "sur", nrow(sta), "\n")

# -----------------------------------------------------------------------------
# 2. Mise en forme longue : deux indicateurs par panneau
#    La similarite de Procrustes et l'indice de Rand, definis entre 0 et 1,
#    sont exprimes en pourcentage pour partager l'echelle des deux autres.
# -----------------------------------------------------------------------------
d_pca <- rbind(
  data.frame(Omitted_time = pca$Omitted_time,
             Indicateur = "Cumulative PC1 + PC2",
             Valeur     = pca$PC1_PC2_cumulative_percent),
  data.frame(Omitted_time = pca$Omitted_time,
             Indicateur = "Procrustes similarity (variables, 2D) \u00d7 100",
             Valeur     = 100 * pca$Variable_Procrustes_similarity_2D))
# La similarite de Procrustes et l'indice de Rand ajuste sont definis entre
# 0 et 1 : le facteur 100 est indique dans le libelle, faute de quoi une
# echelle en pourcentage laisserait croire qu'ils sont deja des pourcentages.

d_hca <- rbind(
  data.frame(Omitted_time = hca$Omitted_time,
             Indicateur = "Membership agreement",
             Valeur     = hca$Membership_agreement_percent),
  data.frame(Omitted_time = hca$Omitted_time,
             Indicateur = "Adjusted Rand index \u00d7 100",
             Valeur     = 100 * hca$Adjusted_Rand_Index))

# -----------------------------------------------------------------------------
#    Panneau (c) : le prefixe "Z" est retire et les chiffres passes en indice
#    ("ZAl2O3" -> Al[2]*O[3]), l'expression etant interpretee par plotmath.
#    Les especes sont classees par stabilite croissante, les plus sensibles se
#    trouvant ainsi en haut du diagramme.
# -----------------------------------------------------------------------------
to_plotmath <- function(x) sub("\\*$", "", gsub("([0-9]+)", "[\\1]*", x))

sta$lab <- to_plotmath(sub("^Z", "", sta$Variable))
sta     <- sta[order(sta$Stability_percent, decreasing = TRUE), ]
sta$lab <- factor(sta$lab, levels = sta$lab)

# Les trois classes de sensibilite du fichier source sont conservees telles
# quelles ; l'ordre fixe evite un classement alphabetique de la legende.
sta$Sensitivity_class <- factor(
  sta$Sensitivity_class,
  levels = c("Very stable", "Moderately sensitive", "Sensitive"))

# -----------------------------------------------------------------------------
# 3. Theme commun
# -----------------------------------------------------------------------------
theme_S4 <- theme_bw(base_size = 14, base_family = FAM) +
  theme(
    axis.title       = element_text(face = "bold", size = 15),
    axis.text        = element_text(size = 13, colour = "black"),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "grey30", linewidth = 0.6),
    legend.title     = element_blank(),
    legend.text      = element_text(size = 11),
    legend.key       = element_blank(),
    legend.key.width = unit(7, "mm"),
    legend.key.height = unit(5, "mm"),
    legend.background = element_rect(fill = "white", colour = "grey50",
                                     linewidth = 0.4),
    plot.margin      = margin(6, 8, 6, 6)
  )

COL2 <- c("#1F4E79", "#B03A2E")

# Les libelles portant desormais la mention "x 100", ils sont trop larges
# pour tenir dans les zones libres des panneaux : la legende est placee
# au-dessus de chaque panneau, comme dans le panneau (c). Aucune donnee n'est
# alors masquee, et la longueur des libelles ne contraint plus la mise en page.
legende_dessus <- function(p) {
  p + theme(legend.position = "top",
            legend.background = element_blank(),
            legend.margin = margin(0, 0, 2, 0),
            legend.box.margin = margin(0, 0, -4, 0))
}

# Les bornes suivent au plus pres les donnees : chaque legende trouve place
# dans une zone deja vide du panneau, sans qu'il soit necessaire d'elargir
# l'echelle. Les graduations s'arretent a 100 %, valeur maximale que les
# indicateurs peuvent atteindre.
panneau <- function(d, ylab, ymin, ymax) {
  ggplot(d, aes(x = Omitted_time, y = Valeur,
                colour = Indicateur, shape = Indicateur, group = Indicateur)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.6) +
    scale_colour_manual(values = setNames(COL2, unique(d$Indicateur))) +
    scale_shape_manual(values  = setNames(c(16, 17), unique(d$Indicateur))) +
    scale_y_continuous(limits = c(ymin, ymax),
                       breaks = scales::breaks_pretty()(c(ymin, 100))) +
    scale_x_discrete(labels = function(l) paste0(sub("^t", "", l), " h")) +
    labs(x = "Omitted sampling time", y = ylab) +
    theme_S4
}

# -----------------------------------------------------------------------------
# 4. Panneaux
#    Les bornes basses different : la HCA descend a 31 % (indice de Rand pour
#    t12), l'ACP reste au-dessus de 80 %. Une echelle commune ecraserait le
#    panneau (a) ; chacun garde donc la sienne.
#    Les bornes suivent au plus pres les donnees : 82.0 % est la plus petite
#    valeur du panneau (a), 31.0 % celle du panneau (b).
# -----------------------------------------------------------------------------
p_pca <- legende_dessus(panneau(d_pca, "Percentage (%)", 80, 101))
p_hca <- legende_dessus(panneau(d_hca, "Percentage (%)", 28, 101))

# -----------------------------------------------------------------------------
#    Panneau (c) : stabilite par espece
# -----------------------------------------------------------------------------
p_sta <- ggplot(sta, aes(x = Stability_percent, y = lab,
                         fill = Sensitivity_class)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%.1f", Stability_percent)),
            hjust = -0.18, size = 4.2, family = FAM) +
  scale_y_discrete(labels = function(l) parse(text = l), limits = rev) +
  scale_x_continuous(limits = c(0, 118), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(name = NULL,
                    values = c("Very stable"          = "#1F4E79",
                               "Moderately sensitive" = "#E8871A",
                               "Sensitive"            = "#B03A2E")) +
  labs(x = "Group stability across the seven LOTO analyses (%)", y = NULL) +
  theme_S4 +
  theme(legend.position = "top",
        legend.background = element_blank(),
        legend.margin = margin(0, 0, 0, 0))

p_all <- (p_pca | p_hca) / p_sta +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 18, face = "bold", family = FAM))

print(p_all)

# -----------------------------------------------------------------------------
# 5. Export : TIFF 1200 dpi + PDF vectoriel
# -----------------------------------------------------------------------------
tif <- file.path(OUTDIR, "Figure_S4.tiff")
pdf <- file.path(OUTDIR, "Figure_S4.pdf")

ggsave(tif, p_all, width = W, height = H, units = UNITS, dpi = DPI, compression = "lzw")
ggsave(pdf, p_all, width = W, height = H, units = UNITS, device = cairo_pdf)

mo <- file.info(tif)$size / 1024^2
message("Ecriture dans : ", OUTDIR)
message(sprintf("  Figure_S4.tiff : %6.1f Mo%s", mo,
                if (mo > 10) "  <-- depasse la limite Elsevier de 10 Mo" else ""))
message(sprintf("  Figure_S4.pdf  : %6.1f Mo", file.info(pdf)$size / 1024^2))

# -----------------------------------------------------------------------------
# Reglages fins :
#   Pour replacer la legende dans le panneau, il faut raccourcir les libelles
#   ("Procrustes similarity x 100" suffit, le detail passant dans la legende
#   de la figure) puis remplacer legende_dessus() par :
#     p + theme(legend.position = "inside",
#               legend.position.inside = c(0.70, 0.72))
#   (p_pca / p_hca) au lieu de (p_pca | p_hca) -> panneaux empiles
#   W <- 190 ; H <- 83 ; UNITS <- "mm"      -> taille finale d'impression
#
# =============================================================================
