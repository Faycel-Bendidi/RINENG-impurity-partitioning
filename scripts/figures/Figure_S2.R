# =============================================================================
#  Figure S2 - RINENG-D-26-06081 (Revision R1)
#  Diagnostics du choix du nombre de groupes de la HCA
#
#  (a) Coude : somme des carres intra-groupes, G = 1 a 10
#  (b) Largeur de silhouette moyenne, G = 2 a 10
#
#  Sources : hca_elbow_diagnostics_uniform_k8.csv
#            hca_silhouette_diagnostics_uniform_k8.csv
#            (scenario retenu : imputation kNN uniforme, k_NN = 8)
# =============================================================================

source(file.path("scripts", "_helpers.R"))
assert_repository_root()
load_required_packages(c("ggplot2", "patchwork"))


# -----------------------------------------------------------------------------
# Chemins et parametres
# -----------------------------------------------------------------------------
TABDIR <- file.path("results", "tables")
OUTDIR <- file.path("results", "figures")

CSV_ELB <- file.path(TABDIR, "hca_elbow_diagnostics_uniform_k8.csv")
CSV_SIL <- file.path(TABDIR, "hca_silhouette_diagnostics_uniform_k8.csv")

W <- 11.0 ; H <- 4.6 ; UNITS <- "in" ; DPI <- 1200

# Times New Roman. Sous Windows la police est disponible d'office ; ailleurs,
# l'alias "serif" du peripherique graphique y renvoie.
if (.Platform$OS.type == "windows") {
  windowsFonts(Times = windowsFont("Times New Roman"))
  FAM <- "Times"
} else {
  FAM <- "serif"
}

for (f in c(CSV_ELB, CSV_SIL))
  if (!file.exists(f))
    stop("Fichier introuvable : ", f,
         "\n  R exige / ou \\\\ dans les chemins, pas un antislash simple.")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Lecture
# -----------------------------------------------------------------------------
elb <- read.csv(CSV_ELB, stringsAsFactors = FALSE)
sil <- read.csv(CSV_SIL, stringsAsFactors = FALSE)

cat("Coude      : G de", min(elb$Number_of_groups_G), "a", max(elb$Number_of_groups_G), "\n")
cat("Silhouette : G de", min(sil$Number_of_groups_G), "a", max(sil$Number_of_groups_G), "\n")
cat("Silhouette moyenne maximale a G =",
    sil$Number_of_groups_G[which.max(sil$Mean_silhouette_width)], "\n")

# -----------------------------------------------------------------------------
# 2. Theme commun aux deux panneaux
# -----------------------------------------------------------------------------
theme_S2 <- theme_bw(base_size = 14, base_family = FAM) +
  theme(
    axis.title       = element_text(face = "bold", size = 15),
    axis.text        = element_text(size = 13, colour = "black"),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "grey30", linewidth = 0.6),
    plot.margin      = margin(6, 8, 6, 6)
  )

# -----------------------------------------------------------------------------
# 3. Panneaux
# -----------------------------------------------------------------------------
p_elb <- ggplot(elb, aes(x = Number_of_groups_G, y = Total_within_cluster_SS)) +
  geom_line(linewidth = 0.7, colour = "black") +
  geom_point(size = 2.2, colour = "black") +
  scale_x_continuous(breaks = 1:10) +
  labs(x = "Number of groups (G)", y = "Total within-cluster sum of squares") +
  theme_S2

p_sil <- ggplot(sil, aes(x = Number_of_groups_G, y = Mean_silhouette_width)) +
  geom_line(linewidth = 0.7, colour = "black") +
  geom_point(size = 2.2, colour = "black") +
  scale_x_continuous(breaks = 2:10, limits = c(2, 10)) +
  labs(x = "Number of groups (G)", y = "Mean silhouette width") +
  theme_S2

# -----------------------------------------------------------------------------
# 4. Assemblage, etiquettes (a) et (b)
# -----------------------------------------------------------------------------
p_all <- (p_elb | p_sil) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 18, face = "bold", family = FAM))

print(p_all)

# -----------------------------------------------------------------------------
# 5. Export : TIFF 1200 dpi + PDF vectoriel
# -----------------------------------------------------------------------------
tif <- file.path(OUTDIR, "Figure_S2.tiff")
pdf <- file.path(OUTDIR, "Figure_S2.pdf")

ggsave(tif, p_all, width = W, height = H, units = UNITS, dpi = DPI, compression = "lzw")
ggsave(pdf, p_all, width = W, height = H, units = UNITS, device = cairo_pdf)

mo <- file.info(tif)$size / 1024^2
message("Ecriture dans : ", OUTDIR)
message(sprintf("  Figure_S2.tiff : %6.1f Mo%s", mo,
                if (mo > 10) "  <-- depasse la limite Elsevier de 10 Mo" else ""))
message(sprintf("  Figure_S2.pdf  : %6.1f Mo", file.info(pdf)$size / 1024^2))

# -----------------------------------------------------------------------------
# MARQUAGE DE G = 3 -- a n'activer qu'apres comparaison avec le dendrogramme
#
# Ajouter dans p_elb et p_sil, AVANT geom_line() pour que le trait passe
# derriere la courbe :
#
#   geom_vline(xintercept = 3, linetype = "dashed",
#              colour = "grey45", linewidth = 0.6) +
#
# Ou, pour marquer seulement le point retenu, apres geom_point() :
#
#   geom_point(data = subset(elb, Number_of_groups_G == 3),
#              size = 4, shape = 21, stroke = 1.1,
#              colour = "red", fill = NA) +
#
# -----------------------------------------------------------------------------
# Reglages fins :
#   size = 2.8 dans geom_point        -> points plus gros
#   W <- 190 ; H <- 80 ; UNITS <- "mm"    -> taille finale d'impression
#   (p_elb / p_sil) au lieu de (p_elb | p_sil) -> panneaux empiles
# =============================================================================
