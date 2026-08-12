# =============================================================================
#  Figure S1 - RINENG-D-26-06081 (Revision R1)
#  Carte des individus de l'ACP : position des sept temps d'echantillonnage
#  dans le plan Dim1-Dim2
#
#  Source : pca_individual_coordinates_uniform_k8.csv
#           (scenario retenu : imputation kNN uniforme, k_NN = 8)
#
#  Les pourcentages de variance sont recalcules a partir des coordonnees, sans
#  saisie manuelle : pour une ACP normee, la valeur propre d'un axe vaut la
#  somme des carres des coordonnees des individus divisee par (n - 1).
#  Le controle croise avec le fichier des variables donne les memes valeurs.
# =============================================================================

source(file.path("scripts", "_helpers.R"))
assert_repository_root()
load_required_packages(c("ggplot2", "ggrepel"))


# -----------------------------------------------------------------------------
# Chemins et parametres
# -----------------------------------------------------------------------------
CSV    <- file.path("results", "tables", "pca_individual_coordinates_uniform_k8.csv")
OUTDIR <- file.path("results", "figures")

W <- 8.0 ; H <- 6.5 ; UNITS <- "in" ; DPI <- 1200

# Times New Roman. Sous Windows la police est disponible d'office ; ailleurs,
# l'alias "serif" du peripherique graphique y renvoie.
if (.Platform$OS.type == "windows") {
  windowsFonts(Times = windowsFont("Times New Roman"))
  FAM <- "Times"
} else {
  FAM <- "serif"
}

if (!file.exists(CSV))
  stop("Fichier introuvable : ", CSV,
       "\n  R exige / ou \\\\ dans les chemins, pas un antislash simple.")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Lecture et variance expliquee
# -----------------------------------------------------------------------------
d <- read.csv(CSV, stringsAsFactors = FALSE)

dims <- grep("^Dim\\.", names(d), value = TRUE)
eig  <- colSums(d[dims]^2) / (nrow(d) - 1)
pct  <- 100 * eig / sum(eig)

cat(sprintf("Dim1 = %.2f %%   Dim2 = %.2f %%   (cumul %.2f %%)\n",
            pct[1], pct[2], pct[1] + pct[2]))

# Les temps sont ordonnes numeriquement, et non alphabetiquement : sans cela
# t12 precederait t0 dans la legende comme dans l'ordre du trace.
d$t_num <- as.numeric(sub("^t", "", d$Time))
d       <- d[order(d$t_num), ]

# -----------------------------------------------------------------------------
# 2. Trace
# -----------------------------------------------------------------------------
p <- ggplot(d, aes(x = Dim.1, y = Dim.2)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey35", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey35", linewidth = 0.5) +
  # Trajectoire chronologique. Les temps ayant ete tries numeriquement plus
  # haut, geom_path suit l'ordre t0 -> t72 et non l'ordre du fichier. Le trait
  # est gris et la fleche unique, portee par le dernier segment : la
  # trajectoire guide la lecture sans passer devant les points.
  geom_path(colour = "grey50", linewidth = 0.6,
            arrow = arrow(length = unit(2.8, "mm"), type = "closed",
                          ends = "last")) +
  geom_point(size = 2.6, colour = "black") +
  geom_text_repel(aes(label = Time),
                  size = 5.2, family = FAM, colour = "black",
                  min.segment.length = Inf, box.padding = 0.45,
                  point.padding = 0.25, force = 1, seed = 1) +
  labs(x = sprintf("Dim1 (%.1f%%)", pct[1]),
       y = sprintf("Dim2 (%.1f%%)", pct[2])) +
  theme_bw(base_size = 14, base_family = FAM) +
  theme(
    axis.title       = element_text(face = "bold", size = 17),
    axis.text        = element_text(size = 14, colour = "black"),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "grey30", linewidth = 0.6),
    plot.margin      = margin(6, 8, 6, 6)
  )

print(p)

# -----------------------------------------------------------------------------
# 3. Export : TIFF 1200 dpi + PDF vectoriel
# -----------------------------------------------------------------------------
tif <- file.path(OUTDIR, "Figure_S1.tiff")
pdf <- file.path(OUTDIR, "Figure_S1.pdf")

ggsave(tif, p, width = W, height = H, units = UNITS, dpi = DPI, compression = "lzw")
ggsave(pdf, p, width = W, height = H, units = UNITS, device = cairo_pdf)

mo <- file.info(tif)$size / 1024^2
message("Ecriture dans : ", OUTDIR)
message(sprintf("  Figure_S1.tiff : %6.1f Mo%s", mo,
                if (mo > 10) "  <-- depasse la limite Elsevier de 10 Mo" else ""))
message(sprintf("  Figure_S1.pdf  : %6.1f Mo", file.info(pdf)$size / 1024^2))

# -----------------------------------------------------------------------------
# Reglages fins :
#   size = 5.8 dans geom_text_repel   -> etiquettes de temps plus grandes
#   size = 3.2 dans geom_point        -> points plus gros
#   box.padding = 0.6                 -> etiquettes plus ecartees des points
#   W <- 190 ; H <- 155 ; UNITS <- "mm"   -> taille finale d'impression
#
# Pour supprimer la trajectoire, retirer le bloc geom_path().
# Pour une fleche sur chaque segment : ends = "both" -> supprimer l'argument
# ends, ou passer a linetype = "dashed" pour un trait discontinu.
# =============================================================================
