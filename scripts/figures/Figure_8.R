# =============================================================================
#  Figure 8 - RINENG-D-26-06081 (Revision R1)
#  Cercle des correlations de l'ACP sur les coefficients de partage D_i(t)
#
#  Source : pca_variable_coordinates_uniform_k8.csv (scenario retenu :
#           imputation kNN uniforme, k_NN = 8)
#
#  Les valeurs propres, les pourcentages de variance et les contributions sont
#  recalcules a partir des coordonnees du fichier, sans saisie manuelle :
#    - valeur propre d'un axe = somme des carres des coordonnees sur cet axe
#    - % de variance          = valeur propre / somme des valeurs propres
#    - contribution au plan   = moyenne des contributions aux deux axes,
#                               ponderee par leurs valeurs propres
#      (definition utilisee par factoextra::fviz_pca_var)
# =============================================================================

source(file.path("scripts", "_helpers.R"))
assert_repository_root()
load_required_packages(c("ggplot2", "dplyr", "ggrepel"))


# -----------------------------------------------------------------------------
# Chemins
# -----------------------------------------------------------------------------
OUTDIR <- file.path("results", "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
CSV     <- file.path("results", "tables", "pca_variable_coordinates_uniform_k8.csv")

if (!file.exists(CSV))
  stop("Fichier introuvable : ", CSV,
       "\n  Verifier le chemin relatif depuis la racine du depot.")

# -----------------------------------------------------------------------------
# 1. Lecture et calcul des grandeurs derivees
# -----------------------------------------------------------------------------
d <- read.csv(CSV, stringsAsFactors = FALSE)

dims <- grep("^Dim\\.", names(d), value = TRUE)
eig  <- colSums(d[dims]^2)          # valeurs propres
pct  <- 100 * eig / sum(eig)        # % de variance expliquee

# Le rang de l'ACP vaut 6 pour 7 individus : la somme des valeurs propres vaut
# 17 x 6/7 = 14.571 et non 17. Les pourcentages sont donc rapportes a cette
# somme, ce qui redonne exactement les valeurs de la Table S4.
cat(sprintf("Dim1 = %.2f %%   Dim2 = %.2f %%   (cumul %.2f %%)\n",
            pct[1], pct[2], pct[1] + pct[2]))

# Contribution de chaque variable au plan Dim1-Dim2
c1 <- 100 * d$Dim.1^2 / eig[1]
c2 <- 100 * d$Dim.2^2 / eig[2]
d$contrib <- (c1 * eig[1] + c2 * eig[2]) / (eig[1] + eig[2])

# Etiquettes : seul le symbole de l'element est affiche. Le prefixe du
# coefficient est retire, qu'il s'agisse de "D" (coefficient de partage) ou de
# "Z" (coefficient centre-reduit) selon le fichier source. La classe [DZ] ne
# retire que le PREMIER caractere : "ZZn" donne bien "Zn", "ZSO3" donne "SO3".
d$Variable <- sub("^[DZ]", "", d$Variable)

if (any(grepl("^[DZ][A-Z]", d$Variable)))
  warning("Certaines etiquettes semblent conserver un prefixe : ",
          paste(d$Variable[grepl("^[DZ][A-Z]", d$Variable)], collapse = ", "))

# Indices typographiques : "Al2O3" devient "Al[2]*O[3]", interprete par
# plotmath via parse = TRUE. Les chiffres sont automatiquement places en
# indice, quel que soit l'element.
to_plotmath <- function(x) sub("\\*$", "", gsub("([0-9]+)", "[\\1]*", x))
d$lab <- to_plotmath(d$Variable)

# Position des etiquettes : chacune est placee dans le prolongement de sa
# propre fleche, juste au-dela de la pointe.
#
# Une fleche plus courte que ses voisines immediates verrait toutefois son
# etiquette retomber SUR la fleche voisine (cas de Cu, K2O, MgO, Zn, Na2O,
# Al2O3). Le rayon retenu est donc le plus grand des rayons observes dans un
# secteur angulaire de +/- ANG_TOL degres autour de la fleche, augmente du
# decalage. Seules les etiquettes concernees sont deplacees ; les autres
# restent a OFFSET de leur pointe.
r        <- sqrt(d$Dim.1^2 + d$Dim.2^2)
ang      <- atan2(d$Dim.2, d$Dim.1) * 180 / pi

OFFSET   <- 0.09     # ecart minimal entre une pointe et son etiquette
ANG_TOL  <- 12       # demi-secteur angulaire, en degres

r_sect <- sapply(seq_along(ang), function(i) {
  diff <- abs(((ang - ang[i] + 180) %% 360) - 180)   # ecart angulaire signe
  max(r[diff < ANG_TOL])                             # inclut la fleche elle-meme
})

# Ajustement manuel, element par element, pour les cas ou le secteur angulaire
# ne suffit pas. Ajouter une entree pour ecarter davantage une etiquette,
# une valeur negative pour la rapprocher.
EXTRA <- c(CaO = 0.05, Cr = 0.05, MgO = 0.05)

extra <- EXTRA[d$Variable]
extra[is.na(extra)] <- 0
extra <- unname(extra)

lab_r   <- r_sect + OFFSET + extra
d$lab_x <- d$Dim.1 / r * lab_r
d$lab_y <- d$Dim.2 / r * lab_r

cat("Etiquettes decalees au-dela d'une fleche voisine :",
    paste(d$Variable[r_sect > r + 1e-9], collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# 2. Cercle unite
# -----------------------------------------------------------------------------
theta  <- seq(0, 2 * pi, length.out = 400)
circle <- data.frame(x = cos(theta), y = sin(theta))

# -----------------------------------------------------------------------------
# 3. Trace
# -----------------------------------------------------------------------------
p <- ggplot() +
  geom_path(data = circle, aes(x, y), colour = "grey55", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  geom_segment(data = d,
               aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2, colour = contrib),
               arrow = arrow(length = unit(2.2, "mm"), type = "closed"),
               linewidth = 0.6) +
  # Etiquettes au bout de leur fleche. parse = TRUE active les indices.
  # force est faible pour que ggrepel n'intervienne que sur les collisions
  # reelles et laisse chaque etiquette pres de sa pointe.
  # min.segment.length = Inf supprime les traits de rappel, qui donnaient
  # l'impression qu'une etiquette designait la fleche voisine.
  geom_text_repel(data = d,
                  aes(x = lab_x, y = lab_y, label = lab, colour = contrib),
                  parse = TRUE,
                  size = 4.6, fontface = "plain", show.legend = FALSE,
                  min.segment.length = Inf, box.padding = 0.12,
                  point.padding = 0, force = 0.35, max.overlaps = Inf,
                  seed = 1) +
  # Degrade de factoextra : bleu-vert -> jaune -> orange-rouge. Ce sont les
  # trois teintes de reference (gradient.cols de fviz_pca_var) ; les nuances
  # intermediaires vertes et orangees resultent de leur interpolation.
  scale_colour_gradientn(colours = c("#00AFBB", "#E7B800", "#FC4E07"),
                         name = "contrib") +
  coord_fixed(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)) +
  scale_x_continuous(breaks = seq(-1, 1, 0.5)) +
  scale_y_continuous(breaks = seq(-1, 1, 0.5)) +
  labs(x = sprintf("Dim1 (%.1f%%)", pct[1]),
       y = sprintf("Dim2 (%.1f%%)", pct[2])) +
  theme_bw(base_size = 14) +
  theme(
    axis.title       = element_text(face = "bold", size = 16),
    axis.text        = element_text(size = 13, colour = "black"),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "grey35", linewidth = 0.5),
    legend.title     = element_text(face = "bold", size = 13),
    legend.text      = element_text(size = 12),
    legend.key.height = unit(12, "mm"),
    plot.margin      = margin(4, 4, 4, 4)
  )

print(p)

# -----------------------------------------------------------------------------
# 4. Export : TIFF 1200 dpi + PDF vectoriel, dans le meme dossier
# -----------------------------------------------------------------------------
W <- 8.0 ; H <- 7.0 ; UNITS <- "in" ; DPI <- 1200

tif <- file.path(OUTDIR, "Figure_8.tiff")
pdf <- file.path(OUTDIR, "Figure_8.pdf")

ggsave(tif, p, width = W, height = H, units = UNITS, dpi = DPI, compression = "lzw")
ggsave(pdf, p, width = W, height = H, units = UNITS, device = cairo_pdf)

mo <- file.info(tif)$size / 1024^2
message("Ecriture dans : ", OUTDIR)
message(sprintf("  Figure_8.tiff : %6.1f Mo%s", mo,
                if (mo > 10) "  <-- depasse la limite Elsevier de 10 Mo" else ""))
message(sprintf("  Figure_8.pdf  : %6.1f Mo", file.info(pdf)$size / 1024^2))

# -----------------------------------------------------------------------------
# Tableau de controle : coordonnees et contributions, a reporter dans le texte
# -----------------------------------------------------------------------------
print(d %>% select(Variable, Dim.1, Dim.2, Dim.3, contrib) %>%
        arrange(desc(contrib)), row.names = FALSE, digits = 3)

# -----------------------------------------------------------------------------
# Reglages fins :
#   size = 5.2 dans geom_text_repel        -> etiquettes plus grandes
#   fontface = "bold" dans geom_text_repel -> etiquettes en gras
#   OFFSET  <- 0.06                        -> etiquettes encore plus pres
#   EXTRA <- c(CaO = 0.08, ...)            -> ecarter une etiquette precise
#   ANG_TOL <- 8                           -> secteur plus etroit : moins
#                                             d'etiquettes decalees, mais
#                                             risque accru d'interception
#   force = 0.6                            -> ggrepel ecarte davantage
#   arrow(length = unit(3, "mm"))          -> pointes de fleches plus grandes
#   W <- 190 ; H <- 148 ; UNITS <- "mm"    -> taille finale d'impression
# =============================================================================
