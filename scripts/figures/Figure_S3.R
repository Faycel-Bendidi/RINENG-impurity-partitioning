# =============================================================================
#  Figure S3 - RINENG-D-26-06081 (Revision R1)
#  Choix et validation de l'imputation kNN
#
#  (a) Validation croisee : RMSE en fonction de k, pour les deux ponderations
#  (b) Sensibilite locale : amplitude relative de chacune des six valeurs
#      imputees lorsque k varie
#
#  Sources : knn_cross_validation_rmse_mae_delta.csv
#            local_k_sensitivity_summary.csv
# =============================================================================

source(file.path("scripts", "_helpers.R"))
assert_repository_root()
load_required_packages(c("ggplot2", "patchwork"))


# -----------------------------------------------------------------------------
# Chemins et parametres
# -----------------------------------------------------------------------------
RESDIR <- file.path("data", "processed")
OUTDIR <- file.path("results", "figures")

CSV_CV  <- file.path(RESDIR, "knn_cross_validation_rmse_mae_delta.csv")
CSV_SEN <- file.path(RESDIR, "local_k_sensitivity_summary.csv")

W <- 11.0 ; H <- 4.6 ; UNITS <- "in" ; DPI <- 1200

# Times New Roman. Sous Windows la police est disponible d'office ; ailleurs,
# l'alias "serif" du peripherique graphique y renvoie.
if (.Platform$OS.type == "windows") {
  windowsFonts(Times = windowsFont("Times New Roman"))
  FAM <- "Times"
} else {
  FAM <- "serif"
}

for (f in c(CSV_CV, CSV_SEN))
  if (!file.exists(f))
    stop("Fichier introuvable : ", f,
         "\n  R exige / ou \\\\ dans les chemins, pas un antislash simple.")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Lecture
# -----------------------------------------------------------------------------
cv  <- read.csv(CSV_CV,  stringsAsFactors = FALSE)
sen <- read.csv(CSV_SEN, stringsAsFactors = FALSE)

# Ponderations en toutes lettres, et ordre fixe pour la legende
cv$Weighting <- factor(cv$Weighting, levels = c("uniform", "distance"),
                       labels = c("Uniform", "Distance"))

best <- cv[cv$Weighting == "Uniform", ]
best <- best[which.min(best$RMSE_SD), ]
cat(sprintf("Minimum de RMSE : ponderation uniforme, k_NN = %d (RMSE = %.4f)\n",
            best$k_NN, best$RMSE_SD))

# -----------------------------------------------------------------------------
# 2. Etiquettes du panneau (b)
#    "DK2O" + 0 h  ->  expression K[2]*O~(0~h), interpretee par plotmath afin
#    que les chiffres passent en indice. Les six cellules sont classees par
#    amplitude relative decroissante.
#    Le tilde tient lieu d'espace insecable dans la syntaxe plotmath.
# -----------------------------------------------------------------------------
to_plotmath <- function(x) sub("\\*$", "", gsub("([0-9]+)", "[\\1]*", x))

sen$label <- sprintf("%s~(%d~h)",
                     to_plotmath(sub("^D", "", sen$Target_variable)),
                     sen$Time_h)
sen$label <- factor(sen$label,
                    levels = sen$label[order(sen$Relative_range_percent)])

# -----------------------------------------------------------------------------
# 3. Titres des axes
#    k porte l'indice NN pour le distinguer sans ambiguite de G, le nombre de
#    groupes de la HCA. Le tiret demi-cadratin est ecrit \u2013 afin que le
#    script reste valide quel que soit l'encodage du fichier.
# -----------------------------------------------------------------------------
XLAB_A <- parse(text = paste0(
  'bold("Number of neighbours (")*bolditalic(k)[bold("NN")]*bold(")")'))

XLAB_B <- parse(text = paste0(
  'bold("Relative range across ")*bolditalic(k)[bold("NN")]',
  '*bold(" = 7\u20139 (%)")'))

# -----------------------------------------------------------------------------
# 4. Theme commun
# -----------------------------------------------------------------------------
theme_S3 <- theme_bw(base_size = 14, base_family = FAM) +
  theme(
    axis.title       = element_text(face = "bold", size = 15),
    axis.text        = element_text(size = 13, colour = "black"),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "grey30", linewidth = 0.6),
    legend.title     = element_text(face = "bold", size = 13),
    legend.text      = element_text(size = 13),
    legend.key       = element_blank(),
    plot.margin      = margin(6, 8, 6, 6)
  )

# -----------------------------------------------------------------------------
# 5. Panneau (a) : validation croisee
#    Le point retenu (uniforme, k_NN = 8) est entoure, sans commentaire ajoute
#    sur le graphe : la legende de figure le precise.
# -----------------------------------------------------------------------------
p_cv <- ggplot(cv, aes(x = k_NN, y = RMSE_SD,
                       colour = Weighting, shape = Weighting)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  geom_point(data = best, shape = 21, size = 4.6, stroke = 1.1,
             colour = "black", fill = NA, show.legend = FALSE) +
  scale_colour_manual(values = c(Uniform = "#1F4E79", Distance = "#B03A2E"),
                      name = "Weighting") +
  scale_shape_manual(values = c(Uniform = 16, Distance = 17),
                     name = "Weighting") +
  scale_x_continuous(breaks = seq(1, max(cv$k_NN), by = 2)) +
  labs(x = XLAB_A, y = "RMSE (standardized units)") +
  theme_S3 +
  theme(legend.position = "inside", legend.position.inside = c(0.22, 0.85),
        legend.background = element_rect(fill = "white", colour = "grey50",
                                         linewidth = 0.4))

# Compatibilite : avant ggplot2 3.5, legend.position recevait directement le
# vecteur de coordonnees.
if (utils::packageVersion("ggplot2") < "3.5.0")
  p_cv <- p_cv + theme(legend.position = c(0.22, 0.85))

# -----------------------------------------------------------------------------
# 6. Panneau (b) : sensibilite locale
# -----------------------------------------------------------------------------
p_sen <- ggplot(sen, aes(x = Relative_range_percent, y = label)) +
  geom_col(fill = "#1F4E79", width = 0.62) +
  geom_text(aes(label = sprintf("%.1f", Relative_range_percent)),
            hjust = -0.25, size = 4.4, family = FAM) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  scale_y_discrete(labels = function(l) parse(text = l)) +
  labs(x = XLAB_B, y = NULL) +
  theme_S3

# -----------------------------------------------------------------------------
# 7. Assemblage, etiquettes (a) et (b)
# -----------------------------------------------------------------------------
p_all <- (p_cv | p_sen) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 18, face = "bold", family = FAM))

print(p_all)

# -----------------------------------------------------------------------------
# 8. Export : TIFF 1200 dpi + PDF vectoriel
# -----------------------------------------------------------------------------
tif <- file.path(OUTDIR, "Figure_S3.tiff")
pdf <- file.path(OUTDIR, "Figure_S3.pdf")

ggsave(tif, p_all, width = W, height = H, units = UNITS, dpi = DPI, compression = "lzw")
ggsave(pdf, p_all, width = W, height = H, units = UNITS, device = cairo_pdf)

mo <- file.info(tif)$size / 1024^2
message("Ecriture dans : ", OUTDIR)
message(sprintf("  Figure_S3.tiff : %6.1f Mo%s", mo,
                if (mo > 10) "  <-- depasse la limite Elsevier de 10 Mo" else ""))
message(sprintf("  Figure_S3.pdf  : %6.1f Mo", file.info(pdf)$size / 1024^2))

# -----------------------------------------------------------------------------
# Reglages fins :
#   MAE au lieu du RMSE dans (a) : remplacer y = RMSE_SD par y = MAE_SD
#   legend.position.inside = c(0.78, 0.85)  -> legende a droite du panneau
#   width = 0.5 dans geom_col               -> barres plus fines
#   W <- 190 ; H <- 80 ; UNITS <- "mm"      -> taille finale d'impression
#   (p_cv / p_sen) au lieu de (p_cv | p_sen) -> panneaux empiles
# =============================================================================
