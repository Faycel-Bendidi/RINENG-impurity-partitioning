# =============================================================================
#  Figure 4 - RINENG-D-26-06081 (Revision R1)
#
#  Une seule image : 10 elements majeurs (5 x 2) au-dessus des 7 elements
#  traces (4 x 2).
#
#  Reponse aux commentaires du Reviewer 5 :
#    - 14 : lisibilite (grille de 17 panneaux reorganisee, panneaux agrandis)
#    -  9 : barres d'erreur (n = 2, +/- 1 ecart-type)
#
#  Reglages : barres en noir, points reduits, caracteres agrandis et en gras,
#             titres avec indices typographiques, trace rompu aux valeurs
#             manquantes, blanc reduit entre chaque titre et son graphe.
# =============================================================================

source(file.path("scripts", "_helpers.R"))
assert_repository_root()
load_required_packages(c("readxl", "ggplot2", "tidyr", "dplyr", "readr", "patchwork"))


# -----------------------------------------------------------------------------
# Chemins
# -----------------------------------------------------------------------------
OUTDIR <- file.path("results", "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
XLSX    <- file.path("data", "external", "Supplementary_Data.xlsx")

if (!file.exists(XLSX))
  stop("Fichier introuvable : ", XLSX,
       "\n  Verifier le chemin relatif depuis la racine du depot.")

# -----------------------------------------------------------------------------
# 1. Lecture du fichier Excel
#    Feuille S2 : 3 blocs empiles (D_i, f, R) ; chaque variable occupe
#    3 colonnes (Test | Mean | SD) et chaque temps occupe 2 lignes, Mean et SD
#    n'etant remplis que sur la premiere des deux.
# -----------------------------------------------------------------------------
TIMES   <- c(0, 12, 24, 36, 48, 60, 72)
HDR_ROW <- 17    # ligne d'en-tete du bloc f ("fP2O5,sol", "fAl2O3,sol", ...)

raw <- as.data.frame(read_excel(XLSX, sheet = "Table S2_coef_mean_SD",
                                col_names = FALSE))

res <- list()
for (cc in seq(2, ncol(raw), by = 3)) {
  nm <- raw[HDR_ROW, cc]
  if (is.na(nm) || nm == "") next
  rows <- HDR_ROW + 2 + 2 * (seq_along(TIMES) - 1)
  res[[length(res) + 1]] <- data.frame(
    Element = as.character(nm),
    t_hours = TIMES,
    Value   = suppressWarnings(as.numeric(raw[rows, cc + 1])),
    SD      = suppressWarnings(as.numeric(raw[rows, cc + 2])),
    stringsAsFactors = FALSE)
}
df_long <- bind_rows(res)

# Retirer UNIQUEMENT les especes entierement sous la limite de detection
# (fCo,sol et fPb,sol). Les valeurs manquantes isolees restent en NA :
# geom_line rompt alors le trace, et rien n'est affiche a ces temps.
df_long <- df_long %>%
  group_by(Element) %>% filter(!all(is.na(Value))) %>% ungroup()

# -----------------------------------------------------------------------------
# 2. Familles et titres des panneaux
#    Syntaxe plotmath : Al[2]*O[3] place 2 et 3 en indice ; bold(...) force le
#    gras, car element_text(face = "bold") n'agit pas sur une expression
#    plotmath.
# -----------------------------------------------------------------------------
major_elements <- c("fAl2O3,sol", "fCaO,sol", "fF,sol", "fFe2O3,sol",
                    "fK2O,sol", "fMgO,sol", "fNa2O,sol", "fP2O5,sol",
                    "fSiO2,sol", "fSO3,sol")
trace_elements <- c("fCd,sol", "fCr,sol", "fCu,sol", "fMn",
                    "fMo,sol", "fNi,sol", "fZn,sol")

lab_map <- c(
  "fP2O5,sol"  = "bold(P[2]*O[5])",
  "fAl2O3,sol" = "bold(Al[2]*O[3])",
  "fCaO,sol"   = "bold(CaO)",
  "fFe2O3,sol" = "bold(Fe[2]*O[3])",
  "fK2O,sol"   = "bold(K[2]*O)",
  "fMgO,sol"   = "bold(MgO)",
  "fNa2O,sol"  = "bold(Na[2]*O)",
  "fSO3,sol"   = "bold(SO[3])",
  "fSiO2,sol"  = "bold(SiO[2])",
  "fF,sol"     = "bold(F)",
  "fCd,sol"    = "bold(Cd)",
  "fCr,sol"    = "bold(Cr)",
  "fCu,sol"    = "bold(Cu)",
  "fMn"        = "bold(Mn)",
  "fMo,sol"    = "bold(Mo)",
  "fNi,sol"    = "bold(Ni)",
  "fZn,sol"    = "bold(Zn)"
)

# Titre de l'axe des y : f avec indice i,sol (i en italique)
YLAB <- expression(bold(italic(f)[list(i, sol)]) ~ bold("(%)"))

# -----------------------------------------------------------------------------
# 3. Fonction de trace, commune aux deux blocs
# -----------------------------------------------------------------------------
make_block <- function(species, couleur, ncol_grid) {
  
  d <- df_long %>% filter(Element %in% species)
  
  # L'ordre des panneaux reste alphabetique sur les noms Excel, comme dans la
  # figure d'origine : on fixe les niveaux du facteur avant de renommer.
  ord <- sort(unique(d$Element))
  d   <- d %>% mutate(Label = factor(lab_map[Element], levels = unname(lab_map[ord])))
  
  ggplot(d, aes(x = t_hours, y = Value, group = Element)) +
    # Barres d'erreur en noir, tracees en premier pour passer sous les points
    geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD),
                  width = 2.5, linewidth = 0.5, colour = "black", na.rm = TRUE) +
    # Trait rompu automatiquement a chaque NA : rien n'est trace la ou la
    # valeur est manquante
    geom_line(linewidth = 0.8, colour = couleur, na.rm = TRUE) +
    geom_point(size = 0.75, colour = couleur, na.rm = TRUE) +
    facet_wrap(~ Label, scales = "free", ncol = ncol_grid,
               labeller = label_parsed) +
    scale_x_continuous(breaks = c(0, 20, 40, 60)) +
    labs(x = "Time (hour)", y = YLAB, title = NULL) +
    theme_bw() +
    theme(
      plot.title = element_blank(),
      
      axis.title = element_text(
        face = "bold",
        size = 20
      ),
      
      # The y-axis title is removed from each block because
      # it is added only once in the final assembly.
      axis.title.y = element_blank(),
      
      axis.text = element_text(
        size = 14,
        colour = "black"
      ),
      
      # Black border around each facet title
      strip.background = element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.7
      ),
      
      strip.text = element_text(
        size = 17,
        margin = margin(t = 1, b = 1)
      ),
      
      # Black border around each individual plotting panel
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.7
      ),
      
      panel.spacing = unit(
        2.5,
        "mm"
      ),
      
      plot.margin = margin(
        2,
        4,
        2,
        2
      ),
      
      panel.grid.major = element_line(
        colour = "grey90"
      ),
      
      panel.grid.minor = element_blank()
    )
}

p_major <- make_block(major_elements, "steelblue", 5)   # 10 panneaux : 5 x 2
p_trace <- make_block(trace_elements, "red",       4)   #  7 panneaux : 4 x 2

# -----------------------------------------------------------------------------
# 4. Assemblage en une seule image, sans etiquette (a) / (b).
#    Le titre de l'axe des x n'apparait qu'une fois, sous le bloc inferieur, et
#    celui de l'axe des y une seule fois a gauche, centre sur toute la hauteur.
# -----------------------------------------------------------------------------
blocs <- (p_major + theme(axis.title.x = element_blank())) / p_trace +
  plot_layout(heights = c(1, 1))

# Titre unique de l'axe des y, pivote de 90 degres
ylab_grob <- grid::textGrob(YLAB, rot = 90, gp = grid::gpar(fontsize = 20))

p_all <- wrap_elements(ylab_grob) | blocs
p_all <- p_all + plot_layout(widths = c(1, 32))   # largeur reservee au titre

print(p_all)

# -----------------------------------------------------------------------------
# 5. Export : TIFF 1200 dpi + PDF vectoriel, dans le meme dossier
# -----------------------------------------------------------------------------
W <- 14.5 ; H <- 10.5 ; UNITS <- "in" ; DPI <- 1200

tif <- file.path(OUTDIR, "Figure_4.tiff")
pdf <- file.path(OUTDIR, "Figure_4.pdf")

ggsave(tif, p_all, width = W, height = H, units = UNITS,
       dpi = DPI, compression = "lzw")
ggsave(pdf, p_all, width = W, height = H, units = UNITS, device = cairo_pdf)

mo <- file.info(tif)$size / 1024^2
message("Ecriture dans : ", OUTDIR)
message(sprintf("  Figure_4.tiff : %6.1f Mo%s", mo,
                if (mo > 10) "  <-- depasse la limite Elsevier de 10 Mo" else ""))
message(sprintf("  Figure_4.pdf  : %6.1f Mo", file.info(pdf)$size / 1024^2))

# -----------------------------------------------------------------------------
# Reglages fins :
#   geom_point(size = 0.6)                             -> points plus petits
#   geom_errorbar(linewidth = 0.65, width = 3.5)       -> barres plus marquees
#   strip.text = element_text(size = 17, margin = margin(t = 0, b = 0))
#                                                      -> titres encore plus
#                                                         serres contre le graphe
#   plot_layout(widths = c(1, 26))                     -> titre de l'axe des y
#                                                         plus eloigne des graphes
#                                                         (augmenter le 1er nombre
#                                                          ou baisser le second)
#
# Taille : 14.5 x 10 pouces = 368 x 254 mm. La largeur utile du journal est de
# 190 mm ; si le TIFF depasse 10 Mo, remplacer la ligne des parametres par
#   W <- 190 ; H <- 131 ; UNITS <- "mm"
# ou televerser le PDF, vectoriel et prefere par Elsevier.
# =============================================================================
