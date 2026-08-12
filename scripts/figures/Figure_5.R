# =============================================================================
#  Figure 5 - RINENG-D-26-06081 (Revision R1)
#
#  Une seule image : 9 elements majeurs (5 x 2) au-dessus des 7 elements
#  traces (4 x 2). Mise en forme strictement identique a la Figure 4.
#
#  Grandeur tracee : R_liq,i / P2O5 (Eq. 5), rapport de la concentration de
#  l'element i a celle de P2O5 dans la phase liquide.
#
#  Reponse aux commentaires du Reviewer 5 :
#    - 14 : lisibilite (grille de 16 panneaux reorganisee, panneaux agrandis)
#    -  9 : barres d'erreur (n = 2, +/- 1 ecart-type)
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
#    HDR_ROW = 33 : ligne d'en-tete du bloc R ("RliqP2O5/P2O5", ...).
# -----------------------------------------------------------------------------
TIMES   <- c(0, 12, 24, 36, 48, 60, 72)
HDR_ROW <- 33

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

# Retirer les especes entierement sous la limite de detection dans le liquide
# (Co et Pb). Le bloc R ne contient aucune autre valeur manquante.
df_long <- df_long %>%
  group_by(Element) %>% filter(!all(is.na(Value))) %>% ungroup()

# Retirer RliqP2O5/P2O5, egal a 1 a tous les temps par construction.
df_long <- df_long %>% filter(Element != "RliqP2O5/P2O5")

# -----------------------------------------------------------------------------
# 2. Familles et titres des panneaux
#    Les elements traces sont rapportes multiplies par 10^4 dans le fichier
#    source ; le facteur est porte par le titre du panneau, l'axe des y restant
#    commun a toute la figure.
# -----------------------------------------------------------------------------
major_elements <- c("RliqAl2O3/P2O5", "RliqCaO/P2O5", "RliqF/P2O5",
                    "RliqFe2O3/P2O5", "RliqK2O/P2O5", "RliqMgO/P2O5",
                    "RliqNa2O/P2O5", "RliqSiO2/P2O5", "RliqSO3/P2O5")
trace_elements <- c("104 x RliqCd/P2O5", "104 x RliqCr/P2O5", "104 x RliqCu/P2O5",
                    "104 x RliqMn/P2O5", "104 x RliqMo/P2O5", "104 x RliqNi/P2O5",
                    "104 x RliqZn/P2O5")

lab_map <- c(
  "RliqAl2O3/P2O5"    = "bold(Al[2]*O[3])",
  "RliqCaO/P2O5"      = "bold(CaO)",
  "RliqF/P2O5"        = "bold(F)",
  "RliqFe2O3/P2O5"    = "bold(Fe[2]*O[3])",
  "RliqK2O/P2O5"      = "bold(K[2]*O)",
  "RliqMgO/P2O5"      = "bold(MgO)",
  "RliqNa2O/P2O5"     = "bold(Na[2]*O)",
  "RliqSiO2/P2O5"     = "bold(SiO[2])",
  "RliqSO3/P2O5"      = "bold(SO[3])",
  "104 x RliqCd/P2O5" = "bold(Cd %*% 10^4)",
  "104 x RliqCr/P2O5" = "bold(Cr %*% 10^4)",
  "104 x RliqCu/P2O5" = "bold(Cu %*% 10^4)",
  "104 x RliqMn/P2O5" = "bold(Mn %*% 10^4)",
  "104 x RliqMo/P2O5" = "bold(Mo %*% 10^4)",
  "104 x RliqNi/P2O5" = "bold(Ni %*% 10^4)",
  "104 x RliqZn/P2O5" = "bold(Zn %*% 10^4)"
)

# Titre de l'axe des y : R avec "liq" en exposant et "i/P2O5" en indice, les
# 2 et 5 etant eux-memes en indice a l'interieur de cet indice (Eq. 5).
# En plotmath, x[indice]^{exposant} combine les deux, et les crochets se
# nichent : P[2] a l'interieur de l'indice donne un indice de second niveau.
YLAB <- expression(bold(italic(R)[italic(i)/P[2]*O[5]]^{"liq"}))

# Ordre alphabetique sur le symbole de l'element, et non sur le nom Excel :
# sans cela "104 x ..." regrouperait mal les panneaux traces.
sym <- function(v) {
  s <- sub("^104 x ", "", v)
  s <- sub("^Rliq", "", s)
  sub("/P2O5$", "", s)
}

# -----------------------------------------------------------------------------
# 3. Fonction de trace, commune aux deux blocs
# -----------------------------------------------------------------------------
make_block <- function(species, couleur, ncol_grid) {
  
  d   <- df_long %>% filter(Element %in% species)
  ord <- species[order(sym(species))]
  d   <- d %>% mutate(Label = factor(lab_map[Element], levels = unname(lab_map[ord])))
  
  ggplot(d, aes(x = t_hours, y = Value, group = Element)) +
    # Barres d'erreur en noir, tracees en premier pour passer sous les points
    geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD),
                  width = 2.5, linewidth = 0.5, colour = "black", na.rm = TRUE) +
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
      
      axis.title.y = element_blank(),
      
      axis.text = element_text(
        size = 14,
        colour = "black"
      ),
      
      strip.background = element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.7
      ),
      
      strip.text = element_text(
        size = 17,
        margin = margin(t = 1, b = 1)
      ),
      
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

p_major <- make_block(major_elements, "steelblue", 5)   # 9 panneaux : 5 x 2
p_trace <- make_block(trace_elements, "red",       4)   # 7 panneaux : 4 x 2

# -----------------------------------------------------------------------------
# 4. Assemblage en une seule image.
#    Le titre de l'axe des x n'apparait qu'une fois, sous le bloc inferieur, et
#    celui de l'axe des y une seule fois a gauche, centre sur toute la hauteur.
# -----------------------------------------------------------------------------
blocs <- (p_major + theme(axis.title.x = element_blank())) / p_trace +
  plot_layout(heights = c(1, 1))

ylab_grob <- grid::textGrob(YLAB, rot = 90, gp = grid::gpar(fontsize = 20))

p_all <- wrap_elements(ylab_grob) | blocs
p_all <- p_all + plot_layout(widths = c(1, 32))   # largeur reservee au titre

print(p_all)

# -----------------------------------------------------------------------------
# 5. Export : TIFF 1200 dpi + PDF vectoriel, dans le meme dossier
# -----------------------------------------------------------------------------
W <- 14.5 ; H <- 10.0 ; UNITS <- "in" ; DPI <- 1200

tif <- file.path(OUTDIR, "Figure_5.tiff")
pdf <- file.path(OUTDIR, "Figure_5.pdf")

ggsave(tif, p_all, width = W, height = H, units = UNITS,
       dpi = DPI, compression = "lzw")
ggsave(pdf, p_all, width = W, height = H, units = UNITS, device = cairo_pdf)

mo <- file.info(tif)$size / 1024^2
message("Ecriture dans : ", OUTDIR)
message(sprintf("  Figure_5.tiff : %6.1f Mo%s", mo,
                if (mo > 10) "  <-- depasse la limite Elsevier de 10 Mo" else ""))
message(sprintf("  Figure_5.pdf  : %6.1f Mo", file.info(pdf)$size / 1024^2))

# -----------------------------------------------------------------------------
# Reglages fins :
#   geom_point(size = 0.6)                             -> points plus petits
#   geom_errorbar(linewidth = 0.65, width = 3.5)       -> barres plus marquees
#   strip.text = element_text(size = 17, margin = margin(t = 0, b = 0))
#                                                      -> titres plus serres
#   plot_layout(widths = c(1, 26))                     -> titre de l'axe des y
#                                                         plus eloigne des graphes
#
# Taille : 14.5 x 10 pouces = 368 x 254 mm. La largeur utile du journal est de
# 190 mm ; si le TIFF depasse 10 Mo, remplacer la ligne des parametres par
#   W <- 190 ; H <- 131 ; UNITS <- "mm"
# ou televerser le PDF, vectoriel et prefere par Elsevier.
# =============================================================================
