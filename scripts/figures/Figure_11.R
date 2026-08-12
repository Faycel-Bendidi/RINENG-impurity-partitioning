# =============================================================================
#  Figure 11 - RINENG-D-26-06081 (Revision R1)
#  Dendrogramme de la classification hierarchique (HCA) des profils Z_i(t)
#
#  Source : partition_coefficients_standardized_Z.csv
#           (coefficients de partage centres-reduits, Eq. 6)
#
#  Methode : distance euclidienne sur les 17 profils cinetiques, chacun decrit
#            par ses 7 valeurs temporelles, et agregation de Ward (ward.D2).
# =============================================================================

source(file.path("scripts", "_helpers.R"))
assert_repository_root()
load_required_packages(c("dendextend"))


# Times New Roman. Sous Windows la police est disponible d'office ; l'alias
# "serif" du peripherique graphique y renvoie. windowsFonts() la declare
# explicitement pour que le TIFF et le PDF l'utilisent aussi.
if (.Platform$OS.type == "windows") {
  windowsFonts(Times = windowsFont("Times New Roman"))
  FAM <- "Times"
} else {
  FAM <- "serif"
}

# -----------------------------------------------------------------------------
# Chemins et parametres
# -----------------------------------------------------------------------------
CSV     <- file.path("data", "processed", "partition_coefficients_standardized_Z.csv")
OUTDIR  <- file.path("results", "figures")

G       <- 3        # nombre de groupes encadres (voir la note en fin de script)
W       <- 9.0 ; H <- 6.0 ; UNITS <- "in" ; DPI <- 1200

LWD       <- 2.4    # epaisseur des branches du dendrogramme
LWD_BOX   <- 1.3    # epaisseur des cadres de groupes (plus fine que LWD)
LAB_GAP   <- 0.03   # ecart entre le bas d'une branche et son etiquette
LAB_SPACE <- 0.42   # hauteur reservee aux etiquettes SOUS la ligne y = 0,
# en fraction de la hauteur du dendrogramme. Les cadres
# descendent jusque-la, de sorte que les etiquettes soient
# entierement contenues dans leur encadre.

# Couleurs des cadres, dans l'ordre des groupes de gauche a droite
COLS    <- c("#E41A1C", "#4DAF4A", "#370EB8", "#00BFFF", "#984EA3")

if (!file.exists(CSV))
  stop("Fichier introuvable : ", CSV,
       "\n  R exige / ou \\\\ dans les chemins, pas un antislash simple.")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Lecture et classification
#    Le fichier a les temps en lignes et les especes en colonnes ; les objets
#    a classer etant les profils cinetiques, la matrice est transposee.
# -----------------------------------------------------------------------------
d <- read.csv(CSV, row.names = 1, check.names = FALSE)
V <- t(as.matrix(d))                      # 17 especes x 7 temps

dst <- dist(V, method = "euclidean")
hc  <- hclust(dst, method = "ward.D2")

cat(sprintf("Correlation cophenetique : %.3f\n", cor(dst, cophenetic(hc))))
grp <- cutree(hc, k = G)
cat("Groupes retenus (G =", G, ") :\n")
for (g in seq_len(G))
  cat(sprintf("  C%d : %s\n", g,
              paste(sub("^Z", "", names(grp)[grp == g]), collapse = ", ")))

# -----------------------------------------------------------------------------
# 2. Etiquettes avec indices typographiques
#    Le prefixe "Z" est retire et les chiffres passes en indice :
#    "ZAl2O3" -> expression Al[2]*O[3].
# -----------------------------------------------------------------------------
sym  <- sub("^Z", "", hc$labels)
expr <- parse(text = sub("\\*$", "", gsub("([0-9]+)", "[\\1]*", sym)))

# -----------------------------------------------------------------------------
# 3. Trace
# -----------------------------------------------------------------------------
draw <- function() {
  op <- par(mar = c(1.5, 5, 1, 1), family = FAM, lwd = LWD)
  on.exit(par(op))
  
  n     <- length(sym)
  ymax  <- max(hc$height)
  ybot  <- -LAB_SPACE * ymax          # bas de la zone reservee aux etiquettes
  
  # plot.hclust calcule ses propres bornes verticales et ignore l'argument
  # ylim : tout ce qui descend sous zero y serait rogne. On passe donc par
  # as.dendrogram(), dont la methode plot accepte xlim et ylim explicitement.
  # leaflab = "none" supprime les etiquettes, reecrites ensuite en plotmath
  # pour obtenir les indices.
  dend <- as.dendrogram(hc)
  
  plot(dend, ylim = c(ybot, ymax), xlim = c(0.5, n + 0.5),
       leaflab = "none", axes = FALSE, yaxt = "n",
       ylab = "Height", cex.lab = 1.5, font.lab = 2,
       edgePar = list(lwd = LWD, col = "black"))
  
  axis(2, at = pretty(c(0, ymax)), cex.axis = 1.2, lwd = LWD)
  
  # Cadres colores. rect.hclust() les ferait descendre jusqu'au bas de la zone
  # de trace sans tenir compte des etiquettes, et n'accepte pas d'argument
  # lwd. Les rectangles sont donc traces directement, avec la meme regle de
  # hauteur que rect.hclust -- mi-hauteur entre les deux dernieres fusions
  # conservees, commune aux G cadres -- et un bas place sous les etiquettes.
  grp_ord <- grp[hc$order]
  eff     <- table(grp_ord)[unique(grp_ord)]       # effectifs, ordre des feuilles
  bornes  <- c(0, cumsum(eff))
  ytop    <- mean(rev(hc$height)[(G - 1):G])
  
  for (i in seq_len(G))
    rect(bornes[i] + 0.66, ybot, bornes[i + 1] + 0.33, ytop,
         border = COLS[i], lwd = LWD_BOX)
  
  # Etiquettes verticales, sous la ligne y = 0 ou s'alignent les feuilles
  ord <- hc$order
  text(x = seq_along(ord), y = -LAB_GAP * ymax,
       labels = expr[ord], srt = 90, adj = 1, cex = 1.15, family = FAM)
}

# Controle : si un appel echoue, le trace resterait incomplet sans message.
tryCatch(draw(),
         error = function(e) stop("Le trace a echoue : ", conditionMessage(e)))

# -----------------------------------------------------------------------------
# 4. Export : TIFF 1200 dpi + PDF vectoriel
# -----------------------------------------------------------------------------
tif <- file.path(OUTDIR, "Figure_11.tiff")
pdf <- file.path(OUTDIR, "Figure_11.pdf")

tiff(tif, width = W, height = H, units = UNITS, res = DPI, compression = "lzw")
draw(); dev.off()

cairo_pdf(pdf, width = W, height = H, family = FAM)
draw(); dev.off()

mo <- file.info(tif)$size / 1024^2
message("Ecriture dans : ", OUTDIR)
message(sprintf("  Figure_11.tiff : %6.1f Mo%s", mo,
                if (mo > 10) "  <-- depasse la limite Elsevier de 10 Mo" else ""))
message(sprintf("  Figure_11.pdf  : %6.1f Mo", file.info(pdf)$size / 1024^2))

# -----------------------------------------------------------------------------
# NOTE SUR LE NOMBRE DE GROUPES
#
# Largeurs de silhouette moyennes obtenues sur ces donnees :
#   G = 2 : 0.501   G = 3 : 0.446   G = 4 : 0.425   G = 5 : 0.426
# La valeur G = 3 est celle retenue dans la Table S5. A G = 4, le quatrieme
# groupe se reduit a Mo seul, ce qui n'apporte aucune information de
# regroupement. Passer G <- 4 ci-dessus pour le verifier.
#
# Reglages fins :
#   cex = 1.3 dans text()        -> etiquettes plus grandes
#   LAB_SPACE <- 0.40            -> cadres plus bas, si une etiquette depasse
#   LAB_GAP <- 0.05              -> etiquettes plus eloignees des branches
#   LWD <- 3 ; LWD_BOX <- 1.6    -> traits plus epais
#   W <- 190 ; H <- 127 ; UNITS <- "mm"   -> taille finale d'impression
# =============================================================================
