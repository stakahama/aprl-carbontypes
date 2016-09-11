
GGTheme <- function() {
  ## theme_set(theme_bw())

  theme_update(strip.text        = element_text(size=14),
               strip.text.x      = element_text(vjust=1),
               strip.text.y      = element_text(vjust=.5, angle=90), #vjust=0
               strip.background  = element_rect(color=NA, fill=NA, linetype=0),
               axis.title        = element_text(size=14),
               axis.text         = element_text(size=12),
               axis.text.x       = element_text(margin=margin(.5, .5, .5, .5, "lines")),
               axis.text.y       = element_text(margin=margin(.5, .5, .5, .5, "lines")),
               #axis.text.y      = element_text(angle=90, hjust=.5),
               axis.ticks.length = unit(-0.3, "lines"),
               legend.background = element_rect(fill="white"),
               legend.key        = element_rect(fill="white"),
               panel.background  = element_rect(fill="white"),
               panel.border      = element_rect(color=1, fill=NA),
               panel.grid.major  = element_blank(),
               panel.grid.minor  = element_blank())
}


colors.FG <- c(
  "ketone"                        =rgb(119, 204, 110, max=255), # "darkkhaki"      ,
  "aldehyde"                      ="darkolivegreen2",
  "COOH"                          ="green",
  "aCOH"                          ="deeppink",
  "phenol"                        ="deeppink4",
  "ether"                         ="lightpink",
  "ester"                         ="darksalmon",
  "anhydride"                     ="darkorchid1",
  "peroxide"                      ="slategray3",
  "hydroperoxide"                 ="slategray2",
  "peroxy acid"                   ="gray48",
  "CONO2"                         ="brown",
  "peroxyacyl nitrate"            ="red",
  "nitro"                         ="darkorange2",
  "aCH"                           ="blue",
  "eCH"                           ="cornflowerblue",
  "rCH"                           ="midnightblue",
  "formaldehyde"                  ="darkolivegreen2",
  "formic acid"                   ="green",
  "peroxy nitrate"                ="gray",
  "oxy radical"                   ="gray",
  "C=O-O"                         ="gray",
  "peroxy acid"                   ="snow3",
  "CO"                            ="#008080"
)

colors.OSC <- setNames(tail(colorRampPalette(c("darkorange", "lightgray", "darkblue"))(9), -1),
                       seq(-4, 3))

labels.FG <- c(
  "alkane CH"="aCH",
  "alcohol"="aCOH",
  "peroxyacylnitrate"="peroxyacyl nitrate",
  "carboxylic acid"="COOH",
  "organonitrate"="CONO2",
  "alkene CH" = "eCH",
  "aromatic CH" = "rCH",
  ##
  "carbonylperoxyacid"="carbonyl peroxy acid",
  "carbonylperoxyacid radical"="carbonyl peroxy acid (*)",
  "carboxylic radical"="carboxyl (*)",
  "oxy radical"="alkoxyl (*)",
  "peroxy radical"="peroxyl (*)",
  "C non quaternary"="tertiary sp2 carbon",
  "C=O-O group no H"="R2C=O-O",
  "C=O-O group single H"="RHC=O-O",
  "C=O-O group two H"="H2C=O-O"
)

labels.method <- c( # export_lambdaC.R, production_nC_tseries.R, figures_nC_reconstruction.R
  "count"="COUNT",
  "solve"="COMPOUND",
  "fit"="MIXTURE",
  "nominal"="NOMINAL"
)

Relabel <- function(x, new)
  ifelse(x %in% names(new), new[x], x)

GGColorHue <- function(n) {
  ## http://stackoverflow.com/a/8197703/143476
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

Capitalize <- function(x)
  `substring<-`(x, 1, 1, toupper(substring(x, 1, 1)))

ExpandLim <- function(x, e=.03)
  x + e*diff(x)*c(-1,1)


## from http://stackoverflow.com/a/25735086/143476
letters.greek <-
  structure(list(grsym = c("α", "β", "γ", "δ", "ε", "ζ",
                           "η", "θ", "ι", "κ", "λ", "μ", "ν", "ξ", "ο", "π", "ρ",
                           "ς", "σ", "τ", "υ", "φ", "χ", "ψ", "ω", "Α", "Β", "Γ",
                           "Δ", "Ε", "Ζ", "Η", "Θ", "Ι", "Κ", "Λ", "Μ", "Ν", "Ξ",
                           "Ο", "Π", "Ρ", "Σ", "Τ", "Υ", "Φ", "Χ", "Ψ", "Ω"),
                 decUTF = c(945, 946, 947, 948, 949, 950, 951, 952, 953, 954,
                            955, 956, 957, 958, 959, 960, 961, 962, 963, 964, 965, 966,
                            967, 968, 969, 913, 914, 915, 916, 917, 918, 919, 920, 921,
                            922, 923, 924, 925, 926, 927, 928, 929, 931, 932, 933, 934,
                            935, 936, 937), hexUTF = structure(c(945L, 946L, 947L, 948L,
                                                                 949L, 950L, 951L, 952L, 953L, 954L, 955L, 956L, 957L, 958L,
                                                                 959L, 960L, 961L, 962L, 963L, 964L, 965L, 966L, 967L, 968L,
                                                                 969L, 913L, 914L, 915L, 916L, 917L, 918L, 919L, 920L, 921L,
                                                                 922L, 923L, 924L, 925L, 926L, 927L, 928L, 929L, 931L, 932L,
                                                                 933L, 934L, 935L, 936L, 937L), class = "hexmode"), htmlSym = c("&alpha;",
                                                                                                                                "&beta;", "&gamma;", "&delta;", "&epsilon;", "&zeta;", "&eta;",
                                                                                                                                "&theta;", "&iota;", "&kappa;", "&lambda;", "&mu;", "&nu;",
                                                                                                                                "&xi;", "&omicron;", "&pi;", "&rho;", "&sigmaf;", "&sigma;",
                                                                                                                                "&tau;", "&upsilon;", "&phi;", "&chi;", "&psi;", "&omega;",
                                                                                                                                "&Alpha;", "&Beta;", "&Gamma;", "&Delta;", "&Epsilon;", "&Zeta;",
                                                                                                                                "&Eta;", "&Theta;", "&Iota;", "&Kappa;", "&Lambda;", "&Mu;",
                                                                                                                                "&Nu;", "&Xi;", "&Omicron;", "&Pi;", "&Rho;", "&Sigma;",
                                                                                                                                "&Tau;", "&Upsilon;", "&Phi;", "&Chi;", "&Psi;", "&Omega;"
                                                                                                                                ), Description = c("GREEK SMALL LETTER ALPHA", "GREEK SMALL LETTER BETA",
                                                                                                                                                   "GREEK SMALL LETTER GAMMA", "GREEK SMALL LETTER DELTA", "GREEK SMALL LETTER EPSILON",
                                                                                                                                                   "GREEK SMALL LETTER ZETA", "GREEK SMALL LETTER ETA", "GREEK SMALL LETTER THETA",
                                                                                                                                                   "GREEK SMALL LETTER IOTA", "GREEK SMALL LETTER KAPPA", "GREEK SMALL LETTER LAMBDA",
                                                                                                                                                   "GREEK SMALL LETTER MU", "GREEK SMALL LETTER NU", "GREEK SMALL LETTER XI",
                                                                                                                                                   "GREEK SMALL LETTER OMICRON", "GREEK SMALL LETTER PI", "GREEK SMALL LETTER RHO",
                                                                                                                                                   "GREEK SMALL LETTER FINAL SIGMA", "GREEK SMALL LETTER SIGMA",
                                                                                                                                                   "GREEK SMALL LETTER TAU", "GREEK SMALL LETTER UPSILON", "GREEK SMALL LETTER PHI",
                                                                                                                                                   "GREEK SMALL LETTER CHI", "GREEK SMALL LETTER PSI", "GREEK SMALL LETTER OMEGA",
                                                                                                                                                   "GREEK CAPITAL LETTER ALPHA", "GREEK CAPITAL LETTER BETA",
                                                                                                                                                   "GREEK CAPITAL LETTER GAMMA", "GREEK CAPITAL LETTER DELTA",
                                                                                                                                                   "GREEK CAPITAL LETTER EPSILON", "GREEK CAPITAL LETTER ZETA",
                                                                                                                                                   "GREEK CAPITAL LETTER ETA", "GREEK CAPITAL LETTER THETA",
                                                                                                                                                   "GREEK CAPITAL LETTER IOTA", "GREEK CAPITAL LETTER KAPPA",
                                                                                                                                                   "GREEK CAPITAL LETTER LAMBDA", "GREEK CAPITAL LETTER MU",
                                                                                                                                                   "GREEK CAPITAL LETTER NU", "GREEK CAPITAL LETTER XI", "GREEK CAPITAL LETTER OMICRON",
                                                                                                                                                   "GREEK CAPITAL LETTER PI", "GREEK CAPITAL LETTER RHO", "GREEK CAPITAL LETTER SIGMA",
                                                                                                                                                   "GREEK CAPITAL LETTER TAU", "GREEK CAPITAL LETTER UPSILON",
                                                                                                                                                   "GREEK CAPITAL LETTER PHI", "GREEK CAPITAL LETTER CHI", "GREEK CAPITAL LETTER PSI",
                                                                                                                                                   "GREEK CAPITAL LETTER OMEGA")), .Names = c("grsym", "decUTF",
                                                                                                                                                                                              "hexUTF", "htmlSym", "Description"), row.names = c(NA, -49L), class = "data.frame")
