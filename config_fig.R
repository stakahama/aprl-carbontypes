
GGTheme <- function() {
  ## theme_set(theme_bw())

  theme_update(strip.text=element_text(size=14),
               strip.text.x=element_text(vjust=1),
               strip.text.y=element_text(vjust=.5, angle=90), #vjust=0
               strip.background=element_rect(color=NA, fill=NA, linetype=0),
               axis.title=element_text(size=14),
               axis.text=element_text(size=12),
               axis.text.x=element_text(margin=margin(.5, .5, .5, .5, "lines")),
               axis.text.y=element_text(margin=margin(.5, .5, .5, .5, "lines")),
               #axis.text.y=element_text(angle=90, hjust=.5),
               axis.ticks.length = unit(-0.3, "lines"),
               panel.background=element_rect(color="white"),
               panel.border=element_rect(color=1, fill=NA),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())
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
  "peroxy acid"                   ="snow3"
)

labels.FG <- c(
  "alkane CH"="aCH",
  "alcohol"="aCOH",
  "peroxyacylnitrate"="peroxyacyl nitrate",
  "carboxylic acid"="COOH",
  "organonitrate"="CONO2",
  "alkene CH" = "eCH",
  "aromatic CH" = "rCH"
)

Relabel <- function(x, new)
  ifelse(x %in% names(new), new[x], x)

GGColorHue <- function(n) {
  ## http://stackoverflow.com/a/8197703/143476
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

