
GGTheme <- function() {
  ## theme_set(theme_bw())
  theme_update(strip.text=element_text(size=14),
               strip.text.x=element_text(vjust=1),
               strip.text.y=element_text(vjust=.5, angle=90), #vjust=0
               strip.background=element_rect(color=NA, fill=NA, linetype=0),
               axis.text=element_text(size=12),
               axis.text.y=element_text(angle=90),
               axis.ticks.length = unit(-0.3, "lines"),
               axis.ticks.margin = unit(0.5, "lines"),
               panel.border=element_rect(color=1,fill=NA),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())
}
