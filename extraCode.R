# Solutions 


# Extract all "core" taxa
core <- occ_abun$otu[occ_abun$fill == 'core']


# Plotting OTUs above and below the natural model prediction

ap = obs.np$freq > (obs.np$pred.upr)
bp = obs.np$freq < (obs.np$pred.lwr)

aboveOTUs=rownames(obs.np)[ap == TRUE]
aboveOTUs=rownames(obs.np)[obs.np$freq > (obs.np$pred.upr)]

aboveOTUs <- aboveOTUs[!is.na(aboveOTUs)]
belowOTUs=rownames(obs.np)[bp == TRUE]
belowOTUs <- belowOTUs[!is.na(belowOTUs)]

ggplot() +
  geom_point(data=occ_abun[occ_abun$otu %in% aboveOTUs,], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='blue', size=2)+
  geom_point(data=occ_abun[occ_abun$otu %in% belowOTUs,], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='darkred', size=2) +
  geom_point(data=occ_abun[!(occ_abun$otu %in% c(aboveOTUs,belowOTUs)),], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='white', size=1) +
  geom_line(color='black', data=obs.np, size=1, aes(y=obs.np$freq.pred, x=log10(obs.np$p)), alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.upr, x=log10(obs.np$p)), alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.lwr, x=log10(obs.np$p)), alpha=.25)+
  labs(x="log10(mean relative abundance)", y="Occupancy") +
  annotate(
    "text", x=-5, y=.8, label=paste("Above prediction: ", round(above.pred, digits = 2), sep=''),
    color="blue", size=4) +
  annotate(
    "text", x=-2.5, y=.1, label=paste("Below prediction: ", round(below.pred, digits = 2), sep=''),
    color="darkred", size=4)
