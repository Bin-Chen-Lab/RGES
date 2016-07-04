library(gplots)
library(RColorBrewer)

colPal <- redgreen(100)
par(mar=c(12, 4, 2, 0.5))
#image(t(druggable_targets_pos), col=redblue(2))
bins = 30
pdf("../fig/drug_sig.pdf")
a = rbind(sample(1:bins, bins), sample(1:bins, bins), sample(1:bins, bins))
image(a, col= colPal,   axes=F, srt=45)
dev.off()
