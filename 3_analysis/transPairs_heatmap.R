
##############
## HEATMAP ##
#############

library(RColorBrewer)

test3 = test2
test3[test3==0] = NA

z <- seq(range.i[1], range.i[2], length = 101)


pdf(paste("~/data/figures/",id,"_heatmap.pdf", sep = ''), 8,10, useDingbats = F)
par(mar=c(0,0,0,0), omi=(c(0,0,0,0.0)+0.2))

layout(matrix(c(1,2), nrow=2, byrow=F), heights=c(0.2, 0.8), widths=c(0.2,0.8))

scale01 <- function(x, low = min(x), high = max(x)) {
  x <- (x - low)/(high - low)
  x
}

#1 -scale
par(mar=c(2,0,2,22))

image(z = matrix(z), col = colorRampPalette(rev(brewer.pal(n = 7, 
                                                           name = "RdBu")))(100), breaks = z, 
      xaxt = "n", yaxt = "n", axes=F)
box(lwd=1)
lv <- pretty(z)
xv <- scale01(as.numeric(lv),0,range.i[2])
axis(1, at = xv, labels = lv, las=1, font=1, lwd = 1, cex=0.8, cex.axis=1, tck= -0.2, padj=-1)
mtext('Likelihood', 3, 0.5, font = 2, cex = 1.2)


par(mar = c(5,6,2,4))
image(1:ncol(test3), 1:nrow(test3), t(test3), axes=F, 
      col=colorRampPalette(rev(brewer.pal(n = 7, 
                                          name = "RdBu")))(100), breaks=z, xlab="",ylab="")

box(lwd = 2)
abline(v = (1:nrow(test2)) + 0.5,
       h = (1:nrow(test2)) + 0.5,
       col = 'black', lwd = 2, lty = 1)


axis(2, at=1:nrow(test3), labels = rownames(test3), cex.axis = lab.cex, cex=0.1, las=1, tick=0, line=-0.5,font=2)

text(x = 1:nrow(test3), y = line2user(line = 2, side = 1), labels = rownames(test3), cex=lab.cex, xpd = NA, srt = lab.angle, font = 2, adj = c(adj.y,adj.x))

dev.off()


