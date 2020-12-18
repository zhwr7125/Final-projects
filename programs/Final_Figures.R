source("D:/Course/BIOS7695/Final projects/programs/7. Dispersion summary.R")

tiff("D:/Course/BIOS7695/Final projects/Figures/Boxplot.tiff",height = 30, width = 20, units="cm",
     compression = "lzw", res = 600)

par(mfrow=c(4,1))
boxplot(seq.counts,xlab ="314 subjects", xaxt = "n",main="Raw counts")
boxplot(counts.between.median,xlab ="314 subjects",main="Normalized counts, method=median", xaxt = "n")
boxplot(counts.between.upper,xlab ="314 subjects",main="Normalized counts, method=upper", xaxt = "n")
boxplot(counts.between.full,xlab ="314 subjects",main="Normalized counts, method=full", xaxt = "n")
dev.off()

tiff("D:/Course/BIOS7695/Final projects/Figures/MAplot.tiff",height = 9, width = 9, units="cm",
     compression = "lzw", res = 400)
meanVarPlot(seq.counts,log=TRUE)
dev.off()

setwd("D:/Course/BIOS7695/Final projects/Figures/Log2change")
tiff(height = 12, width = 8, units="cm",
     compression = "lzw", res = 400)

p1
p2
p3
p4
dev.off()

setwd("D:/Course/BIOS7695/Final projects/Figures/Dispersion")
tiff(height = 8, width = 12, units="cm",
     compression = "lzw", res = 400)
p5
p6
p7
p8
p9
p10
p11
p12
dev.off()
