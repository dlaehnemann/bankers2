library("RColorBrewer")
library("VennDiagram")
venn.plot <- draw.quad.venn( direct.area = TRUE, area.vector = c( 3, 34, 4, 13, 134, 1234, 234, 24, 1, 14, 124, 123, 23, 2, 12 ), 
    category = c( "1", "2", "3", "4" ), 
    lwd = rep(2, 4 ), 
    lty = rep("solid", 4 ), 
    col = c( brewer.pal(4, "Set1" )[1], brewer.pal(4, "Set1" )[2], brewer.pal(4, "Set1" )[3], brewer.pal(4, "Set1" )[4] ), 
    fill = c( brewer.pal(4, "Set1" )[1], brewer.pal(4, "Set1" )[2], brewer.pal(4, "Set1" )[3], brewer.pal(4, "Set1" )[4] ), 
    cat.col = c( brewer.pal(4, "Set1" )[1], brewer.pal(4, "Set1" )[2], brewer.pal(4, "Set1" )[3], brewer.pal(4, "Set1" )[4] ), 
    alpha = rep(0.5, 4 ), 
    label.col = rep("black", 15 ), 
    cex = rep(.9, 15 ), 
    fontface = rep("plain", 15 ), 
    fontfamily = rep("sans", 15 ), 
    cat.pos = c( -15, 15, 0, 0  ), 
    cat.dist = c( 0.22, 0.22, 0.11, 0.11 ), 
    cat.cex = rep(1.5, 4 ), 
    cat.just = rep(list(c(0.5, 0.5)), 4 ), 
    cat.fontface = rep("plain", 4 ), 
    cat.fontfamily = rep("sans", 4 ), 
    rotation.degree = 0, 
    rotation.centre = c(0.5, 0.5  ), 
    print.mode = c("raw", "percent"  ), 
    sigdigs = 2, 
    ind = FALSE, 
    cex.prop = NULL 
)

pdf("vd.test.4.pdf", width = 12, height = 12)
grid.draw(venn.plot)
dev.off()
