library("RColorBrewer")
library("VennDiagram")
venn.plot <- draw.triple.venn( direct.area = TRUE, area.vector = c( 1, 12, 2, 13, 123, 23, 3 ), 
    category = c( "1", "2", "3" ), 
    lwd = rep(2, 3 ), 
    lty = rep("solid", 3 ), 
    col = c( brewer.pal(3, "Set1" )[1], brewer.pal(3, "Set1" )[2], brewer.pal(3, "Set1" )[3] ), 
    fill = c( brewer.pal(3, "Set1" )[1], brewer.pal(3, "Set1" )[2], brewer.pal(3, "Set1" )[3] ), 
    cat.col = c( brewer.pal(3, "Set1" )[1], brewer.pal(3, "Set1" )[2], brewer.pal(3, "Set1" )[3] ), 
    alpha = rep(0.5, 3 ), 
    label.col = rep("black", 7 ), 
    cex = rep(.9, 7 ), 
    fontface = rep("plain", 7 ), 
    fontfamily = rep("sans", 7 ), 
    cat.pos = c( -40, 40, 180  ), 
    cat.dist = c( 0.05, 0.05, 0.025 ), 
    cat.cex = rep(1.5, 3 ), 
    cat.just = rep(list(c(0.5, 0.5)), 3 ), 
    cat.fontface = rep("plain", 3 ), 
    cat.fontfamily = rep("sans", 3 ), 
    rotation.degree = 0, 
    rotation.centre = c(0.5, 0.5  ), 
    print.mode = c("raw", "percent"  ), 
    sigdigs = 2, 
    ind = FALSE, 
    cex.prop = NULL 
)

pdf("vd.test.3.pdf", width = 12, height = 12)
grid.draw(venn.plot)
dev.off()
