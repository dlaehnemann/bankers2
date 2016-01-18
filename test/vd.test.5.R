library("RColorBrewer")
library("VennDiagram")
venn.plot <- draw.quintuple.venn( direct.area = TRUE, area.vector = c( 1, 2, 3, 4, 5, 35, 15, 14, 12, 25, 23, 13, 34, 24, 45, 345, 135, 145, 124, 125, 235, 123, 134, 234, 245, 2345, 1345, 1245, 1235, 1234, 12345 ), 
    category = c( "1", "2", "3", "4", "5" ), 
    lwd = rep(2, 5 ), 
    lty = rep("solid", 5 ), 
    col = c( brewer.pal(5, "Set1" )[1], brewer.pal(5, "Set1" )[2], brewer.pal(5, "Set1" )[3], brewer.pal(5, "Set1" )[4], brewer.pal(5, "Set1" )[5] ), 
    fill = c( brewer.pal(5, "Set1" )[1], brewer.pal(5, "Set1" )[2], brewer.pal(5, "Set1" )[3], brewer.pal(5, "Set1" )[4], brewer.pal(5, "Set1" )[5] ), 
    cat.col = c( brewer.pal(5, "Set1" )[1], brewer.pal(5, "Set1" )[2], brewer.pal(5, "Set1" )[3], brewer.pal(5, "Set1" )[4], brewer.pal(5, "Set1" )[5] ), 
    alpha = rep(0.5, 5 ), 
    label.col = rep("black", 31 ), 
    cex = rep(.9, 31 ), 
    fontface = rep("plain", 31 ), 
    fontfamily = rep("sans", 31 ), 
    cat.pos = c( 0, 287.5, 215, 145, 70  ), 
    cat.dist = rep( 0.1, 5 ), 
    cat.cex = rep(1.5, 5 ), 
    cat.just = rep(list(c(0.5, 0.5)), 5 ), 
    cat.fontface = rep("plain", 5 ), 
    cat.fontfamily = rep("sans", 5 ), 
    rotation.degree = 0, 
    rotation.centre = c(0.5, 0.5  ), 
    print.mode = c("raw", "percent"  ), 
    sigdigs = 2, 
    ind = FALSE, 
    cex.prop = NULL 
)

pdf("vd.test.5.pdf", width = 12, height = 12)
grid.draw(venn.plot)
dev.off()
