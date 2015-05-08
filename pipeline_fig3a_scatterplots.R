if(exists('parent_dir')){
  source(in_parent_dir('init.R'))
}else
  source('init.R')

d = my_fe

keep = apply(d, 1, max) > detect_thresh
d = d[keep,]
fname = "scatterplots.pdf"
pdf(paste('Figure_3a_scatterplots_JB-',date2fname(),'.pdf', sep = ''), width = 12, height = 4)
layout(matrix(1:3, nrow = 1))
for (i in 1:3) {
    xy = d[, c(i , i + 3)]
    
    color_1slope = "#377eb8"
    color_trendLine = "#e41a1c"
    
    #plot(xy, pch = 16, col = rgb(0, 0, 0, 0.1), cex = 0.7, type = "n", xlim = c(0, max(d)), ylim = c(0, max(d)))
    plot(xy, pch = 16, col = rgb(0, 0, 0, 0.1), cex = 0.7, type = "n", xlim = c(min(d), max(d)), ylim = c(min(d), max(d)))
    
    
    
    points(xy, pch = 16, col = rgb(0, 0, 0, 0.1), cex = 0.7)
    
    
    lines(c(detect_thresh, detect_thresh), c(min(d), detect_thresh), col = '#cb181d', lwd = 3)
    lines(c(min(d), detect_thresh), c(detect_thresh, detect_thresh), col = '#cb181d', lwd = 3)
    lines(c(min(d), max(d)), c(min(d), max(d)), lty = 2, col = color_1slope, lwd = 5)
    
    fit <- glm(xy[, 2] ~ xy[, 1])
    co <- coef(fit)
    #abline(fit, col = color_trendLine, lwd = 2, )
    
    
    
    
    r2 = cor(xy)[1, 2]
    r2 = round(r2, digits = 2)
    r2 = format(r2, nsmall = 2)
    #legend("topleft", legend = c("Slope = 1", "Linear Regression"), fill = c(color_1slope, color_trendLine))
    #legend("bottomright", legend = bquote(R^2 == .(rs), list(rs = r2)))
}
dev.off()
print(paste("wrote scatter plots to", fname)) 
