estimateplots <- function(estimates){
  library(ggplot2)
  library(dplyr)
  
  
  "%||%" <- function(a, b) {
    if (!is.null(a)) a else b
  }
  
  GeomFlatViolin <-
    ggproto("GeomFlatViolin", Geom,
            setup_data = function(data, params) {
              data$width <- data$width %||%
                params$width %||% (resolution(data$x, FALSE) * 0.9)
              
              # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
              data %>%
                group_by(group) %>%
                mutate(ymin = min(y),
                       ymax = max(y),
                       xmin = x,
                       xmax = x + width / 2)
              
            },
            
            draw_group = function(data, panel_scales, coord) {
              # Find the points for the line to go all the way around
              data <- transform(data, xminv = x,
                                xmaxv = x + violinwidth * (xmax - x))
              
              # Make sure it's sorted properly to draw the outline
              newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                               plyr::arrange(transform(data, x = xmaxv), -y))
              
              # Close the polygon: set first and last point the same
              # Needed for coord_polar and such
              newdata <- rbind(newdata, newdata[1,])
              
              ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
            },
            
            draw_key = draw_key_polygon,
            
            default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                              alpha = NA, linetype = "solid"),
            
            required_aes = c("x", "y")
    )
  
  raincloud_theme = theme(
    text = element_text(size = 10),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    legend.title=element_text(size=16),
    legend.text=element_text(size=16),
    legend.position = "right",
    plot.title = element_text(lineheight=.8, face="bold", size = 16),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
  
  
  geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                               position = "dodge", trim = TRUE, scale = "area",
                               show.legend = NA, inherit.aes = TRUE, ...) {
    layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomFlatViolin,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
        trim = trim,
        scale = scale,
        ...
      )
    )
  }
  
  #' @rdname ggplot2-ggproto
  #' @format NULL
  #' @usage NULL
  #' @export
  
  brec <- estimates[[2]]
  Sest <- estimates[[1]]

  bmn <- with(brec, tapply(brec$b, brec$name, FUN = function(x) mean(x, na.rm = T)))
  bmdn <- with(brec, tapply(brec$b, brec$name, FUN = function(x) median(x, na.rm = T)))
  amn <- with(brec, tapply(brec$intercept, brec$name, FUN = function(x) mean(x, na.rm = T)))
  amdn <- with(brec, tapply(brec$intercept, brec$name, FUN = function(x) median(x, na.rm = T)))
  Smdn <- with(Sest, tapply(Sinit, name, FUN = function(x) median(x, na.rm = T)))
  Smn <- with(Sest, tapply(Sinit, name, FUN = function(x) mean(x, na.rm = T)))
  
  jpeg('S estimation.jpg')
  g <- ggplot(data = Sest, aes(y = Sinit, x = name, fill = name)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = Sinit, color = name), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
    expand_limits(x = 5.25) +
    guides(fill = FALSE) +
    guides(color = FALSE) +
    scale_color_brewer(palette = "Set3") +
    scale_fill_brewer(palette = "Set3") +
    # coord_flip() +
    theme_bw() +
    xlab("Raingauge used") + ylab("Sest") +
    geom_hline(yintercept = 1.05)+
    raincloud_theme 
  print(g)
  dev.off()
  
  jpeg('b estimation.jpg')
  g <- ggplot(data = brec, aes(y = b, x = name, fill = name)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = b, color = name), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
    expand_limits(x = 5.25) +
    guides(fill = FALSE) +
    guides(color = FALSE) +
    scale_color_brewer(palette = "Set3") +
    scale_fill_brewer(palette = "Set3") +
    # coord_flip() +
    theme_bw() +
    #geom_hline(yintercept = 3.7)
    xlab("Raingauge used") + ylab("b") +
    raincloud_theme 
  print(g)
  
  dev.off()
  
  jpeg('a estimate.jpg')
  g <- ggplot(data = brec, aes(y = intercept, x = name, fill = name)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = intercept, color = name), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
    expand_limits(x = 5.25) +
    guides(fill = FALSE) +
    guides(color = FALSE) +
    scale_color_brewer(palette = "Set3") +
    scale_fill_brewer(palette = "Set3") +
    # coord_flip() +
    theme_bw() +
    xlab("Raingauge used") + ylab("a") +
    raincloud_theme
  print(g)
  
  dev.off()
  
  return(list(bmn,bmdn,amn,amdn,Smn,Smdn))
}
