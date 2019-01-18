

## displayResults.R


## functions that display results

# Define ggplot2 theme
gg_theme <- function(){
    p <- theme(
        plot.title = element_text(size = 20,face = 'bold',
                                  margin = margin(0,0,3,0), hjust = 0.5),
        axis.text = element_text(size = rel(1.05), color = 'black'),
        axis.title = element_text(size = rel(1.45), color = 'black'),
        axis.title.y = element_text(margin = margin(0,10,0,0)),
        axis.title.x = element_text(margin = margin(10,0,0,0)),
        axis.ticks.x = element_line(colour = "black", size = rel(0.8)),
        axis.ticks.y = element_blank(),
        legend.position = "right",
        legend.key.size = unit(1.4, 'lines'),
        legend.title = element_text(size = 12, face = 'bold'),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        panel.grid.major = element_line(colour = "gainsboro"),
        panel.background = element_blank()
    )
    return(p)
}

# end of displayResults.R









