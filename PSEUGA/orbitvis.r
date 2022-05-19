library(plotly)

points <-read.csv("orbit.csv", header = TRUE)

xpoints <- points$x
ypoints <- points$y
zpoints <- points$z


#scatter3d(x = xpoints, y = ypoints, z = zpoints, surface = FALSE, grid = TRUE, ) #needs car library

fig <- plot_ly(points, x = ~x, y = ~y, z = ~z, color = ~group, colors = c('#FF0000','#FFAAAA', '#00FF00','#AAFFAA', '#0000FF','#AAAAFF', '#000000'))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'xpos', range = list(-2,2)),
                                   yaxis = list(title = 'ypos', range = list(-2,2)),
                                   zaxis = list(title = 'zpos', range = list(-2,2)),
                                   camera = list(projection = list(type = 'orthographic'))))


fig


