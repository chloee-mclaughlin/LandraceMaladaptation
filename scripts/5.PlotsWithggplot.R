library(ggplot2)
library(maps) # For map data

# Base map data
world_map <- map_data("world")

# Plot
ggplot() + 
  # Map 
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
              fill = "white", color = "grey15", size = 0.3) +
  # Add vulnerable points
  geom_point(data = temp, aes(x = lon.VulnPix, y = lat.VulnPix), 
             shape = 16, color = "black", size=3) +
  # Add arrows for outline
  geom_segment(data = temp, aes(x = lon.MinVul, y = lat.MinVul, 
                                xend = lon.VulnPix, yend = lat.VulnPix), 
               arrow = arrow(length = unit(0.1, "inches")), 
               color = "black", size = 1.4) +
  # Add arrows colored by Class
  geom_segment(data = temp, aes(x = lon.MinVul, y = lat.MinVul, 
                                xend = lon.VulnPix, yend = lat.VulnPix, 
                                color = Class), 
               arrow = arrow(length = unit(0.1, "inches")), 
               size = 1) +
  #coord_fixed(xlim = c(-115, -35), ylim = c(-40, 35), expand = FALSE) + #Maize
  coord_fixed(xlim = c(70,135), ylim = c(-10,40), expand = FALSE)  + #Rice
  #coord_fixed(xlim = c(-30, 130), ylim = c(-40,45), expand = FALSE)  + #Sorghum
  #coord_fixed(xlim = c(-20, 70), ylim = c(0, 65), expand = FALSE) + #Barley
  
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none") +
  xlab("") +
  ylab("") +
  # Remove axes text
  theme(axis.text = element_blank(),                             
        axis.ticks = element_blank()) +
  # Map the Class colors
  scale_color_identity()  +
  labs(title = "Rice subsp. indica",
       subtitle = "Optimal substitution") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5))  


ggplot() + 
  # Map 
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "grey15", size = 0.3) +
  # Add vulnerable points
  geom_point(data = temp, aes(x = lon.VulnPix, y = lat.VulnPix), 
             shape = 16, color = "black", size=3) +
  # Add arrows for outline
  geom_segment(data = temp, aes(x = lon.Country_MinVul, y = lat.Country_MinVul, 
                                xend = lon.VulnPix, yend = lat.VulnPix), 
               arrow = arrow(length = unit(0.1, "inches")), 
               color = "black", size = 1.4) +
  # Add arrows colored by Class
  geom_segment(data = temp, aes(x = lon.Country_MinVul, y = lat.Country_MinVul, 
                                xend = lon.VulnPix, yend = lat.VulnPix, 
                                color = Class_country), 
               arrow = arrow(length = unit(0.1, "inches")), 
               size = 1) +
  #coord_fixed(xlim = c(-115, -35), ylim = c(-40, 35), expand = FALSE) + #Maize
  coord_fixed(xlim = c(70,135), ylim = c(-10,40), expand = FALSE)  + #Rice
  #coord_fixed(xlim = c(-30, 130), ylim = c(-40,45), expand = FALSE)  + #Sorghum
  #coord_fixed(xlim = c(-20, 70), ylim = c(0, 65), expand = FALSE) + #Barley
  
  
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none") +
  xlab("") +
  ylab("") +
  # Remove axes text
  theme(axis.text = element_blank(),                             
        axis.ticks = element_blank()) +
  # Map the Class colors
  scale_color_identity() +
  labs(title = "Rice subsp. indica",
       subtitle = "Within country best substitution") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5))


