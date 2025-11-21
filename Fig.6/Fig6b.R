#create treatment figure

df <- data.frame(tx = rep(c("once treated", "start between", 
                            "treatment between", "within treatment", 
                            "end between", "untreated"), times=2), timepoint = c(rep("sample 1", times=6),rep("sample 2", times=6)))

Fig6b <- df |> 
  mutate(timepoint = as.numeric(timepoint),
         tx = as.numeric(tx)) |> 
  ggplot(aes(x=timepoint, y=tx))+
  labs(x="", y="", title="Treatment contexts")+
  
  # tx = 3: within treatment
  geom_segment(aes(x = 0.5, y = 6, xend = 3, yend = 6), 
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"), color = "red", size = 0.5)+
  geom_segment(aes(x = 0, y = 6, xend = 0.5, yend = 6), 
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"), color = "black", size = 0.5)+
  
  # tx = 2: untreated
  geom_segment(aes(x = 0, y = 5, xend = 3, yend = 5), 
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"), color = "black", size = 0.5)+

  # tx = 4: treatment between
  geom_segment(aes(x = 1.8, y = 4, xend = 3, yend = 4), 
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"), color = "black", size = 0.5)+
  geom_segment(aes(x = 1.2, y = 4, xend = 1.8, yend = 4), 
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"), color = "red", size = 0.5)+
  geom_segment(aes(x = 0, y = 4, xend = 1.2, yend = 4), 
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"), color = "black", size = 0.5)+

  # tx = 5: start between
  geom_segment(aes(x = 1.5, y = 3, xend = 3, yend = 3), 
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"), color = "red", size = 0.5)+
  geom_segment(aes(x = 0, y = 3, xend = 1.5, yend = 3), 
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"), color = "black", size = 0.5)+
  
  # tx = 6: once treated
  geom_segment(aes(x = 0, y = 2, xend = 0.8, yend = 2), 
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"), color = "red", size = 0.5)+
  geom_segment(aes(x = 0.8, y = 2, xend = 3, yend = 2), 
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"), color = "black", size = 0.5)+
  
  # tx = 1: end between
  geom_segment(aes(x = 1.5, y = 1, xend = 3, yend = 1), 
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"), color = "black", size = 0.5)+
  geom_segment(aes(x = 0, y = 1, xend = 1.5, yend = 1), 
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"), color = "red", size = 0.5)+
  
  #sampling timepoints
  geom_segment(aes(x = 1, y = 0.5, xend = 1, yend = 6.5), color = "grey", linetype = "dashed", size = 0.5)+
  geom_segment(aes(x = 2, y = 0.5, xend = 2, yend = 6.5), color = "grey", linetype = "dashed", size = 0.5)+
  
  scale_x_continuous(breaks = c(1, 2), labels = c("First sample", "Second sample"))+
  scale_y_continuous(breaks = c(seq(1,6,1)), labels = sort(unique(df$tx)))+
  theme(axis.text.y = element_text(color = "black", size = 8), 
        axis.text.x = element_text(angle=45, hjust=1, vjust=1, color = "black", size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        panel.border = element_rect(color = 'black', 
                                    fill = NA, 
                                    size = 1))
Fig6b + theme(plot.margin = margin(t=2, l=1, r=1, b=1, unit = "cm"))
Fig6b