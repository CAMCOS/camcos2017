########rtsne visualization ###
## Input: data, labels(true labels or lables created by clustering output)
rtsne_viz <- function(data, labels) {
  library(ggplot2)
  library(Rtsne)
	rtsne.out<- Rtsne(data)
	embedding <- as.data.frame(rtsne.out$Y)
	embedding$label <- labels
	# plot(rtsne.out$Y, main = '',col = news100.njw$cluster)
	p <- ggplot(embedding, aes(x=V1, y=V2, col= label)) + geom_point(size = 0.8) +
	  guides(colour = guide_legend(override.aes = list(size=3))) +
	  xlab("") +
	  ylab("") +
	  ggtitle("Rtsne 2D Embedding Visulization") +
	  theme_light(base_size=10) +
	  theme(strip.background = element_blank(), strip.text.x = element_blank() , axis.text.x = element_blank(),
	        axis.text.y = element_blank(),axis.ticks = element_blank(), axis.line = element_blank(), panel.border = element_blank())
	p
}
