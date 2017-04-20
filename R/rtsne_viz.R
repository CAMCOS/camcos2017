########rtsne visualization ###
## Input: data, labels(true labels or lables created by clustering output)
#' Title
#'
#' @param data
#' @param labels
#'
#' @return
#' @export
#'
#' @examples
rtsne_viz <- function(data, labels) {
	rtsne.out<- Rtsne::Rtsne(data)
	embedding <- as.data.frame(rtsne.out$Y)
	embedding$label <- labels
	# plot(rtsne.out$Y, main = '',col = news100.njw$cluster)
	p <- ggplot2::ggplot(embedding, aes(x=V1, y=V2, col= label)) + geom_point(size = 0.8) +
	  ggplot2::guides(colour = guide_legend(override.aes = list(size=3))) +
	  ggplot2::xlab("") +
	  ggplot2::ylab("") +
	  ggplot2::ggtitle("Rtsne 2D Embedding Visulization") +
	  ggplot2::theme_light(base_size=10) +
	  ggplot2::theme(strip.background = element_blank(), strip.text.x = element_blank() , axis.text.x = element_blank(),
	        axis.text.y = element_blank(),axis.ticks = element_blank(), axis.line = element_blank(), panel.border = element_blank())
	p
}
