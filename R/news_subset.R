#' Select a subset of the full 20Newsgroupds dataset
#'
#' Label names: comp.graphics comp.os.ms-windows.misc comp.sys.ibm.pc.hardware
#' comp.sys.mac.hardware comp.windows.x rec.autos rec.motorcycles
#' rec.sport.baseball rec.sport.hockey sci.crypt sci.electronics sci.med
#' sci.space misc.forsale talk.politics.misc talk.politics.guns
#' talk.politics.mideast talk.religion.misc alt.atheism soc.religion.christian
#'
#' @param X The newsgroups data object loaded with data(newsgroups).
#' @param filter Either an integer vector specifying the label numbers, or a
#'   character vector specifying the label names. Supports "starts with" partial
#'   matching for label names.
#' @param binary Logical. Make the data values 1's
#' @param vocabulary Logical. Include a word column in the returned data frame
#'   in list element one.
#'
#' @return List of 2. (1) The remapped data. (2) The distinct labels
#'   corresponding to each row.
#' @import dplyr
#' @export
news_subset <- function (X, filter, binary = TRUE, vocabulary = FALSE) {

  #require(dplyr, quietly = TRUE, warn.conflicts = FALSE)

  if (is.numeric(filter)) {
    df <- X %>% filter(label %in% filter)
  } else if (is.character(filter)) {
    data(map)
    tf <- as.logical(rowSums(sapply(filter, startsWith, x = map$subject)))
    df <- X %>% filter(label %in% map[tf,2])
  } else {
    stop("Invalid filter argument. Must be integer vector of label numbers or label name (prefix partial matching supported).")
  }

  if (vocabulary) {
    data(vocab)
    df <- df %>% left_join(vocab, by = "col")
  }
  .remap(df, binary = binary)
}

.remap <- function(x, binary = TRUE) {
  r <- sort(unique(x$row))
  r.df <- data.frame(old.row = r, new.row = 1:length(r))
  c <- sort(unique(x$col))
  c.df <- data.frame(old.col = c, new.col = 1:length(c))
  l <- sort(unique(x$label))
  l.df <- data.frame(old.lab = l, new.lab = 1:length(l))

  cat("Label Mapping: ", fill = TRUE)
  print.data.frame(l.df, row.names = FALSE)

  if (!binary) {
    bx <- x %>%
      left_join(r.df, by = c("row" = "old.row")) %>%
      left_join(c.df, by = c("col" = "old.col")) %>%
      left_join(l.df, by = c("label" = "old.lab")) %>%
      select(-row, -col, -label) %>%
      rename(row = new.row, col = new.col, label = new.lab) %>%
      select(row, col, value, everything())

    labels <- bx %>% distinct(row, label) %>% select(label) %>% unlist(use.names = FALSE)
    list(data = bx, labels = labels)
  } else {
    b <- rep(1, nrow(x))
    fx <- x %>%
      left_join(r.df, by = c("row" = "old.row")) %>%
      left_join(c.df, by = c("col" = "old.col")) %>%
      left_join(l.df, by = c("label" = "old.lab")) %>%
      select(-row, -col, -label, -value) %>%
      mutate(value = b) %>%
      rename(row = new.row, col = new.col, label = new.lab) %>%
      select(row, col, value, everything())

    labels <- fx %>% distinct(row, label) %>% select(label) %>% unlist(use.names = FALSE)
    list(data = fx, labels = labels)
  }

}
