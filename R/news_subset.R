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
    x %>%
      left_join(r.df, by = c("row" = "old.row")) %>%
      left_join(c.df, by = c("col" = "old.col")) %>%
      left_join(l.df, by = c("label" = "old.lab")) %>%
      select(-row, -col, -label) %>%
      rename(row = new.row, col = new.col, label = new.lab) %>%
      select(row, col, value, everything())
  } else {
    b <- rep(1, nrow(x))
    x %>%
      left_join(r.df, by = c("row" = "old.row")) %>%
      left_join(c.df, by = c("col" = "old.col")) %>%
      left_join(l.df, by = c("label" = "old.lab")) %>%
      select(-row, -col, -label, -value) %>%
      mutate(value = b) %>%
      rename(row = new.row, col = new.col, label = new.lab) %>%
      select(row, col, value, everything())
  }
}

news_subset <- function (X, filter, binary = TRUE, vocabulary = FALSE) {
  library(dplyr)
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

#test1 <- news.subset(newsgroups, "comp")
