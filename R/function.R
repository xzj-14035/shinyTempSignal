dateType3 <- function(tree, pattern) {
  date <- numeric(length(tree$tip.label))
  for (i in 1:length(tree$tip.label)) {
    date[i] <- str_extract(tree$tip.label[i], pattern)
  }
  return(date)
}
## convert Date to decimal format
dateNumeric <- function(date, format) {
  if (format == "yy" | format == "yyyy") {
    date <- as.numeric(date)
  }
  if (format == "yyyy-MM-dd") {
    date <- Date2decimal(date)
  }
  if (format == "yyyy/MM/dd") {
    date <- as.Date(date, format="%y/%m/%d")
    date <- Date2decimal(date)
  }
  if (format == "yyyy.MM.dd") {
    date <- as.Date(date, format="%y.%m.%d")
    date <- Date2decimal(date)
  }
  if (format == "MM-dd-yyyy") {
    date <- as.Date(date, format="%m-%d-%y")
    date <- Date2decimal(date)
  }
  if (format == "MM/dd/yyyy") {
    date <- as.Date(date, format="%m/%d/%y")
    date <- Date2decimal(date)
  }
  if (format == "MM.dd.yyyy") {
    date <- as.Date(date, format="%m.%d.%y")
    date <- Date2decimal(date)
  }
  return(date)
}

getdivergence <- function(tree, date, method) {
  ## Root a Tree by Root-to-Tip Regression
  #tree2 <- rtt(tree, date, objective=method)
  tree2 <- tree
  divergence <- data.frame(divergence=numeric(length(tree2$tip.label)))
  dataset <- cbind(tree2$edge, tree2$edge.length)
  dataset <- as.data.frame(dataset)
  names(dataset) <- c("from", "to", "length")
  for (i in 1:length(tree2$tip.label)) {
    edge <- nodepath(tree2, from=length(tree2$tip.label) + 1, to=i)
    for (j in 1:(length(edge) - 1)) {
      a <- edge[j]
      b <- edge[j + 1]
      ##The start and end of the node are matched, and the length is assigned
      len <- dataset[which(dataset$from == a & dataset$to == b), ]$length  
      divergence[i, 1] <- divergence[i, 1] + len
    }
  }
  return(divergence)
}



