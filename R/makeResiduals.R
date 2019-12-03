makeResiduals <- function (x, v, confounder, scaled = TRUE, method = c("res","cc", "pd"))
{
  deletion.method <- match.arg(method) 
  labels <- colnames(x)
  rownames(x) <- 1:nrow(x)
  formeln <- paste0(labels[v], " ~ ", paste(labels[confounder],
                                            collapse = " + "))
  if (any(is.na(x[, c(v, confounder)])))
    x.new <- na.omit(x)
  if (deletion.method == "cc")
    x <- x.new
  daten <- matrix(ncol = length(v), nrow = nrow(x))
  for (node in 1:length(v)) {
    if (deletion.method == "pd") {
      tmp <- lm(as.formula(formeln[node]), data = x)$residuals
      daten[as.numeric(names(tmp)), node] <- tmp
    }
    else {
      daten[, node] <- lm(as.formula(formeln[node]), data = x)$residuals
    }
  }
  colnames(daten) <- labels[v]
  if (scaled == TRUE) {
    daten <- scale(daten)
  }
}
