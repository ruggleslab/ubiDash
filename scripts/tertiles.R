# ++++++++++++++++++++++++++++
# tertiles
# ++++++++++++++++++++++++++++
tertiles <- function(x) {
  vTert <- quantile(x, c(0:3 / 3), na.rm = TRUE)
  cut(
    x,
    vTert,
    include.lowest = TRUE,
    labels = c("Low", "Medium", "High")
  )

}
