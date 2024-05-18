#' @export
niche.overlap <- function(x) {
	# NOTE transpose data.frame (colSums doesn't give correct results!)
	x <- t(x)
	prop.x <- x / rowSums(x)
	# rowSums(prop.x) # should be 1 for all columns!
	index.all <- c(
		# Pianka
		'pianka' = sum(apply(prop.x,2,prod)) / sqrt(prod(rowSums(prop.x^2))),
		# Percentage overlap
		'percentage' = sum(apply(prop.x,2,min)),
		# Simplified Morisita index
		'morisita.simpl' = 2 * sum(apply(prop.x,2,prod))  / sum(rowSums(prop.x^2))
	)

	out <- list(
		'data' = t(x),
		'prop' = t(prop.x),
		'index' = index.all
	)
	out
}
