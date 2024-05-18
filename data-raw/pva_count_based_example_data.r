# PVA DATASETS

# Yellowstone grizzly bear data (From Morris & Doak)
grizzly <- ts(
	c(44,47,46,44,46,45,46,40,39,39,42,39,41,
	40,33,36,34,39,35,34,38,36,37,41,39,51,
	47,57,48,60,65,74,69,65,57,70,81,99,99),
	start=1959)

# Adder data, Mills Box 5.5, p. 110
adder <- ts(
	c(19,18,22,25,23,23,NA,NA,19,17,15,9,6,7,4),
	start=1981)

# Cormorant population in Denmark
cormorant <- ts(
	c(2170,2490,4046,4984,6232,7789,9655,12760,14310,19270,23610,29500,33840,36630,
	37880,38200,40990,36370,38850,39790,42270,39500,40750,37350,40150,40160,37700,
	35850,34310,33090),
	start=1980)

# Peregrines in Sweden (number of breeding pairs in SW Sweden 1972-2007)
peregrine <- ts(
	c(5,4,3,2,1,3,3,3,3,3,2,2,2,2,3,3,3,4,6,9,10,14,
	16,21,23,28,29,28,28,34,39,39,49,49,47,56),
	start=1972)

# Red-cockaded woodpecker (two time series from southern US)
woodpecker_florida <- ts(
	c(65,59,55,38,38,28,32,35,34,27,26,26),
	start=1981)

woodpecker_ncarolina <- ts(
	c(463,419,399,392,433,405,407,460,436,416,415),
	start=1981)

usethis::use_data(grizzly)
usethis::use_data(adder)
usethis::use_data(cormorant)
usethis::use_data(peregrine)
usethis::use_data(woodpecker_florida)
usethis::use_data(woodpecker_ncarolina)
