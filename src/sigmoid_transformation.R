#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')

exp_base <- function(x, base) {
	base^x
}

# Usage:
result <- exp_base(1:5, base=2)

exp_trans <- scales::trans_new(
	name = "exp", 
	transform = function(x) exp_base(x, 1e5),
	inverse = function(x) log(x, base = 1e5))

ggplot(data, aes(x, y)) +
	geom_point() +
	scale_y_continuous(trans = exp_trans())



sigmoid <- function(alpha) {
	1 / (1 + exp(-alpha))
}

inv_sigmoid <- function(alpha) {
	-log((1 / alpha) - 1)
}

sigmoid_trans <- scales::trans_new(
	name = "sigmoid",
	transform = sigmoid,
	inverse = inv_sigmoid
)


numbat_heatmap[[3]] + 
	# scale_alpha_continuous(trans = sigmoid_trans) + 
	# scale_alpha_continuous(trans='log10') + 
	# scale_alpha_continuous(trans = c("log10", "reverse")) + 
	scale_alpha_continuous(trans = exp_trans) +
	NULL


mydf <- data.frame(
	x = seq(-6, 6, length.out = 100),
	y = seq(-6, 6, length.out = 100),
	alpha = seq(-6, 6, length.out = 100)
	)

p1 <- ggplot(mydf, aes(x, y, alpha = alpha)) +
	# stat_function(fun = sigmoid) +
	geom_point() + 
	labs(title = "Linear Scale")

p2 <- p1 + scale_alpha_continuous(trans = sigmoid_trans) +
	labs(title = "Sigmoid Scale")

gridExtra::grid.arrange(p1, p2, ncol = 2)
