#!/home/chialanghsu/miniconda2/bin/Rscript
options(stringsAsFactors = FALSE)
# library(optparse)
library(dplyr)
library(magrittr)
library(drc)

option_list = list(
	make_option(c("-f", "--infile"), type="character", default=NULL,
		help="input file for sequenza", metavar="character"),
	make_option(c("-o", "--outputDir"), type="character", default=NULL,
        help="path to output directory", metavar="character"),
	make_option(c("-s", "--sampleId"), type="character", default=NULL,
        help="sample id", metavar="character"),
	make_option(c("-b", "--blacklist"), type="character", default=NULL,
		help="path to blackpoistion file", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# strip trailing / from all parameters in opt
for(i in 1:length(opt)){
        if(substr(opt[i], nchar(opt[i]), nchar(opt[i])) == '/'){
                opt[i] <- substr(opt[i], 1, nchar(opt[i])-1)
        }
}

print(opt)

##########################
# Command line arguments #
##########################

infile    <- opt$infile
sampleId  <- opt$sampleId
outdir    <- opt$outputDir
blackFile <- opt$blacklist

if (is.null(infile)){
	print_help(opt_parser)
	stop("Missing arguments.\n", call.=FALSE)
}

##########################
# Runing analysis        #
##########################
my.blue <- "#4D82B8"
my.blue.light <- "#B0D7FE"
my.red <- "#D86344"
my.red.light <- "#FEC8B9"
my.grey <- "#F8F2F2"
pt.col <- c(Yes = my.red, No = "grey99")
#b1_exon <- list(E9 = c(1100, 1233),
#				E8 = c(1704, 1901),
#				E7 = c(1981, 2059),
#				E6 = c(2459, 2625),
#				E5 = c(2988, 3142),
#				E4 = c(3429, 3632),
#				E3 = c(3770, 3969),
#				E2 = c(5779, 5934),
#				E1 = c(6322, 6400))

b1_exon <- list(E6 = c(-100, 27),
				E5 = c(390, 544),
				E4 = c(831, 1034),
				E3 = c(1172, 1371),
				E2 = c(3181, 3336),
				E1 = c(3724, 3969))

#b2_exon <- list(E9 = c(1100, 1535),
#				E8 = c(1972, 2169),
#				E7 = c(2249, 2327),
#				E6 = c(2627, 2893),
#				E5 = c(3706, 3860),
#				E4 = c(4147, 4350),
#				E3 = c(4488, 4687),
#				E2 = c(6501, 6656),
#				E1 = c(7044, 7200))

b2_exon <- list(E6 = c(-500, -486),
				E5 = c(327, 481),
				E4 = c(768, 971),
				E3 = c(1109, 1308),
				E2 = c(3122, 3277),
				E1 = c(3665, 3906))

##########################
# Separation Function    #
##########################
separate <- function(data, gene, plot=TRUE){
  # ret_list <- list()
  if (gene == 'CYP11B1') {
    gloc <- 142875852
    # fitting_model <- W1.4()
    plot_location <- "topright"
  }else if (gene == "CYP11B2"){
    gloc <- 142913938
    # fitting_model <- W2.4()
    plot_location <- "topleft"

  }
  x <- data$Pos
  y <- data$Rate + 0.0000001
  model <- drm(y ~ x, fct = L.4(), data = data, separate = TRUE)

  if (plot){
    plot(model, log="", main=gene, add=TRUE, type="none", col="grey30")
  }

	# pseudo.x <- seq(from=0,to=4000,by=10)
  pred <- predict(model, newdata = data.frame(x))

  dif <- vector()
  for (i in 1:length(pred)){
    if (i != length(pred)){
      d <- abs((pred[i+1]-pred[i]) / (x[i+1]-x[i]))
      dif <- c(dif, d)
    }
  }
	max_pos <- x[which(dif == max(dif))]
  max_pos_1 <- x[which(dif == max(dif))+1]

  sep_position <- max_pos + gloc
  sep_position_1 <- max_pos_1 + gloc

  if (plot) {
    abline(v=(max_pos+max_pos_1)/2, lty = 2, col = "grey30")
    legend(plot_location, paste('chr8:', sep_position, "-", sep_position_1, sep=""),
			bty = "n", cex = 0.7, inset = c(-0.01, 0.01))
  }
}


## black-positions will be masked
blackList <- read.table(blackFile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

indata <- read.table(infile, sep = "\t", header = TRUE, check.names = FALSE)
indata$BlackPos <- blackList$BlackPos
fig.out <- file.path(outdir, paste(sampleId, "_4k_plot.png", sep = ""))
png(fig.out, width = 12, height = 12, units = "cm", res = 600, pointsize = 12)
par(mfrow = c(2, 1), mar = c(3, 4, 2, 2))
gene <- "CYP11B1"
bg <- filter(indata, Gene == gene, Major == "No", BlackPos == "No")
fg <- filter(indata, Gene == gene, Major == "Yes", BlackPos == "No")

bg.threshold <- quantile(bg$Rate, 0.995)
fg.threshold <- quantile(fg$Rate[fg$Rate > bg.threshold], 0.5)

bg %$% plot(Pos, Rate, col = "grey90", pch = 16,
		xlab = "", ylab = "Mismatch rate", main = sampleId,
		ylim = c(0, 1), xaxt = "n", yaxt = "n")

abline(h = c(bg.threshold), lty = 2, col = "grey50")

fg %$% points(Pos, Rate, col = ifelse(Count > 20, my.red, my.red.light), pch = 16)

axis(2, las = 2)
mtext(gene, side = 1, line = 1.8)
## add fitting LINE
fit_line <- lm(fg$Rate~fg$Pos)
slope <- abs(coef(fit_line)[2])
if (abs(slope*10000)>1){
	separate(data=fg, gene=gene, plot=TRUE)
}

par(xpd = TRUE)
lines(c(27, 3724), c(-0.15, -0.15) ,lwd = 2, col = "#FF9478")
for (ex in names(b1_exon)){
	rect(b1_exon[[ex]][1], -0.2, b1_exon[[ex]][2], -0.1, col = "#FF9478", border = NA)
	text((b1_exon[[ex]][1] + b1_exon[[ex]][2])/2, -0.3, ex, cex = 0.65)
}


par(xpd = FALSE, mar = c(4.5, 4, 0.5, 2))
gene <- "CYP11B2"
bg <- filter(indata, Gene == gene, Major == "No", BlackPos == "No")
fg <- filter(indata, Gene == gene, Major == "Yes", BlackPos == "No")

bg.threshold <- quantile(bg$Rate, 0.995)
fg.threshold <- quantile(fg$Rate[fg$Rate > bg.threshold], 0.5)

bg %$% plot(Pos, Rate, col = "grey90", pch = 16,
		xlab = "", ylab = "Mismatch rate", main = '',
		ylim = c(0, 1), xaxt = "n", yaxt = "n")

abline(h = c(bg.threshold), lty = 2, col = "grey50")

fg %$% points(Pos, Rate, col = ifelse(Count > 20, my.blue, my.blue.light), pch = 16)

axis(2, las = 2)
mtext(gene, side = 1, line = 1.8)

## add fitting LINE
fit_line <- lm(fg$Rate~fg$Pos)
slope <- abs(coef(fit_line)[2])
if (abs(slope*10000)>1){
	separate(data=fg, gene=gene, plot=TRUE)
}


par(xpd = TRUE)
lines(c(-486, 3665), c(-0.15, -0.15) ,lwd = 2, col = "#60A1E3")
#points()
for (ex in names(b2_exon)){
	rect(b2_exon[[ex]][1], -0.2, b2_exon[[ex]][2], -0.1, col = "#60A1E3", border = NA)
	text((b2_exon[[ex]][1] + b2_exon[[ex]][2])/2, -0.3, ex, cex = 0.65)
}

dev.off()
