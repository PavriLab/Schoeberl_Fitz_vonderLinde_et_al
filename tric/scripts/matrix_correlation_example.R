library(hicrep)
library(corrplot)
args = commandArgs(trailingOnly=TRUE)

if (args[1] == 'run') {

	mat1 = read.table('../TriCplots/TriC12_lane_1_B18AIDKO2_d0_1+3_TriC_interactions_1000_RAW.tab')
	mat2 = read.table('../TriCplots/TriC12_lane_1_B18AIDKO2_d0_2_TriC_interactions_1000_RAW.tab')
	mat3 = read.table('../TriCplots/TriC12_lane_1_B18AIDKO3_d0_1+3_TriC_interactions_1000_RAW.tab')
	mat4 = read.table('../TriCplots/TriC12_lane_1_B18AIDKO3_d0_2_TriC_interactions_1000_RAW.tab')

	#treatments:
	mat2d2_1 = read.table('../TriCplots/TriC12_lane_1_B18AIDKO2_d2_1_TriC_interactions_1000_RAW.tab')
	mat2d2_2 = read.table('../TriCplots/TriC12_lane_1_B18AIDKO2_d2_2_TriC_interactions_1000_RAW.tab')
	mat2d2_3 = read.table('../TriCplots/TriC12_lane_1_B18AIDKO2_d2_3_TriC_interactions_1000_RAW.tab')
	mat2d3_1 = read.table('../TriCplots/TriC12_lane_1_B18AIDKO2_d3_1_TriC_interactions_1000_RAW.tab')
	mat2d3_2 = read.table('../TriCplots/TriC12_lane_1_B18AIDKO2_d3_2_TriC_interactions_1000_RAW.tab')
	mat2d3_3 = read.table('../TriCplots/TriC12_lane_1_B18AIDKO2_d3_3_TriC_interactions_1000_RAW.tab')
	mat3d2_1 = read.table('../TriCplots/TriC12_lane_1_B18AIDKO3_d2_1_TriC_interactions_1000_RAW.tab')
	mat3d2_2 = read.table('../TriCplots/TriC12_lane_1_B18AIDKO3_d2_2_TriC_interactions_1000_RAW.tab')
	mat3d2_3 = read.table('../TriCplots/TriC12_lane_1_B18AIDKO3_d2_3_TriC_interactions_1000_RAW.tab')
	mat3d3_1 = read.table('../TriCplots/TriC12_lane_1_B18AIDKO3_d3_1_TriC_interactions_1000_RAW.tab')
	mat3d3_2 = read.table('../TriCplots/TriC12_lane_1_B18AIDKO3_d3_2_TriC_interactions_1000_RAW.tab')
	mat3d3_3 = read.table('../TriCplots/TriC12_lane_1_B18AIDKO3_d3_3_TriC_interactions_1000_RAW.tab')

	tablefu <- function(mat){
	print('tablefu called')
	range1 = seq(114435000, 114668000, by = 1000) 
	range2 = seq(114436000, 114669000, by = 1000)
	matA = cbind(range2, mat)
	matA = cbind(range1, matA)
	matA = cbind('chr12', matA)
	return(matA)
	}

	correlator <- function(matA, matB, h){
	matA = tablefu(matA)
	matB = tablefu(matB)
	processed <- prep(matA, matB, 1000, h, 234000)
	scc.out = get.scc(processed, resol = 1000, max = 234000)
	print('corr successful')
	return(scc.out$scc)
	# scc.out$std
	}

	# trained = htrain(tablefu(mat1), tablefu(mat2), 1000, 234000, 0:20)

	controls = list(mat1, mat2, mat3, mat4, mat2d2_1, mat2d2_2, mat2d2_3, mat2d3_1, mat2d3_2, mat2d3_3, mat3d2_1, mat3d2_2, mat3d2_3, mat3d3_1, mat3d3_2, mat3d3_3)
	corrtable = matrix(nrow = length(controls), ncol = length(controls))
	for (x in seq(length(controls)-1)) {
	print(x)
	for (y in seq(x+1, length(controls))) {
		print(y)
		cor = correlator(controls[[x]], controls[[y]], 12)
		corrtable[x, y] = cor
	}
	}
	controlnames = c('Exp_1_rep_1_0h', 'Exp_1_rep_2_0h', 'Exp_2_rep_1_0h', 'Exp_2_rep_2_0h', 
	'Exp_1_rep_1_48h', 'Exp_1_rep_2_48h', 'Exp_1_rep_3_48h', 'Exp_1_rep_1_72h', 'Exp_1_rep_2_72h', 'Exp_1_rep_3_72h',
	'Exp_2_rep_1_48h', 'Exp_2_rep_2_48h', 'Exp_2_rep_3_48h', 'Exp_2_rep_1_72h', 'Exp_2_rep_2_72h', 'Exp_2_rep_3_72h')
	corrdf = as.data.frame(corrtable)
	colnames(corrdf) = controlnames
	rownames(corrdf) = controlnames
	corrdf = round(corrdf, 3)

	makeSymm <- function(m) {
	m[lower.tri(m)] <- t(m)[lower.tri(m)]
	return(m)
	}

	corrdf[is.na(corrdf)] <- 0
	corrdf = makeSymm(corrdf)
} else {
	corrdf <- read.csv("Mouse_TriC12_IgG1pool_control_corr.csv", row.names=1, header=TRUE)
}

corrdf = data.matrix(corrdf)

pdf("mouse_TriC12_IgG1pool_corr_plot.pdf")
corrplot(corrdf, type = 'upper', method = 'color', number.cex = 0.8, tl.col = 'black', addCoef.col = 'black', diag = FALSE, mar=c(0,0,1,0), title = 'B18 AIDKO IgG1 pool capture')
dev.off()

write.csv(corrdf, file='Mouse_TriC12_IgG1pool_control_corr.csv')
