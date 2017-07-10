#
# title     : sourceAuxFunctions.R
# purpose   : R script to source all required auxiliary functions
# author    : WD41, LRS/EPFL/PSI
# date      : July 2017
#
source("./r-scripts/GetTimeExpTC.R")
source("./r-scripts/GetTimeExpCO.R")
source("./r-scripts/GetTimeIdxTrc.R")
source("./r-scripts/FlattenTimeIdxExp.R")

source("./r-scripts/CreateInputPCScores.R")
source("./r-scripts/CreateInputBias.R")

source("./r-scripts/CalcPCScores.R")

source("./r-scripts/CalcAveVec.R")
source("./r-scripts/CalcKrigingVarMat.R")
source("./r-scripts/CalcTruncationVarMat.R")
source("./r-scripts/CalcBiasVarMat.R")
source("./r-scripts/CalcLogLikelihood.R")

source("./r-scripts/findpeaks.R")
source("./r-scripts/kneedle.R")
source("./r-scripts/dhalfcauchy.R")