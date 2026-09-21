library(patRoon)

# --------------------------------------------------------- #
#                        ReadMe                          ####
#
# Please ensure that the patRoon script is run in a separate
# RProject and that no other R packages are loaded to avoid
# conflicts between packages.
#
# More information on patRoon can be found in the handbook:
# https://rickhelmus.github.io/patRoon/handbook_bd/index.html
#
# or the patRoon tutorial:
# https://rickhelmus.github.io/patRoon/articles/tutorial.html
#
# --------------------------------------------------------- #




# --------------------------- #
# initialization           ####
# --------------------------- #

workPath <- "Working_Directory_Path"
setwd(workPath)

# Load analysis table
anaInfoPos <- read.csv("samples_pos.csv")
anaInfoNeg <- read.csv("samples_neg.csv")

# --------------------------- #
# features                 ####
# --------------------------- #

# Find all features
fListPos <- findFeatures(anaInfoPos, "openms",
                         noiseThrInt = 2250, chromSNR = 5,
                         chromFWHM = 10, minFWHM = 1,
                         maxFWHM = 20)

fListNeg <- findFeatures(anaInfoNeg, "openms",
                         noiseThrInt = 2250, chromSNR = 5,
                         chromFWHM = 10, minFWHM = 1,
                         maxFWHM = 20)

# Group and align features between analyses
fGroups_pos <- groupFeatures(fListPos, "openms", rtalign = FALSE)
fGroups_neg <- groupFeatures(fListNeg, "openms", rtalign = FALSE)

# Basic rule based filtering
fGroups_pos <- filter(fGroups_pos,
                      preAbsMinIntensity = 100, absMinIntensity = 10000,
                      relMinReplicateAbundance = 1,
                      maxReplicateIntRSD = 0.75,
                      blankThreshold = 5,
                      removeBlanks = TRUE,
                      retentionRange = c(60, Inf),
                      mzRange = NULL)

fGroups_neg <- filter(fGroups_neg,
                      preAbsMinIntensity = 100, absMinIntensity = 10000,
                      relMinReplicateAbundance = 1,
                      maxReplicateIntRSD = 0.75,
                      blankThreshold = 5,
                      removeBlanks = TRUE,
                      retentionRange = c(60, Inf),
                      mzRange = NULL)

# --------------------------- #
# componentization         ####
# --------------------------- #

# Perform automatic generation of components
components_pos <- generateComponents(fGroups_pos, "cliquems", ionization = "positive", ppm = 7) 
components_neg <- generateComponents(fGroups_neg, "cliquems", ionization = "negative", ppm = 7) 


# Only keep results with adducts/isotope annotations and remove exhotic adducts
components_filt_pos <- delete(components_pos, j = function(ct, ...) is.na(ct$isonr) & (is.na(ct$adduct_ion) | !ct$adduct_ion %in% c("[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+H]+")))

components_filt_neg <- delete(components_neg, j = function(ct, ...) is.na(ct$isonr) & (is.na(ct$adduct_ion) | !ct$adduct_ion %in% c("[M-H2O]-", "[M-H]-")))

# Only keep components with at least one feature group also present in the suspect screened data
components_filt_pos <- filter(components_filt_pos, size = c(2, 1E6)) 
components_filt_neg <- filter(components_filt_neg, size = c(2, 1E6))


# manually check/filter components 
checkComponents(components_filt_pos, fGroups_pos, session = "components_check_pos.yml")
# components_filt_pos <- filter(components_filt_pos, checkComponentsSession = "components_check_pos.yml")

checkComponents(components_filt_neg, fGroups_neg, session = "components_check_neg.yml")
# components_filt_neg <- filter(components_filt_neg, checkComponentsSession = "components_check_neg.yml")

# to filter out unwanted adducts/isotopes 
fGroups_pos <- selectIons(fGroups_pos, components_filt_pos, c("[M+H]+"))
fGroups_neg <- selectIons(fGroups_neg, components_filt_neg, c("[M-H]-"))

# Make sets and group features 
fGroups <- makeSet(fGroups_pos, fGroups_neg, groupAlgo = "openms", adducts = c("[M+H]+", "[M-H]-"))


# --------------------------- #
# suspect screening        ####
# --------------------------- #

# Load suspect list
suspList <- read.csv("Path_To_Suspect_List|_File",
                     stringsAsFactors = FALSE)

# Set onlyHits to FALSE to retain features without suspects (eg for full NTA)
fGroups_suspects <- screenSuspects(fGroups, suspList, rtWindow = 12, mzWindow = 0.005, onlyHits = TRUE)

# --------------------------- #
# transformation products  ####
# --------------------------- #

# Load parent suspect list
suspListParents <- read.csv("Suspect_List_File_Path.csv",
                            stringsAsFactors = FALSE)

# Obtain TPs
TPs <- generateTPs("biotransformer", parents = suspListParents, type = "env", generations = 3, calcSims = FALSE)
TPs <- filter(TPs, removeParentIsomers = TRUE)
TPs <- delete(TPs, j = function(tab, par)
{
  parNM <- parents(TPs)[match(par, name)]$neutralMass # parent neutral mass
  # remove if high generation (>1) and mass is <60% of that of the parent
  # generation 1 TPs are always kept and the rest only if their mass is sufficiently high (>=60%)
  return(tab$generation > 1 & (tab$neutralMass / parNM) < 0.6)
})
print(TPs)

# Screen TPs
suspListTPs <- convertToSuspects(TPs, includeParents = FALSE)

fGroups_suspects <- screenSuspects(fGroups_suspects, suspListTPs, rtWindow = 12, mzWindow = 0.005, onlyHits = TRUE, amend = TRUE)

# --------------------------- #
# annotation               ####
# --------------------------- #

# Retrieve MS peak lists
avgMSListParams <- getDefAvgPListParams(clusterMzWindow = 0.005)
mslists <- generateMSPeakLists(fGroups_suspects, "mzr", maxMSRtWindow = 5,
                               precursorMzWindow = 0.5,
                               avgFeatParams = avgMSListParams,
                               avgFGroupParams = avgMSListParams)

# Rule based filtering of MS peak lists. You may want to tweak this. See the manual for more information.
mslists <- filter(mslists, absMSIntThr = NULL,
                  absMSMSIntThr = NULL,
                  relMSIntThr = NULL,
                  relMSMSIntThr = 0.05,
                  topMSPeaks = NULL, topMSMSPeaks = 25)

# Calculate formula candidates
formulas <- generateFormulas(fGroups_suspects, mslists, "genform",
                             relMzDev = 5, elements = "CHNOPSClBr",
                             oc = FALSE, calculateFeatures = FALSE,
                             featThresholdAnn = 0.75, 
                             setThresholdAnn = 0,
                             timeout = 3
)

convertToMFDB(TPs, "TP-database.csv", includeParents = TRUE)

# Calculate compound structure candidates
compounds <- generateCompounds(fGroups_suspects, mslists, "metfrag", dbRelMzDev = 5, fragRelMzDev = 5, fragAbsMzDev = 0.002,
                               database = "csv",
                               extraOpts = list(LocalDatabasePath = "TP-database.csv"),
                               maxCandidatesToStop = 2500, setThresholdAnn = 0)
compounds <- addFormulaScoring(compounds, formulas, updateScore = TRUE)

# Annotate suspects
fGroups <- annotateSuspects(fGroups_suspects, formulas = formulas, compounds = compounds, MSPeakLists = mslists,
                            IDFile = "idlevelrules.yml")
fGroups_suspects <- annotateSuspects(fGroups_suspects, formulas = formulas, compounds = compounds, MSPeakLists = mslists,
                            IDFile = "idlevelrules.yml")


# --------------------------- #
# Filters to remove low ID ####
# --------------------------- #

fGroups_suspects_id3 <- filter(fGroups, maxLevel = 3, onlyHits = TRUE) # remove ID of 4 and 5
fGroups_suspects_id4 <- filter(fGroups, maxLevel = 4, onlyHits = TRUE) # remove ID of 5

# --------------------------- #
# Parent and TP linkage    ####
# --------------------------- #

# You probably want to prioritize the data before componentization. Please see the handbook for more info.
componentsTP <- generateComponents(fGroups_suspects[, suspects = suspList$name], "tp", fGroupsTPs = fGroups_suspects, TPs = TPs, MSPeakLists = mslists,
                                   formulas = formulas, compounds = compounds)

# You may want to configure the filtering step below. See the manuals for more details.
componentsTP <- filter(componentsTP, retDirMatch = TRUE, minSpecSim = NULL, minSpecSimPrec = NULL,
                       minSpecSimBoth = NULL, minFragMatches = NULL, minNLMatches = NULL)

# Only keep linked feature groups
fGroups_linked <- fGroups_suspects[results = componentsTP] # Filter for TP-precursor match
fGroups_unlinked <- fGroups_suspects[, setdiff(names(fGroups_suspects), names(fGroups_linked))] # Filter for no TP-precursor match


# --------------------------- #
# reporting                ####
# --------------------------- #

options(future.globals.maxSize = 2e9)

# Advanced report settings can be edited in the report.yml file.
report(fGroups, MSPeakLists = NULL, formulas = NULL,
       compounds = NULL, components = componentsTP, TPs = TPs,
       settingsFile = "report.yml", openReport = TRUE)
