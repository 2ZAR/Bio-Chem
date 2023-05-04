#Practice#
#Extened Connectivity Finger Print#

> install.packages("rcdk")
> install.packages("BiocManager")
> BiocManager::install("Rcpi")
> library(rcdk)
> library(Rcpi)
> Aspirin = load.molecules('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/2244/record/SDF/?record_type=3d&response_type=display') 
> Acetaminophen = load.molecules('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/1983/record/SDF/?record_type=3d&response_type=display')
> class(Aspirin)
> class(Acetaminophen)
> length(Aspirin)
> length(Acetaminophen)
> class(Aspirin[[1]])
> class(Acetaminophen[[1]])
> names(Aspirin[[1]])
> Aspirin_ecfp <- extractDrugExtended(Aspirin[[1]])
> Acetaminophen_ecfp <- extractDrugExtended(Acetaminophen[[1]])
> Aspirin_ecfp <- extractDrugExtendedComplete(Aspirin[[1]])
> Acetaminophen_ecfp <- extractDrugExtendedComplete(Acetaminophen[[1]])
> Aspirin_ecfp <- extractDrugExtended(Aspirin[[1]])
> Acetaminophen_ecfp <- extractDrugExtended(Acetaminophen[[1]])
> calcDrugFPSim(Aspirin_ecfp, Acetaminophen_ecfp, metric = 'tanimoto')
> Aspirin_ecfp <- extractDrugExtendedComplete(Aspirin[[1]])
> Acetaminophen_ecfp <- extractDrugExtendedComplete(Acetaminophen[[1]])
> calcDrugFPSim(Aspirin_ecfp, Acetaminophen_ecfp, metric = 'tanimoto', fptype = 'complete')
> Aspirin_ecfp <- extractDrugExtended(Aspirin[[1]], size = 100)
> Acetaminophen_ecfp <- extractDrugExtended(Acetaminophen[[1]], size = 100)
> calcDrugFPSim(Aspirin_ecfp, Acetaminophen_ecfp, metric = 'tanimoto')
> Aspirin_ecfp <- extractDrugExtendedComplete(Aspirin[[1]], size = 100)
> Acetaminophen_ecfp <- extractDrugExtendedComplete(Acetaminophen[[1]], size = 100)
> calcDrugFPSim(Aspirin_ecfp, Acetaminophen_ecfp, metric = 'tanimoto', fptype = 'complete')
> Aspirin_ecfp <- extractDrugExtended(Aspirin[[1]], size = 2000)
> Acetaminophen_ecfp <- extractDrugExtended(Acetaminophen[[1]], size = 2000)
> calcDrugFPSim(Aspirin_ecfp, Acetaminophen_ecfp, metric = 'tanimoto', fptype = 'complete')
