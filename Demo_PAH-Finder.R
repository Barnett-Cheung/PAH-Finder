## Please feel free to contact me for any preprocessing and data analysis help using this workflow
## I am more than happy to assist in using PAH-Finder.
## Email: zixuanzhang245@gmail.com

## 2024-04
# This script is a demonstration of applying PAH-Finder for PAH candidate recognition in GC-MS analysis

# Created: 2024-04
# Edited : 2024-06  
# Created by: Zixuan Zhang
# Edited by : Zixuan Zhang
# Version: 0.0.1

## Set the working directiory
setwd("C:/Users/46355/Desktop/test_1")

demo <- read.csv("Demo_MSMS.csv")

# Process each MSMS spectrum and store in a list
List_MS_spectra <- lapply(demo$MSMS, MS_spectrum_extraction)

molecular_ion_mz <- as.numeric(unlist(lapply(List_MS_spectra, Estimate_molecular_ion)))

# Calculate each machine learning features respectively 
# Rename and combine into a dataframe for prediction
RIR_1 <- Calculate_region_intensity_ratio(List_MS_spectra, molecular_ion_mz, 52.5)[, 2]

RIR_2_3 <- Calculate_region_intensity_ratio(List_MS_spectra, molecular_ion_mz, 35)[, c(1,3)]

RIR_4 <- Calculate_region_intensity_ratio(List_MS_spectra, molecular_ion_mz, 10.5)[, 10]

Spectral_features <- Calculate_spectral_features(List_MS_spectra, molecular_ion_mz)

df_ML_features <- cbind(RIR_1, RIR_2_3, RIR_4, Spectral_features)

colnames(df_ML_features) <- c("RIR1","RIR2","RIR3","RIR4",
                              "Skewness", "Entropy", "BasePeakRatio")

RFmodel <- readRDS(file="RFmodel_PAH-Finder.rds")

predictions_prob <- predict(RFmodel, newdata = df_ML_features, type = "prob")
# If you need the predictions to be binary classes (0 or 1) rather than probabilities,
# you can apply a threshold, here we use 0.5 as threshold to distinguish PAH from background
binary_predictions <- ifelse(predictions_prob > 0.5, 1, 0)

ML_prediction <- as.numeric(binary_predictions[ , 2])

DB_frag_check <- unlist(lapply(List_MS_spectra, Identify_doubly_charged_fragments))

Suspect_indices <- which(ML_prediction == 1 & DB_frag_check == T)

Suspect_list <- List_MS_spectra[Suspect_indices]

# Iterate through the Suspect_list and write each element to a CSV file
for (i in seq_along(Suspect_list)) {
  # Create a unique file name for each element
  file_name <- paste0("Suspect_Spectrum_", Suspect_indices[i], ".csv")
  
  # Write the element to a CSV file
  write.csv(Suspect_list[[i]], file = file_name, row.names = FALSE)
}

## In this demo, relative intensity is at the scale of 1000. 
## We set the cutoff as 20 (i.e., 1%) to view the spectra
Int_cutoff <- 10
i <- 1
Suspect_list[[i]][Suspect_list[[i]]$intensity > Int_cutoff, ]

## From the suspect list, I highly recommend to predict formula individually for each spectrum
## PAH formula is characterized by the large ratio of Carbon to Hydrogen that can be easily recognized
## As m/z increase, the space of possible explanation of formula grows significantly larger
## An automatic function for formula prediction is provided, manual validation is needed for robust annotation

## From the spectra, molecular ion peaks for PAHs are usually largest m/z with high intensity.
## Use it for formula prediction with Rdisop 
# Example usage
mass <-  278.1077                                                                                                          
decomposeMass(mass, ppm = 10, minElements = "C10H8", maxElements = "C28H20")

# Isotope pattern: masses and their corresponding intensities
masses <- c(278.1077, 279.1121)
intensities <- c(999, 156.3196)

# Using decomposeIsotopes to predict formula with isotope information
decomposeIsotopes(masses, intensities, ppm = 10, 
                  minElements = "10H8", maxElements = "C28H20")

## In above output, C22H14 is the most reasonable formula that explain the spectrum



##########################
## Optional to use.
# Assuming this is a valid spectrum dataframe
formula_prediction <- lapply(Suspect_list, Predict_formula)

print(formula_prediction)

# Replace NULL values with NA
formula_predictions_filled <- lapply(formula_predictions, function(x) if(is.null(x)) NA else x)

# Now you can unlist, keeping NA for the original NULL values, preserving their indices
unlisted_formula_predictions <- unlist(formula_predictions_filled, use.names = FALSE)

