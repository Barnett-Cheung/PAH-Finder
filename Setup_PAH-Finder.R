## Install required packages if needed
install.packages(c("stringr", "dplyr", "caret", "randomForest", "ROCR", "e1071", "Rdisop"))

library(stringr)
library(dplyr)
library(caret)
library(randomForest)
library(ROCR)
library(e1071)  
library(Rdisop)

# Function to extract MSMS spectrum and create a data frame
MS_spectrum_extraction <- function(spectrum_str) {
  # Split the spectrum string into parts (m/z-intensity pairs)
  parts <- unlist(str_split(spectrum_str, "; "))
  
  # Initialize vectors to store m/z and intensity values
  mz_values <- numeric(length(parts))
  intensity_values <- numeric(length(parts))
  
  # Parsing each part to extract m/z and intensity
  for (i in seq_along(parts)) {
    part <- trimws(parts[i])
    split_part <- str_split(part, ",", simplify = TRUE)
    
    mz_values[i] <- as.numeric(split_part[1])
    intensity_values[i] <- as.numeric(split_part[2])
  }
  
  # Create a data frame with m/z, intensity, and group
  spectrum_df <- data.frame(
    mz = mz_values,
    intensity = intensity_values,
    group = rep("fragment", length(parts)), # group is set as 'fragment'
    stringsAsFactors = FALSE
  )
  
  return(spectrum_df)
}

## For MS-DIAL 
# Function to extract and normalize MSMS spectrum, creating a data frame
MS_spectrum_extraction_MSDIAL <- function(spectrum_str) {
  # Split the spectrum string into parts (m/z-intensity pairs)
  parts <- unlist(str_split(spectrum_str, " "))
  
  # Initialize vectors to store m/z and intensity values
  mz_values <- numeric(length(parts))
  intensity_values <- numeric(length(parts))
  
  # Parsing each part to extract m/z and intensity
  for (i in seq_along(parts)) {
    part <- trimws(parts[i])
    split_part <- str_split(part, ":", simplify = TRUE)
    
    mz_values[i] <- as.numeric(split_part[1])
    intensity_values[i] <- as.numeric(split_part[2])
  }
  
  # Normalize the intensities by Base Peak
  max_intensity <- max(intensity_values)
  normalized_intensities <- intensity_values / max_intensity * 999
  
  # Create a data frame with m/z, normalized intensity, and group
  spectrum_df <- data.frame(
    mz = mz_values,
    intensity = normalized_intensities,
    group = rep("fragment", length(parts)), # group is set as 'fragment'
    stringsAsFactors = FALSE
  )
  
  return(spectrum_df)
}

# Function to estimate molecular weight from a spectrum dataframe, considering top five m/z values
Estimate_molecular_ion <- function(spectrum_df) {
  # Ensure the spectrum data frame contains the necessary columns
  if (!("mz" %in% names(spectrum_df)) || !("intensity" %in% names(spectrum_df))) {
    stop("The spectrum data frame must contain 'mz' and 'intensity' columns.")
  }
  
  # Filter out fragments with intensity below a threshold, e.g., 20
  spectrum_filtered <- spectrum_df %>%
    filter(intensity > 20)
  
  # Select the top five fragments based on their m/z values
  top_fragments <- spectrum_filtered %>%
    arrange(desc(mz)) %>%
    slice_head(n = 5)
  
  # From those, choose the one with the largest intensity
  potential_molecular_ion <- top_fragments %>%
    arrange(desc(intensity)) %>%
    slice(1) %>%
    pull(mz)
  
  if (length(potential_molecular_ion) > 0) {
    return(potential_molecular_ion)
  } else {
    return(NA) # Return NA if no potential molecular ion is found
  }
}

# Helper function for calculating entropy
Entropy <- function(x) {
  # Normalize x so that it sums to 1 (probabilities)
  p <- x / sum(x)
  # Calculate entropy
  entropy <- -sum(p * log(p))
  return(entropy)
}

## Include precursor
Calculate_region_intensity_ratio <- function(spectra_list, molecular_weights, bin_step) {
  # Calculate the number of bins based on the bin step
  bins <- seq(0, 105, by = bin_step)
  processed_data <- list()
  ## processed_inchikeys <- character()  # For storing INCHIKEYs of processed spectra
  
  for (i in seq_along(spectra_list)) {
    spectrum <- spectra_list[[i]]
    spectrum <- spectrum[spectrum$group == "fragment", ]
    spectrum$mz <- as.numeric(spectrum$mz)
    mw <- as.numeric(molecular_weights[i])
    
    if (is.na(mw) || mw == 0) {
      next
    }
    
    # Calculate the total ion intensity
    total_ion_intensity <- sum(spectrum$intensity)
    
    spectrum$percent_of_mw <- (spectrum$mz / mw) * 100
    bin_intensities <- rep(0, length(bins) - 1)
    
    for (j in seq_along(bins)[-length(bins)]) {
      lower_bound <- bins[j]
      upper_bound <- bins[j + 1]
      bin_sum <- sum(spectrum$intensity[spectrum$percent_of_mw >= lower_bound & spectrum$percent_of_mw < upper_bound])
      bin_intensities[j] <- bin_sum / total_ion_intensity  # Normalize by total ion intensity
    }
    
    processed_data[[i]] <- bin_intensities
  }
  
  result_df <- do.call(rbind, processed_data)
  colnames(result_df) <- paste(bins[-length(bins)], "-", bins[-1] - bin_step, "%", sep="")
  ## result_df$INCHIKEY <- processed_inchikeys
  
  return(result_df)
}

Calculate_spectral_features <- function(spectra_list, molecular_weights) {
  features <- data.frame(

    Skewness = numeric(length(spectra_list)),
    Spectrum_Entropy = numeric(length(spectra_list)),
    BasePeak_Ratio = numeric(length(spectra_list))
    
  )
  
  error_indices <- c()  # To store indices of spectra with errors
  
  for (i in seq_along(spectra_list)) {
    tryCatch({
      spectrum <- spectra_list[[i]]
      mw <- molecular_weights[i]
      
      # Exclude precursor cluster (assuming ±6 Dalton range)
      # Check if all intensities are below 50 after excluding the precursor range
      if (sum(spectrum$intensity > 50, na.rm = TRUE) <= 1){
        
        features[i, ] <- 0
        
      } else {
        
        # Filter out rows where intensity is NA
        spectrum <- spectrum[!is.na(spectrum$intensity), ]
        spectrum <-  spectrum[spectrum$intensity > 20, ]
        
        # Skewness
        features$Skewness[i] <- skewness(spectrum$mz)
        
        # Spectrum Entropy
        features$Spectrum_Entropy[i] <- Entropy(spectrum$intensity)
        
        # Max Intensity to Max m/z Ratio
        max_intensity_mz <- spectrum$mz[which.max(spectrum$intensity)]
        max_mz <- max(spectrum$mz)
        features$BasePeak_Ratio[i] <- max_intensity_mz / max_mz
        
      }
    }, error = function(e) {
      # Capture the index of the spectrum where the error occurred
      error_indices <- c(error_indices, i)
      message(sprintf("Error in processing spectrum at index %d: %s", i, e$message))
    })
  }
  
  if (length(error_indices) > 0) {
    warning(sprintf("Errors occurred in processing spectra at indices: %s", paste(error_indices, collapse=", ")))
  }
  
  return(features )
}

### Updated function to identify doubly charged fragments 
Identify_doubly_charged_fragments <- function(spectrum_df) {
  # Ensure the dataframe is ordered by mz and filter out low-intensity peaks
  spectrum_df <- spectrum_df[order(spectrum_df$mz), ]
  spectrum_df <- spectrum_df[spectrum_df$intensity > 10,]
  
  # Calculate the differences between adjacent mz values
  mz_diffs <- diff(spectrum_df$mz)
  
  # Tolerance for measurement precision
  tolerance <- 0.02
  # Identify potential doubly charged fragments
  double_charge_indices <- which(mz_diffs > (0.5 - tolerance) & mz_diffs < (0.5 + tolerance))
  
  # Initially assume no doubly charged fragments are found
  has_doubly_charged <- FALSE
  
  # Loop through indices to check for corresponding singly charged fragments
  for (index in double_charge_indices) {
    # Calculate expected mz for singly charged fragment
    expected_singly_mz <- spectrum_df$mz[index] * 2
    
    # Check the entire spectrum for a match with expected_singly_mz allowing for 1 Dalton difference
    # due to isotopic variation and tolerance
    if (any(abs(spectrum_df$mz - expected_singly_mz) < 1 + tolerance)) {
      # If found, confirm the presence of doubly charged fragments
      has_doubly_charged <- TRUE
      break # Stop checking once a doubly charged fragment is confirmed
    }
  }
  
  # Return TRUE if doubly charged fragments are confirmed, otherwise FALSE
  return(has_doubly_charged)
}

# Function to use isotopic distribution to predict formula
Predict_formula <- function(spectrum_df, ppm = 20, minElements = "C8H6", maxElements = "C28H20") {
  # Sort the dataframe by m/z
  spectrum_df <- spectrum_df[order(spectrum_df$mz), ]
  
  # Assuming the main peak is the one with the highest intensity
  main_peak <- spectrum_df[which.max(spectrum_df$intensity), ]
  
  # Look for an isotopic peak (0.5 to 2 m/z difference with decent intensity)
  isotopic_peaks <- subset(spectrum_df, 
                           mz > main_peak$mz & mz <= main_peak$mz + 2 & 
                             intensity > main_peak$intensity * 0.05 & intensity < main_peak$intensity * 1.0 )
  
  # Initialize variables for results
  formula_result <- NULL
  dbe_result <- -Inf  # Initialize DBE to very low value
  
  # If there seems to be a reasonable isotopic distribution
  if (nrow(isotopic_peaks) > 0) {
    masses <- c(main_peak$mz, isotopic_peaks$mz[1])
    intensities <- c(main_peak$intensity, isotopic_peaks$intensity[1])
    
    # Use decomposeIsotopes to predict formula
    isotopes_result <- decomposeIsotopes(masses, intensities, ppm = ppm, 
                                         minElements = minElements, maxElements = maxElements)
    
    # Assuming a way to calculate or retrieve DBE from isotopes_result
    # This is a placeholder; you'll need to adjust based on actual output of decomposeIsotopes
    if (length(isotopes_result$formula) > 0) {
      formula_result <- isotopes_result$formula[1]
      dbe_result <- isotopes_result$DBE[1]  # Example to extract DBE
    }
  } 
  
  # If isotopic distribution isn't reasonable or no high DBE result from isotopes
  if (is.null(formula_result) || dbe_result <= 0) {  # Adjust the DBE condition as needed
    # Use decomposeMass for the main peak's mass
    mass_result <- decomposeMass(main_peak$mz, ppm = ppm, 
                                 minElements = minElements, maxElements = maxElements)
    
    # Choose the formula with the highest DBE from decomposeMass results
    # Placeholder for DBE extraction and comparison; adjust according to actual function's output
    if (length(mass_result$formula) > 0) {
      formula_result <- mass_result$formula[1]  # Example, assuming first result is most significant
    }
  }
  
  return(formula_result)
}
