
#...............................................................................
#                                                                              .
#  Statistical analysis for the paper “"Etiologies of rapidly progressive      .
#  dementias: A systematic review and meta-analysis of causes in worldwide     .
#  and Latin America”                                                          .
#                                                                              .
#  Author: Jonathan Cubas-Guillen                                              .
#                                                                              .
#  Date: 5th of August 2024                                                    .
#                                                                              .
#...............................................................................


# Load necessary libraries
library(meta)      # For meta-analysis
library(readxl)    # For reading Excel files
library(tidyverse) # For data manipulation
library(rstatix)   # For descriptive statistics and tests
library(ggplot2)   # For plotting
library(metafor)   # For meta-analysis functions
library(dmetar)    # For additional meta-analysis tools
library(gt)        # For creating tables
library(grid)      # For additional plotting

# Read data from an Excel file
BASE <- read_excel("base.xlsx")

# Filter out a specific author from the dataset
BASE2 <- BASE %>%
  filter(Author != "Grau-Rivera, 2015*")

## NOS (Newcastle-Ottawa Scale) Quality Assessment

# Perform quality assessment using the NOS tool
NOS <- rob(
  item1 = BASE2$Selection,        # Selection bias
  item2 = BASE2$Comparability,     # Comparability of study groups
  item3 = BASE2$Exposure,          # Exposure assessment
  cat1 = "**",                     # Category for high quality
  symb1 = "2",                     # Symbol for high quality
  col1 = "white",                  # Color for high quality
  cat2 = "",                       # Category for moderate quality
  symb2 = "+",                     # Symbol for moderate quality
  cat3 = "*",                      # Category for low quality
  symb3 = "1",                     # Symbol for low quality
  col3 = "white",                  # Color for low quality
  overall = BASE2$Total,           # Overall quality score
  cat.overall = c("***", "**"),    # Categories for overall quality
  symb.overall = c("3", "2"),      # Symbols for overall quality
  col.overall = c("white", "white")# Colors for overall quality
)

# Rename columns in the NOS object for clarity
NOS <- rename(NOS, S = A, C = B, E = C)

### NEURODEGENERATIVE DISEASES ANALYSIS

# Perform meta-analysis for proportions of neurodegenerative diseases
meta.nd <- metaprop(
  event = nND,        # Number of events (cases)
  n = nTotal,         # Total number of subjects
  studlab = Author,   # Study labels
  data = BASE2        # Data frame
)

# Create a forest plot for neurodegenerative diseases
forest(meta.nd,
       leftcols = c("studlab", "event", "n"),  # Columns to display on the left
       leftlabs = c("Author", "Cases", NA),    # Labels for the columns
       rob = NOS,                             # Risk of bias object
       rob.text = "NOS",                      # Text to display for risk of bias
       rob.xpos = 1.45,                       # X position for risk of bias text
       rob.legend = FALSE,                    # Do not display legend for risk of bias
       prediction = TRUE                     # Show prediction intervals
)
grid.text("Neurodegenerative diseases", x = 0.5, y = 0.8, gp = gpar(fontsize = 20))

### PUBLICATION BIAS ANALYSIS (NEURODEGENERATIVE DISEASES)

# Create a funnel plot to assess publication bias
funnel(meta.nd, studlab = TRUE, xlim = c(-7, 3))
title("Publication Bias Neurodegenerative")

# Perform publication bias test using Egger's test
metabias(meta.nd, method.bias = "Egger", k.min = 10)
eggers.test(meta.nd)

# Create a funnel plot with contours to detect publication bias
col.contour <- c("gray75", "gray85", "gray95")
meta::funnel(meta.nd,
             studlab = TRUE,
             xlim = c(-6, 3),
             contour = c(0.9, 0.95, 0.99),
             col.contour = col.contour
)
title("Contour-enhanced funnel plots - Neurodegenerative")

# Add a legend for contour plot
legend(
  x = 1, y = 0.01,
  legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
  fill = col.contour
)

## HETEROGENEITY ANALYSIS

# Update meta-analysis to include prediction intervals
prediction.nd <- update(meta.nd, prediction = TRUE)
summary(prediction.nd)

# Identify potential outliers
identify_outliers(BASE2, nND)

# Conduct influence analysis
meta.nd.inf <- metainf(meta.nd, random = TRUE)
summary(meta.nd.inf)

# Find outliers in meta-analysis
find.outliers(meta.nd)

## PLOTS WITHOUT OUTLIERS FOR NEURODEGENERATIVE DISEASES

# Filter out specific authors considered outliers
BASE.nd.outliers <- BASE2 %>%
  filter(Author != "Chandra, 2017") %>%
  filter(Author != "Day, 2018")

# Perform meta-analysis without outliers
meta.nd.outliers <- metaprop(
  event = nND,
  n = nTotal,
  studlab = Author,
  data = BASE.nd.outliers
)

# Create a forest plot without outliers
forest(meta.nd.outliers,
       leftcols = c("studlab", "event", "n"),
       leftlabs = c("Author", "Cases", NA),
       prediction = TRUE
)
grid.text("Neurodegenerative diseases", x = 0.5, y = 0.95, gp = gpar(fontsize = 20))

### INFLUENCE ANALYSIS

# Perform influence analysis and display results
meta.nd.inf <- InfluenceAnalysis(meta.nd, random = TRUE)
meta.nd.inf

# Convert influence analysis results to a dataframe and display it
tabla.nd <- as.data.frame(meta.nd.inf)
gt(tabla.nd)

# Create influence plots
plot(meta.nd.inf, "influence")

## LEAVE-ONE-OUT ANALYSIS FOR NEURODEGENERATIVE DISEASES

# Plot the influence of individual studies on effect size (ES)
plot(meta.nd.inf, "es")

# Plot the influence of individual studies on heterogeneity (I2)
plot(meta.nd.inf, "i2")

### BAUJAT PLOT

# Create a Baujat plot to assess the influence of each study on the meta-analysis
baujat(meta.nd, cex = 1, bg = "blue")
title("Baujat - Neurodegenerative")

# Fit a random-effects model to the meta-analysis data
nd.rma <- rma(
  yi = meta.nd$TE,              # Effect sizes
  sei = meta.nd$seTE,           # Standard errors of effect sizes
  method = meta.nd$method.tau,  # Method for estimating between-study variance
  test = "knha"                 # Test for heterogeneity
)

# Perform Global Outlier Sensitivity Analysis using GOSH
gosh.nd <- gosh(nd.rma)

# Plot the results of the GOSH analysis
plot(gosh.nd)

# Perform diagnostics on the GOSH results
gosh.nd.diag <- gosh.diagnostics(gosh.nd,
                                 km.params = list(centers = 2),      # Parameters for k-means clustering
                                 db.params = list(
                                   eps = 0.08,                      # Epsilon for DBSCAN algorithm
                                   MinPts = 50                      # Minimum number of points for DBSCAN clustering
                                 )
)

### LATIN AMERICA ANALYSIS FOR NEURODEGENERATIVE DISEASES

# Perform meta-analysis for neurodegenerative diseases with a focus on Latin American studies
meta.nd.latam <- metaprop(
  event = nND,        # Number of neurodegenerative disease cases
  n = nTotal,         # Total number of subjects
  studlab = Author,   # Study labels
  subgroup = LatinAmerica, # Subgroup for Latin America
  data = BASE2        # Data frame
)

# Create a forest plot for Latin American studies on neurodegenerative diseases
forest(meta.nd.latam,
       leftcols = c("studlab", "event", "n"),  # Columns to display on the left
       leftlabs = c("Author", "Cases", NA),    # Labels for the columns
       prediction = TRUE                      # Show prediction intervals
)
grid.text("Neurodegenerative diseases", x = 0.5, y = 0.95, gp = gpar(fontsize = 20))

## DEFINITION ANALYSIS FOR NEURODEGENERATIVE DISEASES

# Perform meta-analysis for neurodegenerative diseases based on different definitions
meta.nd.definition <- metaprop(
  event = nND,        # Number of neurodegenerative disease cases
  n = nTotal,         # Total number of subjects
  studlab = Author,   # Study labels
  subgroup = Definition, # Subgroup for definitions
  data = BASE2        # Data frame
)

# Create a forest plot for studies based on different definitions of neurodegenerative diseases
forest(meta.nd.definition,
       leftcols = c("studlab", "event", "n"),  # Columns to display on the left
       leftlabs = c("Author", "Cases", NA),    # Labels for the columns
       prediction = TRUE                      # Show prediction intervals
)
grid.text("Neurodegenerative diseases", x = 0.5, y = 0.95, gp = gpar(fontsize = 20))


### CJD ANALYSIS

# Perform a meta-analysis for Creutzfeldt-Jakob Disease (CJD) using proportion data
meta.cjd <- metaprop(
  event = nCJD,        # Number of CJD cases
  n = nTotal,         # Total number of subjects
  studlab = Author,   # Study labels
  data = BASE2        # Data frame containing the study data
)

# Create a forest plot for CJD data
forest(meta.cjd,
  leftcols = c("studlab", "event", "n"),  # Columns to display on the left
  leftlabs = c("Author", "Cases", NA),    # Labels for the columns
  rob = NOS,          # Risk of Bias (NOS) scores to include in the plot
  rob.text = "NOS",   # Label for the risk of bias
  rob.xpos = 1.45,    # X-position of the risk of bias text
  rob.legend = FALSE, # Whether to display the risk of bias legend
  prediction = TRUE  # Show prediction intervals
)
grid.text("Prion diseases", x = 0.5, y = 0.8, gp = gpar(fontsize = 20))

### CJD PUBLICATION BIAS

# Create a funnel plot to assess publication bias for CJD
funnel(meta.cjd, studlab = TRUE, xlim = c(-7, 3))
title("Publication Bias CJD")

# Perform Egger's test for publication bias
metabias(meta.cjd, method.bias = "Egger", k.min = 9)
eggers.test(meta.cjd)

# Create a contour-enhanced funnel plot to visualize publication bias more clearly
col.contour <- c("gray75", "gray85", "gray95")
meta::funnel(meta.cjd,
  xlim = c(-6, 3),
  contour = c(0.9, 0.95, 0.99),
  studlab = TRUE,
  col.contour = col.contour
)
title("Contour-enhanced funnel plots - Prion diseases")

# Add a legend to the contour-enhanced funnel plot
legend(
  x = 1, y = 0.01,
  legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
  fill = col.contour
)

# Re-run Egger's test for publication bias
metabias(meta.cjd, method.bias = "Egger", k.min = 9)
regtest(meta.cjd)

## HETEROGENEITY AND CJD

# Update the meta-analysis object to include prediction intervals
prediction.cjd <- update(meta.cjd, prediction = TRUE)
summary(prediction.cjd)

# Identify outliers in the CJD data
identify_outliers(BASE2, nCJD)

# Perform influence analysis to assess the impact of each study on the meta-analysis results
meta.cjd.inf <- metainf(meta.cjd, random = TRUE)
summary(meta.cjd.inf)

# Find outliers in the meta-analysis results
find.outliers(meta.cjd)

# Create a Baujat plot to assess the influence of each study on the meta-analysis
baujat(meta.cjd, cex = 1, bg = "blue")
title("Baujat - CJD")

## PLOTS WITHOUT OUTLIERS FOR CJD

# Exclude specific outlier studies from the data
BASE.cjd.outliers <- BASE2 %>%
  filter(Author != "Grau-Rivera, 2015")

# Perform meta-analysis without the outliers
meta.cjd.outliers <- metaprop(
  event = nCJD,        # Number of CJD cases
  n = nTotal,         # Total number of subjects
  studlab = Author,   # Study labels
  data = BASE.cjd.outliers  # Data frame with outliers excluded
)

# Create a forest plot for CJD data without outliers
forest(meta.cjd.outliers,
  leftcols = c("studlab", "event", "n"),  # Columns to display on the left
  leftlabs = c("Author", "Cases", NA),    # Labels for the columns
  prediction = TRUE  # Show prediction intervals
)
grid.text("Prion diseases", x = 0.5, y = 0.95, gp = gpar(fontsize = 20))

### INFLUENCE ANALYSIS

# Perform influence analysis to assess the impact of each study on the meta-analysis results
meta.cjd.inf <- InfluenceAnalysis(meta.cjd, random = TRUE)

# Convert influence analysis results to a data frame for further examination
meta.cjd.inf
tabla.cjd <- as.data.frame(meta.cjd.inf)

# Display the influence analysis results using a gt table
gt(tabla.cjd)

# Plot the influence of each study on the meta-analysis results
plot(meta.cjd.inf, "influence")

## LEAVE-ONE-OUT ANALYSIS FOR CJD

# Plot the effect size (ES) and heterogeneity (I2) for the leave-one-out analysis
plot(meta.cjd.inf, "es")
plot(meta.cjd.inf, "i2")


### CJD ANALYSIS

# Perform a meta-analysis for Creutzfeldt-Jakob Disease (CJD) using proportion data
meta.cjd <- metaprop(
  event = nCJD,        # Number of CJD cases
  n = nTotal,         # Total number of subjects
  studlab = Author,   # Study labels
  data = BASE2        # Data frame containing the study data
)

# Create a forest plot for CJD data
forest(meta.cjd,
       leftcols = c("studlab", "event", "n"),  # Columns to display on the left
       leftlabs = c("Author", "Cases", NA),    # Labels for the columns
       rob = NOS,          # Risk of Bias (NOS) scores to include in the plot
       rob.text = "NOS",   # Label for the risk of bias
       rob.xpos = 1.45,    # X-position of the risk of bias text
       rob.legend = FALSE, # Whether to display the risk of bias legend
       prediction = TRUE  # Show prediction intervals
)
grid.text("Prion diseases", x = 0.5, y = 0.8, gp = gpar(fontsize = 20))

### CJD PUBLICATION BIAS

# Create a funnel plot to assess publication bias for CJD
funnel(meta.cjd, studlab = TRUE, xlim = c(-7, 3))
title("Publication Bias CJD")

# Perform Egger's test for publication bias
metabias(meta.cjd, method.bias = "Egger", k.min = 9)
eggers.test(meta.cjd)

# Create a contour-enhanced funnel plot to visualize publication bias more clearly
col.contour <- c("gray75", "gray85", "gray95")
meta::funnel(meta.cjd,
             xlim = c(-6, 3),
             contour = c(0.9, 0.95, 0.99),
             studlab = TRUE,
             col.contour = col.contour
)
title("Contour-enhanced funnel plots - Prion diseases")

# Add a legend to the contour-enhanced funnel plot
legend(
  x = 1, y = 0.01,
  legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
  fill = col.contour
)

# Re-run Egger's test for publication bias
metabias(meta.cjd, method.bias = "Egger", k.min = 9)
regtest(meta.cjd)

## HETEROGENEITY AND CJD

# Update the meta-analysis object to include prediction intervals
prediction.cjd <- update(meta.cjd, prediction = TRUE)
summary(prediction.cjd)

# Identify outliers in the CJD data
identify_outliers(BASE2, nCJD)

# Perform influence analysis to assess the impact of each study on the meta-analysis results
meta.cjd.inf <- metainf(meta.cjd, random = TRUE)
summary(meta.cjd.inf)

# Find outliers in the meta-analysis results
find.outliers(meta.cjd)

# Create a Baujat plot to assess the influence of each study on the meta-analysis
baujat(meta.cjd, cex = 1, bg = "blue")
title("Baujat - CJD")

## PLOTS WITHOUT OUTLIERS FOR CJD

# Exclude specific outlier studies from the data
BASE.cjd.outliers <- BASE2 %>%
  filter(Author != "Grau-Rivera, 2015")

# Perform meta-analysis without the outliers
meta.cjd.outliers <- metaprop(
  event = nCJD,        # Number of CJD cases
  n = nTotal,         # Total number of subjects
  studlab = Author,   # Study labels
  data = BASE.cjd.outliers  # Data frame with outliers excluded
)

# Create a forest plot for CJD data without outliers
forest(meta.cjd.outliers,
       leftcols = c("studlab", "event", "n"),  # Columns to display on the left
       leftlabs = c("Author", "Cases", NA),    # Labels for the columns
       prediction = TRUE  # Show prediction intervals
)
grid.text("Prion diseases", x = 0.5, y = 0.95, gp = gpar(fontsize = 20))

### INFLUENCE ANALYSIS

# Perform influence analysis to assess the impact of each study on the meta-analysis results
meta.cjd.inf <- InfluenceAnalysis(meta.cjd, random = TRUE)

# Convert influence analysis results to a data frame for further examination
meta.cjd.inf
tabla.cjd <- as.data.frame(meta.cjd.inf)

# Display the influence analysis results using a gt table
gt(tabla.cjd)

# Plot the influence of each study on the meta-analysis results
plot(meta.cjd.inf, "influence")

## LEAVE-ONE-OUT ANALYSIS FOR CJD

# Plot the effect size (ES) and heterogeneity (I2) for the leave-one-out analysis
plot(meta.cjd.inf, "es")
plot(meta.cjd.inf, "i2")


### CJD AND LATAM ANALYSIS

# Perform a meta-analysis for Creutzfeldt-Jakob Disease (CJD) by Latin American subgroup
meta.cjd.latam <- metaprop(
  event = nCJD,        # Number of CJD cases
  n = nTotal,         # Total number of subjects
  studlab = Author,   # Study labels
  subgroup = LatinAmerica,  # Subgroup analysis for Latin America
  data = BASE2        # Data frame containing the study data
)

# Create a forest plot for CJD data by Latin American subgroup
forest(meta.cjd.latam,
       leftcols = c("studlab", "event", "n"),  # Columns to display on the left
       leftlabs = c("Author", "Cases", NA),    # Labels for the columns
       prediction = TRUE  # Show prediction intervals
)
grid.text("Prion diseases", x = 0.5, y = 0.95, gp = gpar(fontsize = 20))


### CJD BY DEFINITION

# Perform a meta-analysis for Creutzfeldt-Jakob Disease (CJD) by definition subgroup
meta.cjd.definition <- metaprop(
  event = nCJD,        # Number of CJD cases
  n = nTotal,         # Total number of subjects
  studlab = Author,   # Study labels
  subgroup = Definition,  # Subgroup analysis by definition
  data = BASE2        # Data frame containing the study data
)

# Create a forest plot for CJD data by definition subgroup
forest(meta.cjd.definition,
       leftcols = c("studlab", "event", "n"),  # Columns to display on the left
       leftlabs = c("Author", "Cases", NA),    # Labels for the columns
       prediction = TRUE  # Show prediction intervals
)
grid.text("Prion diseases", x = 0.5, y = 0.95, gp = gpar(fontsize = 20))

### AUTOIMMUNE ENCEPHALITIS ANALYSIS

# Perform a meta-analysis for Autoimmune Encephalitis (AI)
meta.ai <- metaprop(
  event = nAI,         # Number of AI cases
  n = nTotal,          # Total number of subjects
  studlab = Author,    # Study labels
  data = BASE2         # Data frame containing the study data
)

# Create a forest plot for AI data
forest(meta.ai,
       leftcols = c("studlab", "event", "n"),  # Columns to display on the left
       leftlabs = c("Author", "Cases", NA),    # Labels for the columns
       rob = NOS,               # Risk of bias data
       rob.text = "NOS",        # Label for risk of bias
       rob.xpos = 1.45,         # Position of the risk of bias label
       rob.legend = FALSE,      # Do not show the legend for risk of bias
       prediction = TRUE        # Show prediction intervals
)
grid.text("Autoimmune Encephalitis", x = 0.5, y = 0.8, gp = gpar(fontsize = 20))


#### HETEROGENEITY AND AI

# Update the meta-analysis to include prediction intervals
prediction.ai <- update(meta.ai, prediction = TRUE)
summary(prediction.ai)

# Identify outliers in the dataset based on the number of AI cases
identify_outliers(BASE2, nAI)

# Perform influence analysis to identify influential studies
meta.ai.inf <- metainf(meta.ai, random = TRUE)
summary(meta.ai.inf)

# Identify outliers in the meta-analysis
find.outliers(meta.ai)

## LEAVE-ONE-OUT ANALYSIS AND AI

# Conduct a leave-one-out sensitivity analysis
meta.ai.inf <- InfluenceAnalysis(meta.ai, random = TRUE)

# Plot the leave-one-out analysis results (effect sizes and I2 statistics)
plot(meta.ai.inf, "es")
plot(meta.ai.inf, "i2")


## PLOTS WITHOUT OUTLIERS AI

# Filter out the outliers from the dataset
BASE.ai.outliers <- BASE2 %>%
  filter(Author != "Chandra, 2017")

# Perform a meta-analysis on the filtered dataset
meta.ai.outliers <- metaprop(
  event = nAI,         # Number of AI cases
  n = nTotal,          # Total number of subjects
  studlab = Author,    # Study labels
  data = BASE.ai.outliers  # Filtered data frame
)

# Create a forest plot for AI data without outliers
forest(meta.ai.outliers,
       leftcols = c("studlab", "event", "n"),  # Columns to display on the left
       leftlabs = c("Author", "Cases", NA),    # Labels for the columns
       prediction = TRUE  # Show prediction intervals
)
grid.text("Autoimmune Encephalitis", x = 0.5, y = 0.95, gp = gpar(fontsize = 20))


# Summarize the influence analysis
summary(meta.ai.inf)


### PUBLICATION BIAS AND AI

# Create a funnel plot to assess publication bias for AI studies
funnel(meta.ai, studlab = TRUE, xlim = c(-7, 3))
title("Publication Bias AI")

# Perform Egger's test to statistically assess publication bias
metabias(meta.ai, method.bias = "Egger", k.min = 9)
eggers.test(meta.ai)

# Enhance the funnel plot with contour lines representing significance levels
col.contour <- c("gray75", "gray85", "gray95")
meta::funnel(meta.ai,
             xlim = c(-6, 3),
             contour = c(0.9, 0.95, 0.99),
             studlab = TRUE,
             col.contour = col.contour
)
title("Contour-enhanced funnel plots - Autoimmune Encephalitis")

# Add a legend to the funnel plot to explain contour levels
legend(
  x = 1, y = 0.01,
  legend = c  ("p < 0.1", "p < 0.05", "p < 0.01"),
  fill = col.contour
)


### AI AND LATIN AMERICA

# Perform a subgroup meta-analysis for studies from Latin America
meta.ai.latam <- metaprop(
  event = nAI,          # Number of AI cases
  n = nTotal,           # Total number of subjects
  studlab = Author,     # Study labels
  subgroup = LatinAmerica,  # Subgroup for Latin America
  data = BASE2          # Data frame containing the study data
)

# Create a forest plot for the Latin American subgroup
forest(meta.ai.latam,
       leftcols = c("studlab", "event", "n"),  # Columns to display on the left
       leftlabs = c("Author", "Cases", NA),    # Labels for the columns
       prediction = TRUE  # Show prediction intervals
)
grid.text("Autoimmune Encephalitis", x = 0.5, y = 0.95, gp = gpar(fontsize = 20))


### AI AND DEFINITION

# Perform a subgroup meta-analysis based on different definitions used in the studies
meta.ai.definition <- metaprop(
  event = nAI,          # Number of AI cases
  n = nTotal,           # Total number of subjects
  studlab = Author,     # Study labels
  subgroup = Definition,  # Subgroup for different definitions
  data = BASE2          # Data frame containing the study data
)

# Create a forest plot for the Definition subgroup
forest(meta.ai.definition,
       leftcols = c("studlab", "event", "n"),  # Columns to display on the left
       leftlabs = c("Author", "Cases", NA),    # Labels for the columns
       prediction = TRUE  # Show prediction intervals
)
grid.text("Autoimmune Encephalitis", x = 0.5, y = 0.95, gp = gpar(fontsize = 20))
