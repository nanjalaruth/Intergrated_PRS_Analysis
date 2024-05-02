#!/usr/local/bin/Rscript

#Install packages
if(!require(pacman)) install.packages("pacman")

pacman::p_load(
  tidyverse,
  data.table)

#Association analysis of Protsig with DM,BP markers
trait <- fread("${lipid_trait}")
prs_pcs <- fread("${ug_prs_pcs}")
traits <- trait %>% 
  dplyr::filter(!(is.na(sangerid)|sangerid == "")) %>% 
  rename(FID = sangerid) %>% 
  setnames(., names(.), c(gsub("-","_",colnames(.))))
traits_cov <- merge(prs_pcs, traits, by="FID")

#Prediction - Lipid traits
predictor <- "s1"
outcome_variables <- c("TC", "HbA1c", "LDL_C" , "HDL_C", "systolic", "TG")
models <- list()
for (outcome_variable in outcome_variables) {
  formula <- as.formula(paste0(outcome_variable, "~", predictor, "+ sex + age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
  model <- lm(formula, data=traits_cov)
  models[[outcome_variable]] <- model
  cat("Linear model for", outcome_variable, ":\\n")
  print(summary(model))
  conf_intervals <- confint(model)
  cat("Confidence Intervals:\\n")
  print(conf_intervals)
  cat("\\n")
}

# Generate a forest plot
results <- data.frame(
  Outcome = character(),
  Estimate = numeric(),
  LowerCI = numeric(),
  UpperCI = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each model to extract coefficients for 's1' and store in the dataframe
for (outcome_variable in outcome_variables) {
  model <- models[[outcome_variable]]
  
  # Check if the 's1' term exists in the model summary
  term <- "s1"  # This is the main effect of interest
  if (term %in% rownames(summary(model)\$coefficients)) {
    est <- coef(summary(model))[term, "Estimate"]
    ci <- confint(model, term, level = 0.95)
    
    # Append results to dataframe
    results <- rbind(results, data.frame(
      Outcome = outcome_variable,
      Estimate = est,
      LowerCI = ci[1],
      UpperCI = ci[2]
    ))
  } else {
    # Handle the case where the main effect is not present
    results <- rbind(results, data.frame(
      Outcome = outcome_variable,
      Estimate = NA,
      LowerCI = NA,
      UpperCI = NA
    ))
  }
}

# Adjust the plot with modified 'fatten' parameter
vitc_forest_plot <- ggplot(results, aes(x=reorder(Outcome,Estimate), 
                                        y=Estimate, ymin=LowerCI, ymax=UpperCI)) + 
  geom_linerange(size=5,position=position_dodge(width = 0.5),color="blue4") + 
  geom_hline(yintercept=0, lty=2) + 
  geom_point(size=3, shape=21, colour="white", stroke = 0.5,
             position=position_dodge(width = 0.5)) +
  coord_flip() + 
  theme_classic() +
  theme(axis.title.y = element_blank()) + 
  theme( panel.background = element_blank(), panel.grid = element_blank()) +
  theme(text=element_text(size = 14)) +ylab("Beta(95%CI)")


pdf("${blood_trait}_forest_plot.pdf", width = 6, height = 4)
print(vitc_forest_plot)
dev.off()

