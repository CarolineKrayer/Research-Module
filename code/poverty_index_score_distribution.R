#####################################################################################
######### This file generates the distribution plots in Figure 5.1 and E.1 ##########
#####################################################################################

# Create subset of first SISBEN dataset and
# plot distribution of Poverty Index Score
# between 1994 and 2003.

library(foreign)
library(dplyr)
library(ggplot2)

# Housekeeping
data_path <- "../data/"
out_path <- "../Figures/"

# Create subset of first SISBEN dataset
# (only observations of Poverty Index Score by year).
sisben <- read.dta(paste0(data_path, "sisben_aejep.dta"))
sisben_small <- with(sisben, data.frame(score=puntaje, date=fencuesta))
year <- format(as.Date(sisben_small$date, format="%d/%m/%Y"),"%Y")
sisben_small <- mutate(sisben_small, year)
sisben_small <- select(sisben_small, -matches("date"))
sisben_small <- filter(sisben_small, year >= 1994 & year <= 2003)
rm(sisben, year)

# Alternatively, load pre-processed data, "sisben_small.csv".
#sisben_small <- read.table(paste0(data_path, "sisben_small.csv"), header = TRUE)
#sisben_small <- data[, -3]
#names(sisben_small)[1] <- "score"
#names(sisben_small)[2] <- "year"

# Plot of Poverty Index Score distribution by year.
for(i in (1994:2003)){
   data <- filter(sisben_small, sisben_small$year == i)
   histplot <- ggplot(data=data, aes(data$score)) +
               geom_histogram(binwidth = 0.5, aes(y=(..count..)/sum(..count..)),
                              fill="blue") +
               scale_y_continuous(labels=scales::percent, limits=c(0, 0.06)) +
               labs(title=i) +
               labs(x="Poverty index score", y="Percent")

   final <- histplot +
            geom_vline(xintercept=47, color="red", size=1) +
            theme_minimal() +
            theme(plot.margin=margin(0.4, 0.4, 0.4, 0.4, "cm"),
                  text = element_text(size=18))

   ggsave(final, file=paste0(out_path, "plot_", i, ".png"), height=4, width=8)
}
