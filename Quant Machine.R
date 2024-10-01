library(readxl)
library(openxlsx)
library(janitor)
library(tidyverse)

#Part I - Metabolite Group Comparison

# Input the total extracted volume and the amount injected from that volume.
extraction_volume <- 100
injection_volume <- 5
# Input below or import a list of weights and group classifications for each sample in order of the
# Sample ID. Also, input the name of the control group for calculating the fold change.
my_weights <- c(10, 10, 10, 10, 10, 10)
my_groups <- c("A", "A", "B", "B", "C", "C")
unique_groups <- unique(my_groups)
control_group <- "A"

# Assigns file path of original workbook to variable, workbook, and creates a new workbook
# to write edited sheets to.
workbook <- "workbook filepath"
new_workbook <- createWorkbook("New Workbook")

metabolites_data <- c()
metabolite_names <- c()
anova_Fs <- c()

for (i in 1:length(excel_sheets(workbook))){
  # Opens specific sheet and removes extraneous rows.
  sheet_name <- excel_sheets(workbook)[i]
  sheet <- data.frame(read_excel(workbook, sheet = sheet_name))
  sheet <- tail(sheet, -2)
  sheet <- row_to_names(sheet, 1)
  sheet <- head(sheet, -6)
  # Filters out rows with no Sample ID.
  sheet <- sheet[!is.na(sheet$`Sample ID`), ]
  # Returns specific columns.
  sheet <- select(sheet, Filename, `Sample ID`, `Sample Name`, `Response Ratio`, `Calc Amt`)
  # Assigns the weight and group classification for each sample.
  sheet$weight <- my_weights
  sheet$group <- my_groups
  # Calculates new columns.
  sheet$Total_Amount <- (as.numeric(sheet$`Calc Amt`) * extraction_volume) / injection_volume
  sheet$Normalized_Amount <- sheet$Total_Amount / sheet$weight
  control_avg <- mean(sheet$Normalized_Amount[which(sheet$group==control_group)], na.rm = TRUE)
  sheet$FC <- sheet$Normalized_Amount / control_avg
  # Adds worksheet to new workbook and writes data to it.
  addWorksheet(new_workbook, sheet_name)
  writeDataTable(new_workbook, sheet_name, sheet)
  
  mean_fcs <- c()
  sd_fcs <- c()
  all_group_fcs <- c()
  p_values <- c()
  
  # Calculates mean fold change, standard deviation, and p_value between the control and other groups.
  for (i in unique_groups){
    group_fcs <- sheet$FC[which(sheet$group==i)]
    mean_fc <- mean(group_fcs, na.rm = TRUE)
    mean_fcs <- c(mean_fcs, mean_fc)
    sd_fc <- sd(group_fcs, na.rm = TRUE)
    sd_fcs <- c(sd_fcs, sd_fc)
    all_group_fcs <- rbind(all_group_fcs, data.frame(fcs = I(list(group_fcs))))
    t_test <- t.test(sheet$FC[which(sheet$group==control_group)], group_fcs)
    p_value <- t_test$p.value
    p_values <- c(p_values, p_value)
  }
  
  # Compiles all data for each analyte into a dataframe and adds it to a list
  metabolite_name <- sheet_name
  metabolite_names <- c(metabolite_names, metabolite_name)
  metabolite_data <- data.frame(group = unique_groups, mean_fc = mean_fcs, sd_fc = sd_fcs, all_fcs = all_group_fcs, p_value = p_values)
  metabolites_data <- c(metabolites_data, list(metabolite_data))
  
  # Perfomrms ANOVA and adds it to list
  fcs_table <- unnest(metabolite_data, fcs)[c("group", "fcs")]
  fcs_anova <- oneway.test(fcs ~ group, fcs_table)
  anova_F <- fcs_anova$statistic
  anova_Fs <- c(anova_Fs, anova_F)
}

# Assigns metabolite name to each dataframe
metabolites_data <- setNames(metabolites_data, metabolite_names)
anova_data <- setNames(anova_Fs, metabolite_names)

# Part II - Plotting Metabolites of Interest

# Choose which metabolite to plot
my_metabolite <- "NAD+"
# Exclude groups
group_exclusion <- c("")
# Choose order of groups. By default the order of groups will be based on which groups appear first within the sample list.
group_order <- unique_groups[unique_groups != group_exclusion]
# Color of each group
my_colors <- c("red", "orange", "yellow", "green", "blue", "violet")
# Choose the order of annotation significance bars. By default, the annotations for each comparison will go from comparing the leftmost group to the control
# starting at the bottom then working its way to the next group to the right and annotating upwards.
my_comparisons <- group_order[group_order != control_group]

# Formats the comparisons for geom_signif
list_comparisons <- list(c(control_group, my_comparisons[1]))
for (i in my_comparisons[-1]){
  list_comparisons <- append(list_comparisons, list(c(control_group, i)))
}

# Finds and formats data for ggplot
my_metabolite_index <- which(metabolite_names == my_metabolite)
my_metabolite_data <- metabolites_data[[my_metabolite_index]]
my_metabolite_data <- filter(my_metabolite_data, !group %in% group_exclusion)
my_metabolite_data$group <- factor(my_metabolite_data$group, levels = group_order)

# Creates annotations for each comparison
my_p_values <- my_metabolite_data$p_value[match(my_comparisons, my_metabolite_data$group)]
my_annotations <- case_when(my_p_values >= 0.05 ~ "ns",
                            my_p_values < 0.05 & my_p_values >= 0.01 ~ "*",
                            my_p_values < 0.01 & my_p_values >= 0.001 ~ "**",
                            my_p_values < 0.001 ~ "***")

# Unnests data for automatic plotting and geom_point
unnested_data <- unnest(my_metabolite_data, fcs)

# Finds points of reference for automatic plotting
y_label <- paste(my_metabolite, "fold change")
top_sig <- max(my_metabolite_data$mean_fc, na.rm = TRUE) + max(my_metabolite_data$sd_fc, na.rm = TRUE)
top_fc <- max(unnested_data$fcs, na.rm = TRUE)
top_y <- max(c(top_sig, top_fc))
y_annotations <- c()

for (i in 1:length(my_annotations)){
  y_annotation <- top_y + ((i-1)*(0.1*top_y))
  y_annotations <- c(y_annotations, y_annotation)
}

setwd("filepath for plot")

tiff("filename.tiff", res = 300, height = 5, width = 5, units = "in")
ggplot(data=my_metabolite_data, aes(x=group, y=mean_fc, fill = group)) +
  geom_col(width=0.6) +
  scale_fill_manual(values = setNames(my_colors, my_metabolite_data$group)) +
  geom_errorbar(aes(ymin = mean_fc - sd_fc, ymax = mean_fc + sd_fc), width=.1) +
  geom_point(data = unnested_data, aes(x=group,y=fcs, fill="black"), position = position_jitterdodge(jitter.width = 0.2)) +
  guides(fill="none") +
  geom_signif(comparisons = list_comparisons, map_signif_level = TRUE, y_position = y_annotations, annotations = my_annotations) +
  xlab(NULL) +
  ylab(y_label) +
  scale_y_continuous(limits = c(0, ceiling(max(y_annotations)+(top_y*0.1))), expand = c(0, 0), breaks = seq(0,ceiling(max(y_annotations)+(top_y*0.1)),1)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), 
        axis.line = element_line(color = "#3D4852"),
        axis.ticks = element_line(color = "#3D4852"))

dev.off()



# Saves new workbook with calculations.
new_workbook_file <- "fiepath.xlsx"
saveWorkbook(new_workbook, file=new_workbook_file, overwrite = TRUE)




