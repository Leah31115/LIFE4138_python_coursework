from pathlib import Path
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.cluster import hierarchy
import sys

"""
This script will only accept two deseq .tsv files. Required data columns include: "gene_id", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"
All generated files and plots will be saved automatically into a newly made folder titled "gene_analysis_output" within your working directory.
The number of upregulated and downregulated genes for both treatments will be calculated.
A statistical summary of p-values and log fold changes for all genes are saved as separate tsv files for both conditions titled:
- condition1_statistical_summary.tsv
- condition2_statistical_summary.tsv

The top significantly upregulated and downregulated genes for both conditions is saved as a file titled top_differentially_expressed_genes.tsv

The plots generated include: volcano plots, MA plots and histograms of p-values. Heatmap and heirarchal cluster is generated for the top differentially expressed genes.
"""
# Obtain the file path of the TSV files
Tk().withdraw() # Prevent using tkinter's full GUI since it's not needed
# Open A_vs_D_file
print("Please select a DESeq DGEA file (condition 1)")
AvD_file_path_tk = askopenfilename(filetypes=[("tab separated values", ".tsv")]) # Open a dialog box to obtain the file A vs D Deseq TSV filepath as a string. Only allows the user to select .tsv files to prevent errors
AvD_file_path = Path(AvD_file_path_tk) # Handle file path using pathlib. Used later to make a directory

# Handle error if no file is selected
try:
    f = open(AvD_file_path_tk)
except FileNotFoundError:
    print("No file selected"),
    sys.exit()

# Open G_vs_D file
print("Please select another DESeq DGEA file you wish to compare against (condition 2)")
GvD_file_path_tk = askopenfilename(filetypes=[("tab separated values", ".tsv")]) # Open a dialog box to obtain the file G vs D Deseq TSV filepath as a string. Only allows the user to select .tsv files to prevent errors

# Handle error if no file is selected
try:
    f = open(GvD_file_path_tk)
except FileNotFoundError:
    print("No file selected"),
    sys.exit()

# Make a new directory within the parent folder containing the TSV files to store the output files from this script
new_dir = AvD_file_path.parent / "gene_analysis_output" # Make new directory in parent folder
new_dir.mkdir(exist_ok=True) # If the directory already exists, then another will not be made

# Load in TSV files
AvD_data = pd.read_csv(AvD_file_path_tk, sep='\t').dropna() # Remove any NA values
GvD_data = pd.read_csv(GvD_file_path_tk, sep='\t').dropna() # Remove any NA values

# Make a negative log10 adjusted p-value column (required later for volcano plots)
AvD_data["neg_log10_padj"] = -np.log(AvD_data["padj"])
GvD_data["neg_log10_padj"] = -np.log(GvD_data["padj"])
# Make log_baseMean column (required later for MA plots)
AvD_data["log_baseMean"] = -np.log(AvD_data["baseMean"])
GvD_data["log_baseMean"] = -np.log(GvD_data["baseMean"])


# Number of significantly upregulated and downregulated genes (based on p-value and fold change thresholds) for each comparison.
# Assume significance occurs when adjusted p-value < 0.05 and log2foldchange >= 2.00
# For A vs D Deseq
# Upregulated genes
AvD_up_reg = AvD_data[(AvD_data["log2FoldChange"] >= 2) & (AvD_data["padj"] < 0.05)]
count_AvD_upreg = len(AvD_up_reg)
print(f"\nThere are {count_AvD_upreg} upregulated genes for condition 1.")

# Downregulated genes
AvD_down_reg = AvD_data[(AvD_data["log2FoldChange"] <= -2) & (AvD_data["padj"] < 0.05)]
count_AvD_downreg = len(AvD_down_reg)
print(f"There are {count_AvD_downreg} downregulated genes for condition 1.")

# For G vs D Deseq
# Upregulated genes
GvD_up_reg = GvD_data[(GvD_data["log2FoldChange"] >= 2) & (GvD_data["padj"] < 0.05)]
count_GvD_upreg = len(GvD_up_reg)
print(f"There are {count_GvD_upreg} upregulated genes for condition 2.")

# Downregulated genes
GvD_down_reg = GvD_data[(GvD_data["log2FoldChange"] <= -2) & (GvD_data["padj"] < 0.05)]
count_GvD_downreg = len(GvD_down_reg)
print(f"There are {count_GvD_downreg} downregulated genes for condition 2.")


# Summary of p-values and log fold changes across all genes for each condition.
# For A vs D Deseq
subset_AvD_data = AvD_data[["log2FoldChange", "pvalue", "padj"]]
summary_AvD_data = subset_AvD_data.describe()
round_summary_AvD_data = summary_AvD_data.round({"log2FoldChange": 3}) # Round log2FoldChange to 3dp
print("\nHere is the summary of p-values and log-fold changes for all genes for condition 1:")
print(round_summary_AvD_data)
# Save statistical summary
round_summary_AvD_data.to_csv(new_dir / "condition1_statistical_summary.tsv", index=False, header=True, sep='\t')

# For G vs D Deseq
subset_GvD_data = GvD_data[["log2FoldChange", "pvalue", "padj"]]
summary_GvD_data = subset_GvD_data.describe()
round_summary_GvD_data = summary_GvD_data.round({"log2FoldChange": 3}) # Round log2FoldChange to 3dp
print("\nHere is the summary of p-values and log-fold changes for all genes for condition 2:")
print(round_summary_GvD_data)
# Save statistical summary
round_summary_GvD_data.to_csv(new_dir / "condition1_statistical_summary.tsv", index=False, header=True, sep='\t')

# Set threshold for the top significant genes
# Attempt at making a dynamic threshold which is flexible with other datasets
top_padj_threshold = (round_summary_GvD_data["log2FoldChange"].iloc[6] * round_summary_GvD_data["log2FoldChange"].iloc[2] * 10)**2 # 75% * mean * 10 and squared to make number positive

# Generate volcano plots to visualize the significance and magnitude of changes in gene expression for each condition.
# For A vs D Deseq
# Volcano plot which highlights up-regulated and down-regulated genes
plt.scatter(x=AvD_data["log2FoldChange"], y=AvD_data["neg_log10_padj"], color="grey", label="Not significant", s=1.5)
plt.scatter(x=AvD_up_reg["log2FoldChange"], y=AvD_up_reg["neg_log10_padj"], color="orange", label="Up-regulated")
plt.scatter(x=AvD_down_reg["log2FoldChange"], y=AvD_down_reg["neg_log10_padj"], color="blue", label="Down-regulated")
plt.title("Condition 1 Volcano Plot", fontsize=15)
plt.xlabel("Log2 Fold Change")
plt.ylabel("-Log(adjusted pvalue)")
plt.axhline(top_padj_threshold, color="black", linestyle="--") # top significantly expressed genes
plt.axvline(2, color="black", linestyle="--")
plt.axvline(-2, color="black", linestyle="--")
plt.legend()

# Label only the extreme gene points for clarity
# Annotate gene names at the extreme horizontal ends along the x-axis
# Extreme upregulated genes
AvD_up_reg_lab = AvD_data[(AvD_data["log2FoldChange"] >= 5) & (AvD_data["padj"] < 0.01)]

# Extreme downregulated genes
AvD_down_reg_lab = AvD_data[(AvD_data["log2FoldChange"] <= -5) & (AvD_data["padj"] < 0.01)]

# Annotate gene names at the extreme horizontal ends along the x-axis
for i, gene in enumerate(AvD_up_reg_lab["gene_id"]):
    plt.annotate(gene, (AvD_up_reg_lab["log2FoldChange"].iloc[i], AvD_up_reg_lab["neg_log10_padj"].iloc[i]), size=7)

for i, gene in enumerate(AvD_down_reg_lab["gene_id"]):
    plt.annotate(gene, (AvD_down_reg_lab["log2FoldChange"].iloc[i], AvD_down_reg_lab["neg_log10_padj"].iloc[i]), size=7)


# Annotate gene names at the extreme top, along the y-axis, with large neg_log10_padj values
# Filter significant up-regulated and down-regulated genes with large negative log10 adjusted p values
AvD_reg_lab_top = AvD_data[AvD_data["neg_log10_padj"] >= top_padj_threshold]

# Annotate gene names along the y-axis
for i, gene in enumerate(AvD_reg_lab_top["gene_id"]):
    plt.annotate(gene, (AvD_reg_lab_top["log2FoldChange"].iloc[i], AvD_reg_lab_top["neg_log10_padj"].iloc[i]), size=7)

plt.savefig(new_dir / "condition1_volcano_plot.png") # save fig must be before show else the images will be empty
plt.show()


# Generate volcano plot for G vs D Deseq
# Highlight up-regulated and down-regulated genes
plt.scatter(x=GvD_data["log2FoldChange"], y=GvD_data["neg_log10_padj"], color="grey", label="Not significant", s=1.5)
plt.scatter(x=GvD_up_reg["log2FoldChange"], y=GvD_up_reg["neg_log10_padj"], color="orange", label="Up-regulated")
plt.scatter(x=GvD_down_reg["log2FoldChange"], y=GvD_down_reg["neg_log10_padj"], color="blue", label="Down-regulated")
plt.title("Condition 2 Volcano Plot", fontsize=15)
plt.xlabel("Log2 Fold Change")
plt.ylabel("-Log(adjusted pvalue)")
plt.axhline(top_padj_threshold, color="black", linestyle="--") # top significantly expressed genes
plt.axvline(2, color="black", linestyle="--")
plt.axvline(-2, color="black", linestyle="--")
plt.legend()

# Label only the extreme gene points for clarity
# Annotate gene names at the extreme horizontal ends along the x-axis
# Extreme upregulated genes
GvD_up_reg_lab = GvD_data[(GvD_data["log2FoldChange"] >= 5) & (GvD_data["padj"] < 0.01)]

# Extreme downregulated genes
GvD_down_reg_lab = GvD_data[(GvD_data["log2FoldChange"] <= -5) & (GvD_data["padj"] < 0.01)]

# Annotate gene names at the extreme horizontal ends along the x-axis
for i, gene in enumerate(GvD_up_reg_lab["gene_id"]):
    plt.annotate(gene, (GvD_up_reg_lab["log2FoldChange"].iloc[i], GvD_up_reg_lab["neg_log10_padj"].iloc[i]), size=7)

for i, gene in enumerate(GvD_down_reg_lab["gene_id"]):
    plt.annotate(gene, (GvD_down_reg_lab["log2FoldChange"].iloc[i], GvD_down_reg_lab["neg_log10_padj"].iloc[i]), size=7)

# Annotate gene names at the extreme top, along the y-axis, with large neg_log10_padj values
# Filter significant up-regulated and down-regulated genes with large negative log10 p-adjusted values
GvD_reg_lab_top = GvD_data[GvD_data["neg_log10_padj"] >= top_padj_threshold]

# Annotate gene names along the y-axis
for i, gene in enumerate(GvD_reg_lab_top["gene_id"]):
    plt.annotate(gene, (GvD_reg_lab_top["log2FoldChange"].iloc[i], GvD_reg_lab_top["neg_log10_padj"].iloc[i]), size=7)

plt.savefig(new_dir / "condition2_volcano_plot.png")
plt.show()


# MA plot to display the relationship between the log fold change and mean expression.
# For A vs D Deseq
plt.scatter(x=AvD_data["log_baseMean"],y=AvD_data['log2FoldChange'], color="grey", label="Non-significant")
plt.scatter(x=AvD_up_reg["log_baseMean"], y=AvD_up_reg["log2FoldChange"], color="orange", label="Up-regulated")
plt.scatter(x=AvD_down_reg["log_baseMean"], y=AvD_down_reg["log2FoldChange"], color="blue", label="Down-regulated")
plt.xlabel("Log Mean expression")
plt.ylabel("Log2 Fold Change")
plt.title("MA plot of the relationship between deseq gene log mean expression and log fold change for condition 1", fontsize=13)
plt.legend()
plt.savefig(new_dir / "condition1_MAplot.png")
plt.show()

# MA plot for G vs D Deseq
plt.scatter(x=GvD_data["log_baseMean"],y=GvD_data['log2FoldChange'], color="grey", label="Non-significant")
plt.scatter(x=GvD_up_reg["log_baseMean"], y=GvD_up_reg["log2FoldChange"], color="orange", label="Up-regulated")
plt.scatter(x=GvD_down_reg["log_baseMean"], y=GvD_down_reg["log2FoldChange"], color="blue", label="Down-regulated")
plt.xlabel("Log Mean expression")
plt.ylabel("Log2 Fold Change")
plt.title("MA plot of the relationship between deseq gene log mean expression and log fold change for condition 2", fontsize=13)
plt.legend()
plt.savefig(new_dir / "condition2_MAplot.png")
plt.show()


# Histogram of p-values to assess the distribution of statistical significance.
# For A vs D Deseq
plt.hist(AvD_data['padj'], bins=20, color="blue", ec="black")
plt.xlabel("Adjusted p-values")
plt.ylabel("Frequency")
plt.title("Histogram of condition 1 adjusted p-value frequencies", fontsize=15)
plt.savefig(new_dir / "condition1_gene_pvalue_histogram.png")
plt.show()

# For G vs D Deseq
plt.hist(GvD_data['padj'], bins=20, color="blue", ec="black")
plt.xlabel("Adjusted p-values")
plt.ylabel("Frequency")
plt.title("Histogram of condition 2 adjusted p-value frequencies", fontsize=15)
plt.savefig(new_dir / "condition2_gene_pvalue_histogram.png")
plt.show()


# A table of significantly up-regulated and down-regulated genes with their corresponding fold changes, p-values, and adjusted p-values
# Merge the up-regulated and down-regulated genes for both conditions for the values with the highest statistical significance (found at the top of the volcano plot)
# Features at the edges are most promising candidates for further investigation but these are not included in this table
# Assign treatment group categories as new columns to the data frames for genes with high statistical significance
subset_AvD_reg_top = AvD_reg_lab_top[["gene_id", "log2FoldChange", "pvalue", "padj"]]
subset_AvD_reg_top = subset_AvD_reg_top.copy() # Made a copy to stop the warning which did not like editing (adding a column) after the dataframe had already been sliced
subset_AvD_reg_top["treatment_group"] = 'Condition_1' 

subset_GvD_reg_top = GvD_reg_lab_top[["gene_id", "log2FoldChange", "pvalue", "padj"]]
subset_GvD_reg_top = subset_GvD_reg_top.copy() # Made a copy to stop the warning which did not like editing (adding a column) after the dataframe had already been sliced
subset_GvD_reg_top["treatment_group"] = 'Condition_2' 

# Combine both significantly up-regulated and down-regulated genes into one dataframe for both conditions
data_frames = [subset_AvD_reg_top, subset_GvD_reg_top]
sig_genes = pd.concat(data_frames)
round_sig_genes = sig_genes.round({"log2FoldChange": 3}) # Round log2FoldChange to 3dp
print("\nHere are the significant up-regulated and down-regulated genes for both conditions:")
print(round_sig_genes)
# Save significant up-regulated and down-regulated genes for both conditions as a tsv file
round_sig_genes.to_csv(new_dir / "top_differentially_expressed_genes.tsv", index=False, header=True, sep='\t')


# Heatmap and Heirarchal clustering of the top differentially expressed genes, using log2FoldChange, to illustrate gene expression patterns and relationships across the conditions.
# All mean values from both files are the same so used log2FoldChange instead
# Each row contains significantly expressed genes for at least one condition. Significant genes taken from sig_genes dataframe
unique_genes = sig_genes["gene_id"].unique() # Remove duplicates present within both AvD and GvD (YNR075W for the case of my files)

# Make a dataframe of significantly expressed gene names with their log2FoldChange values for each condition
cluster_values = {"gene_id": unique_genes,
                "Condition1_log2FoldChange": [], "Condition2_log2FoldChange": []}

# Extract the gene information for both treatment groups
for name in cluster_values["gene_id"]:
    AvD_current_row = AvD_data[AvD_data["gene_id"] == name] # find row with the gene_id = name in AvD data
    AvD_log2FoldChange_value = AvD_current_row["log2FoldChange"].iloc[0] # obtain the value from the located row
    cluster_values["Condition1_log2FoldChange"].append(AvD_log2FoldChange_value) # append value to the list for the "AvD_log2FoldChange" key within the dictionary

    GvD_current_row = GvD_data[GvD_data["gene_id"] == name] # find row with the gene_id = name in GvD data
    GvD_log2FoldChange_value = GvD_current_row["log2FoldChange"].iloc[0] # obtain the value from the located row
    cluster_values["Condition2_log2FoldChange"].append(GvD_log2FoldChange_value) # append value to the list for the "GvD_log2FoldChange" key within the dictionary

# Convert dictionary into dataframe
cluster_values_df = pd.DataFrame.from_dict(cluster_values)
# Index the gene_id
cluster_values_df = cluster_values_df.set_index("gene_id")

# Plot heatmap
sns.heatmap(cluster_values_df, annot=True)
plt.ylabel("Gene")
plt.title("Heatmap of the top differentially expressed genes for both treatments")
plt.savefig(new_dir / "significant_gene_expression_heatmap.png")
plt.show()

# Heirarchal clustering
# Transpose the dataframe for gene clustering
cluster_values_df = cluster_values_df.T # Each significant gene will now be a column instead of a row

# Plot heirarchal cluster
dendrogram = hierarchy.dendrogram(hierarchy.linkage(cluster_values_df.iloc[:, 1:].T,method="average"), labels=cluster_values_df.columns[1:])
plt.title("Gene Expression Hierachal Clustering")
plt.xlabel("Genes")
plt.ylabel("Distance")
plt.tight_layout()
plt.xticks(rotation=90)
plt.savefig(new_dir / "significant_gene_expression_heirarchal_cluster.png")
plt.show()