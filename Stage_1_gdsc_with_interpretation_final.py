#!/usr/bin/env python
# coding: utf-8

# <center>
#   <h1 style='color:green;'> Genomics of Drug Sensitivity in Cancer (GDSC)</h1>
#     <h2>Exploratory Data Analysis</h2>
#     <h3>HackBio AI for Genomics Internship</h3>
#     <h4>Team valine-asparagine-alanine</h4>
#     <h4>Uncovering drug response patterns, cancer-type sensitivity profiles, and genomic biomarkers of pharmacological outcome</h4>
# </center>
# 
# 
# **Dataset:** <a href = 'https://github.com/HackBio-Internship/public_datasets/raw/refs/heads/main/GDSC.xlsx'>GDSC Public Repository </a>
# 

# ## Objective
# 
# This notebook performs an exploratory analysis of the GDSC dataset to identify:
# 
# - Drug sensitivity patterns
# - Cancer-specific therapeutic responses
# - Highly selective anti-cancer drugs
# - Genomic determinants of drug response
# - Precision oncology insights
# 
# 
# 
# ## Scientific Background
# 
# The Genomics of Drug Sensitivity in Cancer (GDSC) project links:
# - Cancer cell lines
# - Genomic features
# - Drug response measurements
# 
# The aim is to understand how tumour molecular profiles influence therapeutic sensitivity.
# 
# 

# import the necessary libraries and the GDSC dataset

# # Setup & Data Loading

# In[10]:


import numpy as np
import pandas as pd


# In[11]:


gdsc = pd.read_excel("GDSC.xlsx")


# In[12]:


gdsc.head()


# The Key Variables in the Datasets include *LN_IC50, AUC, Z_SCORE* (sensitivity metrics); *DRUG_NAME, TARGET, TARGET_PATHWAY* (therapeutic agents); *CELL_LINE_NAME, TCGA_DESC, Microsatellite instability Status* (Cancer Models); and CNA, Gene Expression, Methylation (Omics)

# In[13]:


gdsc.shape
gdsc.columns
gdsc.info()
gdsc.describe()
gdsc["DRUG_NAME"].nunique()
gdsc["TCGA_DESC"].nunique()
gdsc["CELL_LINE_NAME"].nunique()


# In[14]:


missing = gdsc.isnull().sum().sort_values(ascending=False)
missing.head(20)


# There is no missing values

# In[15]:


gdsc = gdsc.drop_duplicates()
gdsc


# There are duplicates in the datasets

# In[99]:


plt.hist(gdsc['LN_IC50'], bins=50)
plt.xlabel('LN_IC50')
plt.ylabel('Frequency')
plt.title('Distribution of Drug Sensitivity')
plt.show()


# # Drug Sensitivity Patterns

# Here, we aim to find out which drugs are highly potent, which ones fail universally, and which ones have highly specific targets (high variability).

# ## Visualisation of distribution of drug sensitivity across samples

# In[72]:


# Calculate average sensitivity and variance for each drug
drug_stats = gdsc.groupby('DRUG_NAME').agg(
    Median_LN_IC50=('LN_IC50', 'median'),
    Median_AUC=('AUC', 'median'),
    Var_LN_IC50=('LN_IC50', 'var'),
    Cell_Line_Count=('CELL_LINE_NAME', 'count')
).reset_index()

# Top 5 most effective drugs (Lowest mean LN_IC50)
most_effective = drug_stats.sort_values(by='Median_LN_IC50').head(5)
print("Top 5 Most Effective Drugs (Lowest LN_IC50):")
print(most_effective[['DRUG_NAME', 'Median_LN_IC50', 'Median_AUC']])

# Top 5 least effective drugs (Highest mean LN_IC50)
least_effective = drug_stats.sort_values(by='Median_LN_IC50', ascending=False).head(5)
print("\nTop 5 Least Effective Drugs (Highest LN_IC50):")
print(least_effective[['DRUG_NAME', 'Median_LN_IC50', 'Median_AUC']])

# Top 5 most variable drugs (Context-dependent/Targeted treatments)
most_variable = drug_stats.sort_values(by='Var_LN_IC50', ascending=False).head(5)
print("\nTop 5 Most Variable Drugs (Potential Personalized Medicine Candidates):")
print(most_variable[['DRUG_NAME', 'Var_LN_IC50']])


# **INTERPRETATION**
# 
# **1. Which Drugs Appear to be the Most Effective?**
# 
# The top five most effective drugs are: Romidepsin, Bortezomib, Sepantronium bromide, Docetaxel, and Daporinad.
# 
# **2. Which Drugs Show the Least Effectiveness?**
# 
# The five least effective drugs are: Ascorbate (Vitamin C), N-acetyl cysteine, Glutathione, Alpha-lipoic acid, and Temozolomide.
# 
# Four out of five of these compounds *(ascorbate, N-acetyl cysteine, glutathione, and alpha-lipoic acid)* are antioxidants or dietary supplements. Because they lack localized cytotoxic machinery, they do not induce cell death on their own.
# 
# **3. Are There Drugs with Highly Variable Responses Across Cell Lines?**
# 
# The five most variable drugs are Gemcitabine, AZD5991, Daporinad, Docetaxel, and BI-2536.
# 
# The high response variance observed in drugs like ***Gemcitabine*** and ***AZD5991*** serves as the perfect statistical springboard for genomic biomarker tracking.

# ## Cancer types analysis
# Do certain cancer types respond better to specific drugs?

# ## Visualisation of drug sensitivity across cancer types

# In[77]:


heatmap_data = cancer_drug.pivot(
    index="TCGA_DESC",
    columns="DRUG_NAME",
    values="LN_IC50"
).sort_index(ascending=False)

plt.figure(figsize=(60,10))

sns.heatmap(
    heatmap_data,
    cmap="coolwarm",
    linewidths=0.5
)

plt.title("Drug Sensitivity Across Cancer Types")
plt.xlabel("Drug")
plt.ylabel("Cancer Type")

plt.show()


# **Interpretation**
# 
# The drug sensitivity across the cancer types was visualised using a heatmap. The drugs were plotted against the cancer types to see the pattern. Blue represents higher sensitivity(lower IC50), and red represents lower sensitivity (higher IC50). The white cells represent that the drug was not tested for that cancer. The heatmap is more clumsy due to the larger number of drugs. So, filtration of top sensitive drugs is required for a much cleaner plot. 

# In[88]:


heatmap_data = cancer_drug.pivot(
    index="TCGA_DESC",
    columns="DRUG_NAME",
    values="LN_IC50"
)

sorted_drugs = heatmap_data.median(axis=0).sort_values(ascending=True).index
heatmap_data_sorted = heatmap_data[sorted_drugs]

sorted_cancers = heatmap_data_sorted.median(axis=1).sort_values(ascending=True).index
heatmap_data_final = heatmap_data_sorted.loc[sorted_cancers]

plt.figure(figsize=(60, 12))

sns.heatmap(
    heatmap_data_final,
    cmap="coolwarm",       
    linewidths=0.2,        
    xticklabels=True,      
    yticklabels=True,      
    cbar_kws={'label': 'LN_IC50 (Lower/Blue = More Sensitive)'}
)

plt.title("Landscape of Genomic Drug Sensitivity Profiling Across Human Cancer Types Sorted by Pan-Cancer Efficacy Gradient")
plt.xlabel("Drug Name (Most Sensitive → Least Sensitive)")
plt.ylabel("Cancer Type (TCGA Code)")

plt.xticks(rotation=90, fontsize=9, ha='center')
plt.yticks(fontsize=12)

plt.tight_layout()
plt.show()


# **Romidepsin** and **sepantronium bromide** stand out due to their remarkably low median, having the absolute lowest median concentration threshold.  
# This indicates that they have high broad-spectrum cytotoxicity, completely eliminating cell viability across a massive portion of the treated concentration range.

# ## Filteration of top sensitive drugs and visualisation through heatmap

# In[101]:


top_drugs = (
    cancer_drug.groupby("DRUG_NAME")["LN_IC50"]
    .std()
    .sort_values(ascending=False)
    .head(30)
    .index
)
filtered_heatmap = heatmap_data[top_drugs]
plt.figure(figsize=(25,10))

sns.heatmap(
    filtered_heatmap,
    cmap="coolwarm"
)

plt.xticks(rotation=90)

plt.title("Top Selective Drugs Across Cancer Types")

plt.show()


# After filtering the top sensitive drugs using standard deviation values, a heatmap was plotted again. This heatmap shows drug sensitivity across the cancer types. The drugs Tozasertib, AZD5991, Gemcitabine, Dasatinib and Nilotinib show more variable sensitivity across the cancer types compared to other drugs. That is, they have lower IC50 values in certain cancer types than others. This shows that certain drugs are more sensitive in certain cancer types. 

# ## Identification of the most selective drugs across Cancer types

# In[89]:


drug_selectivity = (
    cancer_drug.groupby("DRUG_NAME")["LN_IC50"]
    .std()
    .sort_values(ascending=False)
)


# In[90]:


top_variable = drug_selectivity.head(15)

plt.figure(figsize=(10,6))

sns.barplot(
    x=top_variable.values,
    y=top_variable.index
)

plt.xlabel("Standard Deviation of LN_IC50")
plt.ylabel("Drug")
plt.title("Most Selective Drugs Across Cancer Types")

plt.show()


# The most selective drugs across the cancer types were identified and visualised using a bar plot. The plot shows the top 15 selective drugs( these drugs have more sensitivity in certain cancer types). This was identified using the standard deviation of IC50 values, which shows the variability of the drug response across the cancers. The tozasertib drug shows a higher standard deviation (more variability in response), meaning that the drug works well only in certain cancer types.

# ## Visualisation of cancer-specific response to the drug Tozasertib

# In[102]:


Tozasertib = gdsc[gdsc["DRUG_NAME"] == "Tozasertib"]
Tozasertib_summary = (
    Tozasertib.groupby("TCGA_DESC")["LN_IC50"]
    .mean()
    .sort_values()
)

Tozasertib_summary


# In[92]:


order = (
    Tozasertib.groupby("TCGA_DESC")["LN_IC50"]
    .median()
    .sort_values()
    .index)
plt.figure(figsize=(14,6))

sns.boxplot(
    data=Tozasertib,
    x="TCGA_DESC",
    y="LN_IC50",
    order=order
)

plt.xticks(rotation=90)

plt.title("Cancer-Specific Response to Tozasertib")

plt.show()


# Taking the most sensitive drug tozasertib as an example, its response across the cancer types was visualised using a box plot. We can conclude that the drug works well with low IC50 values in the DLBC and ALL cancer types compared to other types. 

# ## Identification of Most Drug-Sensitive Cancer Types

# In[93]:


cancer_type_sensitivity = (
    gdsc.groupby("TCGA_DESC")["LN_IC50"]
    .mean()
    .sort_values()
)
top_sensitive = cancer_type_sensitivity.head(20)

top_sensitive
plt.figure(figsize=(10,6))

sns.barplot(
    x=top_sensitive.values,
    y=top_sensitive.index
)

plt.xlabel("Mean LN_IC50")

plt.title("Most Drug-Sensitive Cancer types")

plt.show()


# ## Identification of Most Selective response Cancer Types

# In[94]:


cancer_type_sensitivity = (
    gdsc.groupby("TCGA_DESC")["LN_IC50"]
    .std()
    .sort_values(ascending=False)
)
top_sensitive = cancer_type_sensitivity.head(20)

top_sensitive
plt.figure(figsize=(10,6))

sns.barplot(
    x=top_sensitive.values,
    y=top_sensitive.index
)

plt.xlabel("Std LN_IC50")

plt.title("Most selective response Cancer types")

plt.show()


# Most drug-sensitive (broadly sensitive to drugs) and drug-selective (strong response to only certain drugs) cancer types were identified using the mean and standard deviation of IC50 values. The cancer type CLL is both sensitive and selective to drugs. LAML, DLBC, and ALL are also more drug-sensitive cancer types.

# ## Interpretation
# Yes, certain cancer types like CLL, LAML, DLBC and ALL respond better to certain drugs. Drugs like AZD5991, Tozasertib, and Gemcitabine are more selective across all the cancer types

# ## Cancer cell line analysis 
# Do certain cancer cell lines respond better to specific drugs?

# ## Identification of the best and the worst cell line response to drugs

# In[95]:


cellline_drug_response = (
    gdsc.groupby(["CELL_LINE_NAME", "DRUG_NAME"])["LN_IC50"]
    .mean()
    .reset_index()
)


# In[96]:


best_cellline_per_drug = (
    cellline_drug_response.loc[
        cellline_drug_response.groupby("DRUG_NAME")["LN_IC50"].idxmin()
    ]
)

best_cellline_per_drug.head(10)


# In[97]:


worst_cellline_per_drug = (
    cellline_drug_response.loc[
        cellline_drug_response.groupby("DRUG_NAME")["LN_IC50"].idxmax()
    ]
)

worst_cellline_per_drug.head(10)


# Analysing the drug response in the cell lines, the best and worst responding cell line for each drug were identified. 

# ## Vizualisation of Cell-line specific drug response

# In[98]:


heatmap_cellline = cellline_drug_response.pivot(
    index="CELL_LINE_NAME",
    columns="DRUG_NAME",
    values="LN_IC50"
)
plt.figure(figsize=(18,12))

sns.heatmap(
    heatmap_cellline,
    cmap="coolwarm"
)

plt.title("Cell Line-Specific Drug Response")

plt.show()


# The drug response against the cell lines was visualised using a heatmap. Blue represents higher sensitivity(lower IC50), and red represents lower sensitivity (higher IC50). The white cells represent that the drug was not tested for that cell line. The heatmap is more clumsy due to the larger number of drugs. So, filtration of the top sensitive drugs is required for a much cleaner plot. 

# ## Filteration of top sensitive drugs and visualisation

# In[33]:


top_drugs = (
    cellline_drug_response.groupby("DRUG_NAME")["LN_IC50"]
    .std()
    .sort_values(ascending=False)
    .head(25)
    .index
)

filtered_heatmap_cellline = heatmap_cellline[top_drugs]


# In[34]:


plt.figure(figsize=(18,12))

sns.heatmap(
    filtered_heatmap_cellline,
    cmap="coolwarm"
)

plt.title("Cell Line-Specific Drug Response")

plt.show()


# After filtering the top sensitive drugs using standard deviation values, a heatmap was plotted again with the top 25 sensitive drugs. This heatmap shows drug sensitivity across the cell lines. The drugs Gemcitabine, AZD5991, Daporinad, BI-2536 and Dasatinib show the most variable sensitivity across the cell lines compared to other drugs. That is, they have lower IC50 values in certain cell lines than others. This shows that certain drugs are more sensitive in certain cell lines. 

# ## Identification of most selective drugs across cell lines

# In[35]:


drug_variability = (
    cellline_drug_response.groupby("DRUG_NAME")["LN_IC50"]
    .std()
    .sort_values(ascending=False)
)


# In[36]:


top_variable_cellline = drug_variability.head(25)

plt.figure(figsize=(10,6))

sns.barplot(
    x=top_variable_cellline.values,
    y=top_variable_cellline.index
)

plt.xlabel("Standard Deviation of LN_IC50")
plt.ylabel("Drug")
plt.title("Most Selective Drugs Across Cell lines")

plt.show()


# The most selective drugs across the cell lines were identified and visualised using a bar plot. The plot shows the top 5 selective drugs( these drugs have more sensitivity in certain cell lines). This was identified using the standard deviation of IC50 values, which shows the variability of the drug response across the cell lines. The Gemcitabine drug shows a higher standard deviation (more variability in response), meaning that the drug works well only in certain cell lines.

# ## Visualisation of cell line-specific response to the drug Gemcitabine

# In[37]:


Gemcitabine = gdsc[gdsc["DRUG_NAME"] == "Gemcitabine"]
gem_summary = (
    Gemcitabine.groupby("CELL_LINE_NAME")["LN_IC50"]
    .mean()
    .sort_values()
)

top_sensitive = gem_summary.head(20)
top_resistant = gem_summary.tail(20)

selected_celllines = pd.concat([top_sensitive, top_resistant])
filtered_gem = Gemcitabine[
    Gemcitabine["CELL_LINE_NAME"].isin(selected_celllines.index)
]


# In[38]:


plt.figure(figsize=(14,6))

sns.boxplot(
    data=filtered_gem,
    x="CELL_LINE_NAME",
    y="LN_IC50",
    order=selected_celllines.index
)

plt.xticks(rotation=90)

plt.title("Most Sensitive and Resistant Cell Lines to Gemcitabine")

plt.show()


# Taking the most sensitive drug gemcitabine as an example, its response across the cell lines was visualised using a box plot. Since there are more than 700 cell lines, the top sensitive and the resistant cell lines were filtered and plotted against IC50 values of gemcitabine. We can conclude that the drug works well with low IC50 values in the CTB-1 and MC-IXC cancer types compared to other types. 

# ## Visualisation of top drug-sensitive cell lines 

# In[39]:


cellline_sensitivity = (
    gdsc.groupby("CELL_LINE_NAME")["LN_IC50"]
    .mean()
    .sort_values()
)
top_sensitive = cellline_sensitivity.head(20)

top_sensitive


# In[40]:


plt.figure(figsize=(10,6))

sns.barplot(
    x=top_sensitive.values,
    y=top_sensitive.index
)

plt.xlabel("Mean LN_IC50")

plt.title("Most Drug-Sensitive Cell Lines")

plt.show()


# ## Visualisation of most selective cell lines

# In[41]:


cellline_sensitivity = (
    gdsc.groupby("CELL_LINE_NAME")["LN_IC50"]
    .std()
    .sort_values(ascending=False)
)
top_sensitive = cellline_sensitivity.head(20)

top_sensitive
plt.figure(figsize=(10,6))

sns.barplot(
    x=top_sensitive.values,
    y=top_sensitive.index
)

plt.xlabel("Std LN_IC50")

plt.title("Most selective response Cell Lines")

plt.show()


# Most drug-sensitive (broadly sensitive to drugs) and drug-selective (strong response to only certain drugs) cell lines were identified using the mean and standard deviation of IC50 values. The cell lines MOLM-13 and MDA-MB-175-VII were the most sensitive and selective cell lines to drugs. DOHH-2 and SU-DHL-5 are also more drug-sensitive cell lines.

# ## Interpretation
# Yes, certain cell lines respond better to specific drugs. The cell lines MOLM-13 and MDA-MB-175-VII were the most sensitive and selective cell lines to drugs. DOHH-2 and SU-DHL-5 are also more drug-sensitive cell lines. Drugs like Gemcitabine, AZD5991 and Daporinad were more sensitive across all cancer cell lines.

# ## Cluster plot showing the clustering among the top sensitive and resistant cell lines with drug response
# What patterns exist in the drug response across cell lines?
# 

# In[103]:


cellline_sensitivity = (
    gdsc.groupby("CELL_LINE_NAME")["LN_IC50"]
    .mean()
    .sort_values()
)
top_sensitive = cellline_sensitivity.head(25)

top_resistant = cellline_sensitivity.tail(25)

selected_celllines = pd.concat([
    top_sensitive,
    top_resistant
]).index
cellline_matrix = gdsc.pivot_table(
    index="CELL_LINE_NAME",
    columns="DRUG_NAME",
    values="LN_IC50",
    aggfunc="mean"
)
top_drugs = (
    gdsc.groupby("DRUG_NAME")["LN_IC50"]
    .std()
    .sort_values(ascending=False)
    .head(25)
    .index
)
filtered_matrix = cellline_matrix.loc[
    selected_celllines,
    top_drugs
]
filtered_matrix = filtered_matrix.fillna(
    filtered_matrix.mean()
)
sns.clustermap(
    filtered_matrix,
    cmap="coolwarm",
    figsize=(14,10)
)


# A cluster heatmap was plotted for drug response across cell lines to see the pattern. The top sensitive and resistant cell lines were only used to get a much cleaner plot. The plot clearly shows some grouping of red and blue representing the resistant and sensitive cell lines, respectively. Clustering within the cell lines shows that they respond similarly. Likewise, the clustering among the drugs shows they work similarly. 

# # Genomic, transcriptomic and epigenomic influence on drug response
# Are mutations (CNA), gene expression or methylation level influencing drug sensitivity at all?

# In[43]:


gdsc[["CNA", "Gene Expression", "Methylation"]].head()


# In[44]:


gdsc[["CNA", "Gene Expression", "Methylation"]].dtypes


# ## Box Plot showing the effect of CNA on drug response

# In[45]:


plt.figure(figsize=(8,6))

sns.boxplot(
    data=gdsc,
    x="CNA",
    y="LN_IC50"
)

plt.title("Effect of CNA on Drug Sensitivity")

plt.show()


# In[46]:


gdsc.groupby("CNA")["LN_IC50"].describe()


# ## Statistical test to identify the significant difference between the two groups

# In[47]:


from scipy.stats import ttest_ind


# In[48]:


cna_yes = gdsc[gdsc["CNA"] == "Y"]["LN_IC50"]

cna_no = gdsc[gdsc["CNA"] == "N"]["LN_IC50"]
ttest_ind(cna_yes, cna_no)


# ## Interpretation
# 
# Statistical test gives a p-value = 0.005497744925679454 (<0.05), explaining that there is a statistical difference between the CNA yes group and the CNA no group. The box plot clearly shows that there is a clear difference between the IC50 values of the two groups. Since the IC50 values of the yes group (2.82) are much lower than the no group(3.12), the Copy number (CNA) influences the drug response. A higher copy number can give a better drug response.

# ## Box Plot showing the effect of Gene Expression on drug response

# In[49]:


plt.figure(figsize=(8,6))

sns.boxplot(
    data=gdsc,
    x="Gene Expression",
    y="LN_IC50"
)

plt.title("Gene Expression vs Drug Sensitivity")

plt.show()
gdsc.groupby("Gene Expression")["LN_IC50"].describe()


# ## Statistical test to identify the significant difference between the two groups

# In[50]:


expr_high = gdsc[gdsc["Gene Expression"] == "Y"]["LN_IC50"]

expr_low = gdsc[gdsc["Gene Expression"] == "N"]["LN_IC50"]

ttest_ind(expr_high, expr_low)


# ## Interpretation
# 
# Statistical test gives a p-value = 2.1731496680513867e-22 (<0.05), explaining that there is a statistical difference between the gene expression yes group and the no group. The box plot clearly shows that there is a clear difference between the IC50 values of the two groups. Since the IC50 values of the yes group (2.8) are much lower than the no group(2.38), the gene expression influences the drug response. 

# ## Box Plot showing the effect of Methylation on drug response

# In[51]:


plt.figure(figsize=(8,6))

sns.boxplot(
    data=gdsc,
    x="Methylation",
    y="LN_IC50"
)

plt.title("Methylation vs Drug Sensitivity")

plt.show()
gdsc.groupby("Methylation")["LN_IC50"].describe()


# ## Statistical test to identify the significant difference between the two groups

# In[52]:


Methyl_yes = gdsc[gdsc["Methylation"] == "Y"]["LN_IC50"]

Methyl_no = gdsc[gdsc["Methylation"] == "N"]["LN_IC50"]

ttest_ind(Methyl_yes, Methyl_no)


# ## Interpretation 
# 
# Statistical test gives a p-value = 0.5858384515515012 (>0.05), explaining that there is no statistical difference between the methylation yes group and the no group. The box plot clearly shows that there is no clear difference between the IC50 values of the two groups. Since the IC50 values of the yes group (2.81) are the same as the no group(2.79), the Metylation does not influence the drug response. 

# ## Analysis of Drug response across the pathways

# In[53]:


gdsc['TARGET_PATHWAY'].value_counts()


# In[54]:


pathway_response = gdsc.groupby(
    'TARGET_PATHWAY'
)['LN_IC50'].mean().sort_values()

pathway_response


# In[55]:


plt.figure(figsize=(10,6))

sns.barplot(
    x=pathway_response.values,
    y=pathway_response.index
)

plt.xlabel("Mean LN_IC50")
plt.ylabel("Target Pathway")
plt.title("Drug Sensitivity Across Target Pathways")

plt.show()


# The mean IC50 values were plotted against the pathways. The mitosis pathway shows the lowest IC50 value, explaining that mitosis-inhibiting drugs are more effective. While JNK and p38 signalling pathways show a higher IC50 value (i.e., JNK and p38 signalling pathway drugs show resistance).

# In[56]:


heatmap_data = gdsc.pivot_table(
    values='LN_IC50',
    index='TCGA_DESC',
    columns='TARGET_PATHWAY',
    aggfunc='mean'
)


# In[57]:


plt.figure(figsize=(14,8))

sns.heatmap(
    heatmap_data,
    cmap='coolwarm',
    linewidths=0.5
)

plt.title("Cancer Type vs Pathway Drug Sensitivity")

plt.show()


# A heatmap was used to visualise the drug responses in the cancer types against pathways. This clearly shows that in most of the cancer types, the mitosis pathway targeting drugs work better (are more sensitive) than other pathway targeting drugs. We can also see that in the cancer LCML, ABL signal-targeting drugs have a more sensitive response than other drugs.

# ## Pathway-Level Ranking

# In[58]:


pathway_summary = gdsc.groupby(
    'TARGET_PATHWAY'
).agg({
    'LN_IC50': ['mean', 'median', 'std', 'count']
})

pathway_summary


# Mitosis pathway targeting drugs are more sensitive (mean IC50 = -1.39), while metabolism targeting drugs are more variable with a higher std =3.9.  

# In[121]:


mit = gdsc[gdsc["TARGET_PATHWAY"] == "Mitosis"]

mit_sorted = mit.groupby("DRUG_NAME")[["LN_IC50", "AUC"]].median().reset_index().sort_values(by="LN_IC50")

print(mit_sorted)


# **INTERPRETATION**
# 
# Targeting the **Mitosis** pathway is biologically unique. Unlike target-specific kinase therapies that require specific genetic mutations to work, mitosis inhibitors target the physical machinery of cell division.
# 
# The Broad-Spectrum Cluster **(Docetaxel, Vinblastine, Paclitaxel, Vinorelbine)** disrupt tubulin structures to halt cell division mechanically. This translates to the universal, high-potency blue bands on the left side of the pan-cancer landscape heatmap.
# 
# The Selective Cluster **(Alisertib, Tozasertib, ZM447439)** targets small molecules that inhibit regulatory cell-cycle enzymes (Aurora Kinases). Because they rely on specific signaling vulnerabilities rather than physical force, they require higher concentrations across the panel, resulting in positive median values and the widespread orange/red baseline resistance seen across the rest of the heatmap.

# In[ ]:




