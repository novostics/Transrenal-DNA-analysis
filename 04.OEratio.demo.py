'''
#### Step1, calculate the number of DNA fragments located on the regions of interest
#### for example, ChIP seq, H3K27ac
dic=~/hg19_db/ChIPseq/Version2/All_peak
for tissue in B_cell T_cell NK_cell Liver Neutrophil Erythroblast Megakaryocyte
do
cat sample.info | while read sam ol2snp
do
    cat ${dic}/${tissue}.H3K27ac.allpeaks.bed | awk '{print $1"\t"$2"\t"$3}' | bedtools intersect -a ../01.bed/$sam.singlebase.bed -b - -u | cut -f 1,2,3 | awk -v sam=$sam 'BEGIN{OFS="\t"}END{print sam,NR}'
done > ALL_fragments.H3K27ac.${tissue}.xls
done

#### regions of interest, bed format: ${dic}/${tissue}.H3K27ac.allpeaks.bed
#### sample's bed file: ../01.bed/$sam.singlebase.bed, we use the center location to represent every fragment's location
####
'''

#### Step2, python analysis
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import kruskal, mannwhitneyu
from statannot import add_stat_annotation

# Define the cell types and their respective constants, which were calculated by the simulation results
cell_types = {
    "B_cell": 251076,
    "T_cell": 324224,
    "NK_cell": 284660,
    "Liver": 541587,
    "Neutrophil": 168674,
    "Placenta": 453191,
    "Erythroblast": 117378,
    "Megakaryocyte.spe": 27362
}

# Define the total number for the constant denominator, which were calculated by the simulation results
total_number = 13108692

# Read the labels file
labels_file = 'sample.info' ### sample name and label
labels_df = pd.read_csv(labels_file, sep='\t', header=None, names=['sampleID', 'label', 'ot1', 'ot2'])

# Read the second file
file2 = 'ALL_fragments.All.xls' ## Number of total DNA fragment per sample
df2 = pd.read_csv(file2, sep='\t', header=None, names=['sampleID', 'number'])
df2.set_index('sampleID', inplace=True)

# Dictionary to store Kruskal-Wallis p-values
kruskal_p_values = {}

# Loop through each cell type
for cell_type, constant in cell_types.items():
    # Read the first file for the current cell type
    file1 = f'ALL_fragments.H3K27ac.{cell_type}.xls'
    df1 = pd.read_csv(file1, sep='\t', header=None, names=['sampleID', 'number'])
    df1.set_index('sampleID', inplace=True)

    # Perform the initial division
    initial_ratio = df1['number'] / df2['number']

    # Divide by the given constant to get the OE ratio
    oe_ratio = initial_ratio / (constant / total_number)

    # Convert the OE ratio Series to a DataFrame
    oe_ratio_df = oe_ratio.reset_index()
    oe_ratio_df.columns = ['sampleID', 'OE_ratio']

    # Merge with the labels DataFrame
    final_df = pd.merge(oe_ratio_df, labels_df, on='sampleID')

    # Save the output to a file
    output_file = f'OE_ratio_output.CDK.H3K27ac.{cell_type}.txt'
    final_df['OE_ratio'] = final_df['OE_ratio']
    final_df[['sampleID', 'OE_ratio']].to_csv(output_file, sep='\t', index=False, header=False)

    # Plot the OE ratio using a boxplot
    plt.figure(figsize=(4, 4))
    ax = sns.boxplot(x='label', y='OE_ratio', data=final_df, showfliers=False, order=["Control", "Proteinuria"], palette=["grey", "darkred"],width=0.5)
    plt.title(f'O/E Ratio for {cell_type}-marked \n H3K27ac genomic regions')
    plt.xlabel('Label')
    plt.ylabel('O/E Ratio')
    # plt.ylim(0.9, 2)

    # Perform one-sided Mann-Whitney tests and adjust p-values
    comparisons = [
        ("Control", "Proteinuria")
    ]

    pvalues = []
    for comp in comparisons:
        group1 = final_df[final_df['label'] == comp[0]]['OE_ratio']
        group2 = final_df[final_df['label'] == comp[1]]['OE_ratio']
        stat, p = mannwhitneyu(group1, group2, alternative='greater')
        pvalues.append(p)  # Adjust for one-sided test

    # Add statistical annotations
    add_stat_annotation(ax, data=final_df, x='label', y='OE_ratio',
                        box_pairs=comparisons,
                        perform_stat_test=False, pvalues=pvalues,
                        test_short_name='Mann-Whitney', text_format='full',
                        loc='inside', verbose=2,
                        line_height=0.02, text_offset=1)
    
     # Save the plot to a PDF file
    pdf_output_file = f'OE_ratio_plot.CDK.H3K27ac.{cell_type}.pdf'
    plt.savefig(pdf_output_file)

    # Show the plot
    plt.show()