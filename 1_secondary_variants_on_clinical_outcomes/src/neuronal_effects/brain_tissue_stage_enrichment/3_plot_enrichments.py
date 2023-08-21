#!/bin/python
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages

# Plot OR for gene enrichment in probands and carrier parents
df = pd.read_csv('Result_tables/2_fishers_results.csv')

regions = ['cerebellum', 'hippocampus', 'amygdala', 'striatum', 'thalamus', 'medial frontal cortex', 'orbitofrontal cortex',
           'dorsolateral frontal cortex', 'ventrolateral frontal cortex', 'primary motor cortex',
           'primary somatosensory cortex', 'inferior parietal cortex', 'primary auditory cortex',
           'superior temporal cortex', 'inferior temporal cortex', 'primary visual cortex']
times = ['early fetal', 'early mid-fetal', 'late mid-fetal', 'late fetal', 'early infancy', 'late infancy',
         'early childhood', 'middle and late childhood', 'adolescence', 'young adulthood', 'middle adulthood']
rels = ['proband', 'carrier_parent', 'noncarrier_parent']
rel_colors = {'carrier_parent':'#FF5757', 'noncarrier_parent':'#559291', 'proband':'k'}
types = ['missense', 'lof', 'splice', 'DEL', 'DUP', 'STR']
type_colors = dict(zip(types, [(0.83765537, 0.30784814, 0.3723105), (0.66635279, 0.22717328, 0.43008427), (0.90732341, 0.4939774, 0.38990532), (0.253935, 0.265254, 0.529983), '#4c2a85',
                    (0.5494734341032215, 0.7473587230220409, 0.6345610675232811)]))

# Function to plot a region
def plot_lines(region, rel, type, color):
    # Get relevant lines
    sub_df = df[(df.relationship==rel) & (df.region==region) & (df.type==type)]
    estimate = sub_df['log2OR'].to_list()
    upper = sub_df['log2OR_upper'].to_list()
    lower = sub_df['log2OR_lower'].to_list()
    
    xs = [i for i in range(len(estimate))]

    plt.fill_between(xs, lower, upper, color = color, alpha = 0.2)
    plt.plot(xs, estimate, color = color)
    
    # Add * for p-values
    ps = sub_df['bonferroni_p'].to_list()

    if rel=='proband':
        buffer = 0
    elif rel=='carrier_parent':
        buffer = 0.2
    else:
        buffer = -0.2

    for i, p in enumerate(ps):
        if p < 0.05:
            plt.text(xs[i]+buffer, max(estimate)*1.1, '*', color = color)

# Plot all variants in all regions
pdf = PdfPages('Figures/3_AllVariants_by_region.pdf')

for region in regions:
    plt.plot([0, len(times)-1], [0, 0], linestyle = ':', color = 'k')
    
    for rel in rels:
        plot_lines(region, rel, 'All', rel_colors[rel])
      
    plt.title(region)
    
    ax = plt.gca()
    ax.set_xticks([i for i in range(len(times))])
    ax.set_xticklabels(times, rotation = 90)
    plt.ylim(-2.5, 2)
    
    pdf.savefig()
    plt.close()

pdf.close()

# Plot each class of variants in each region
# Separate by relationship
def custom_legend(colors):
    lines = []
    for i in colors:
        lines.append(Line2D([0], [0], color=i, lw=2))
    return(lines)
    
for rel in rels:
    pdf = PdfPages('Figures/3_'+rel+'_VariantType_by_region.pdf')
    
    for region in regions:
        plt.plot([0, len(times)-1], [0, 0], linestyle = ':', color = 'k')
        
        for type in types:
            plot_lines(region, rel, type, type_colors[type])
        plt.title(rel+' '+region)
        
        ax = plt.gca()
        ax.set_xticks([i for i in range(len(times))])
        ax.set_xticklabels(times, rotation = 90)
        plt.ylim(-8.5, 4)
    
        # Legend
        custom_lines = custom_legend([type_colors[i] for i in types])
        plt.legend(custom_lines, types,
           bbox_to_anchor = (1, 1))
    
        pdf.savefig()
        plt.close()
    pdf.close()
    

