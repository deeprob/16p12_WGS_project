#!/bin/python
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42

# Plot OR for gene enrichment in probands only
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
def plot_lines(region, rel, type, color, ci = True):
    # Get relevant lines
    sub_df = df[(df.relationship==rel) & (df.region==region) & (df.type==type)]
    estimate = sub_df['log2OR'].to_list()
    upper = sub_df['log2OR_upper'].to_list()
    lower = sub_df['log2OR_lower'].to_list()
    
    xs = [i for i in range(len(estimate))]

    if ci:
        plt.fill_between(xs, lower, upper, color = color, alpha = 0.08)
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

def custom_legend(colors):
    lines = []
    for i in colors:
        lines.append(Line2D([0], [0], color=i, lw=2))
    return(lines)

pdf = PdfPages('Figures/4_all_regions_proband_only.pdf')

# Plot all regions in one graph
cmaps = ['Blues_r', 'BrBG', 'BrBG_r', 'BuGn_r', 'BuPu_r', 'CMRmap', 'GnBu_r', 'Greens_r', 'Greys_r', 'OrRd_r', 'Oranges_r', 'PRGn', 'PRGn_r', 'PiYG', 'PiYG_r', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd_r',
        'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds_r', 'Spectral', 'Spectral_r', 'YlGnBu_r', 'YlGn_r', 'YlOrBr_r', 'YlOrRd_r', 'afmhot', 'autumn', 'autumn_r', 'bone', 'brg', 'brg_r', 'bwr',
        'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'gist_earth', 'gist_heat', 'gist_ncar', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gray', 'hot',
        'hsv', 'hsv_r', 'inferno', 'jet', 'jet_r', 'magma', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'pink', 'plasma', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'summer', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r',
        'tab20c', 'tab20c_r', 'terrain', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'winter', 'winter_r']

for cmap in cmaps:
    print(cmap)
    
    plt.plot([0, len(times)-1], [0, 0], linestyle = ':', color = 'k')
    map = cm.get_cmap(cmap, 17)
    for i, region in enumerate(regions):
        plot_lines(region, 'proband', 'All', map(i))
    ax = plt.gca()
    ax.set_xticks([i for i in range(len(times))])
    ax.set_xticklabels(times, rotation = 90)
    # Legend
    custom_lines = custom_legend([map(i) for i in range(16)])
    plt.legend(custom_lines, regions, bbox_to_anchor = (1, 1))
    plt.title(cmap)
    pdf.savefig()
    plt.close()

    # Plot without CI
    plt.plot([0, len(times)-1], [0, 0], linestyle = ':', color = 'k')
    for i, region in enumerate(regions):
        plot_lines(region, 'proband', 'All', map(i), ci = False)
    ax = plt.gca()
    ax.set_xticks([i for i in range(len(times))])
    ax.set_xticklabels(times, rotation = 90)
    # Legend
    custom_lines = custom_legend([map(i) for i in range(16)])
    plt.legend(custom_lines, regions, bbox_to_anchor = (1, 1))
    plt.title(cmap)
    pdf.savefig()
    plt.close()

pdf.close()