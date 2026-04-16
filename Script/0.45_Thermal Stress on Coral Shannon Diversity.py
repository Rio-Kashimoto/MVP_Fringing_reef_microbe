#Author: Rio Kashinmoto
#Last edit: April 15th, 2026
#Project: MVP Coral Microbiome Analysis (Porites & Pocillopora)
#Description:This script plots the daily fringing reef benthic water temperature alongside accumulated thermal stress (Degree Heating Weeks, DHW). It overlays the mean Shannon diversity index for Pocillopora and Porites corals before (2017) and after (2019) a major thermal stress event.

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# 1. Load Temperature Data
df_2017 = pd.read_csv('LTER01_2017_Jan_Dec_DailyMean_Temp_1155_1455.csv')
df_2018 = pd.read_csv('LTER01_2018_Jan_Dec_DailyMean_Temp_1155_1455.csv')
df_2019 = pd.read_csv('LTER01_2019_Jan_Dec_DailyMean_Temp_1155_1455.csv')

df_temp = pd.concat([df_2017, df_2018, df_2019]).reset_index(drop=True)
df_temp['date'] = pd.to_datetime(df_temp['date'])
df_temp = df_temp.sort_values('date').reset_index(drop=True)
df_temp['rolling_temp'] = df_temp['mean_temperature_c'].rolling(window=7, min_periods=1, center=True).mean()

# Calculate DHW (Accumulated Thermal Stress)
df_temp['month'] = df_temp['date'].dt.month
monthly_means = df_temp.groupby('month')['mean_temperature_c'].mean()
MMM = monthly_means.max()
df_temp['hotspot'] = df_temp['mean_temperature_c'] - MMM
df_temp['hotspot_local'] = df_temp['hotspot'].apply(lambda x: x if x >= 0.0 else 0)
df_temp['dhw'] = df_temp['hotspot_local'].rolling(window=84, min_periods=1).sum() / 7.0

# 2. Load and Average Diversity Data
df_div = pd.read_csv('POCPOR2017_2019_alphaDiv.csv')
df_div['Group'] = df_div['Samples'].apply(lambda x: "_".join(x.split('_')[:2]))

# Filter: Revert back to the specific 6 groups preferred
requested_groups = [
    'POCFRG2017_1', 'POCFRG2019_1', 
    'PORFRG2017_0', 'PORFRG2017_1', 'PORFRG2019_0', 'PORFRG2019_1'
]
df_div = df_div[df_div['Group'].isin(requested_groups)]

# Group by Host and Year to get the new averages
df_div['Host_Year'] = df_div['Host'] + "_" + df_div['Year'].astype(str)
summary_df = df_div.groupby('Host_Year')['Shannon'].mean().reset_index()

# Define plotting dates and UPDATED colors (forestgreen, firebrick)
plot_info = {
    'Pocillopora_2017': {'date': pd.to_datetime('2017-09-05'), 'color': 'forestgreen'},
    'Porites_2017':     {'date': pd.to_datetime('2017-09-25'), 'color': 'firebrick'},
    'Pocillopora_2019': {'date': pd.to_datetime('2019-08-05'), 'color': 'forestgreen'},
    'Porites_2019':     {'date': pd.to_datetime('2019-08-25'), 'color': 'firebrick'}
}

# 3. Create Plot
fig, ax1 = plt.subplots(figsize=(14, 7))
ax2 = ax1.twinx()

# Plot Temp on Left Axis (UPDATED TO BLACK)
color_temp = 'black'
ax1.plot(df_temp['date'], df_temp['rolling_temp'], color=color_temp, linewidth=2.5, label='Benthic Temp')
ax1.set_ylabel("Fringing Reef Benthic Water Temp (°C)", fontsize=16, fontweight='bold', color=color_temp, labelpad=10)
# Make actual tick numbers larger and bold
ax1.tick_params(axis='y', labelcolor=color_temp, labelsize=14)
for label in ax1.get_yticklabels():
    label.set_fontweight('bold')

# Plot DHW on Right Axis
color_dhw = '#de2d26'
ax2.plot(df_temp['date'], df_temp['dhw'], color=color_dhw, linewidth=1.5, alpha=0.8, label='Thermal Stress')
ax2.fill_between(df_temp['date'], 0, df_temp['dhw'], color=color_dhw, alpha=0.3)
ax2.set_ylabel("Accumulated Thermal Stress (DHW, °C-weeks)", fontsize=16, fontweight='bold', color=color_dhw, labelpad=10)
# Make actual tick numbers larger and bold
ax2.tick_params(axis='y', labelcolor=color_dhw, labelsize=14)
for label in ax2.get_yticklabels():
    label.set_fontweight('bold')
    
ax2.set_ylim(0, 5)

# Plot Bubbles on Left Axis
size_scale = 1000 
for index, row in summary_df.iterrows():
    h_y = row['Host_Year']
    shannon_val = row['Shannon']
    p_info = plot_info[h_y]
    
    # Find temperature at this exact shifted date
    closest_idx = (df_temp['date'] - p_info['date']).abs().idxmin()
    temp_at_date = df_temp.loc[closest_idx, 'rolling_temp']
    
    # Plot bubble mapped to Shannon diversity
    ax1.scatter(p_info['date'], temp_at_date, 
                s=shannon_val * size_scale, 
                c=p_info['color'], 
                edgecolors='black', linewidth=1.5, alpha=0.85, zorder=10)

# Titles & formatting
plt.title("Impact of Thermal Stress on Averaged Coral Shannon Diversity (2017–2019)", fontsize=16, fontweight='bold', pad=20)

ax1.set_xlabel("Time (Month/Year)", fontsize=16, fontweight='bold', labelpad=10)
ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
# Make the X-axis date labels larger and bold
ax1.tick_params(axis='x', labelsize=14)
for label in ax1.get_xticklabels():
    label.set_fontweight('bold')

# Legends
import matplotlib.patches as mpatches

# Legend 1: Species Colors
poc_patch = mpatches.Patch(color='forestgreen', label='Pocillopora (POC)')
por_patch = mpatches.Patch(color='firebrick', label='Porites (POR)')
legend1 = ax1.legend(handles=[poc_patch, por_patch], loc='lower right', title="Host Coral", fontsize=11, title_fontsize=12)
ax1.add_artist(legend1)

# Legend 2: Bubble Size Reference
legend_shannon_vals = [0.5, 1.0, 1.5]
dummy_scatters = []
for val in legend_shannon_vals:
    dummy_scatters.append(ax1.scatter([], [], s=val * size_scale, c='gray', alpha=0.6, edgecolors='black'))

legend2 = ax1.legend(dummy_scatters, [f"H' = {v}" for v in legend_shannon_vals], 
                     loc='upper right', title="Shannon Div. (Cycle Size)", 
                     labelspacing=1.5, borderpad=1, fontsize=11, title_fontsize=12)

ax1.grid(True, linestyle='--', alpha=0.6)
ax1.set_ylim(min(df_temp['mean_temperature_c']) - 0.5, 31.0)
ax1.set_xlim(pd.to_datetime('2017-01-01'), pd.to_datetime('2019-12-31'))

plt.tight_layout()
plt.savefig("final_updated_bubble_plot_v4.png", dpi=300)
