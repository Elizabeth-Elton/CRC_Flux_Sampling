import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from genNN import loadModel

# Define file paths and parameters
fld = 'results_samplingData'  # Folder where the trained model and results are saved
fld_csv = 'data/sampling_results_rp'
n_nn = '25'  # Neural network identifier
model_path = f'{fld}/projector/model_{n_nn}.keras'

# List of model names
model_names = [
    'base_kras_caf', 'caf_80_HPI', 'caf_80_HK', 'caf_60_GAPDH', 'caf_60_ALDO', 
    'caf_40_LDH', 'caf_20_LDH', 'caf_100_HPI', 'caf_100_G6PDH', 
    'caf_100_PYK', 'caf_100_ENO', 'caf_100_PGAM', 'caf_100_HK',
    'base_kras_crc', 'crc_80_HPI', 'crc_80_HK', 'crc_60_GAPDH', 'crc_60_ALDO', 
    'crc_40_LDH', 'crc_20_LDH', 'crc_100_HPI', 'crc_100_G6PDH', 
    'crc_100_PYK', 'crc_100_ENO', 'crc_100_PGAM', 'crc_100_HK'
]

# Initialize lists to store simulations and labels
sims = []
labels = []
uniqueLabels = []

# Loop through each model and load the data
for model_name in model_names:
    file_path = f'{fld_csv}/{model_name}.csv'  # Adjust path as needed
    s = np.loadtxt(file_path, delimiter=',')  # Load simulation data
    
    # Assign labels for each simulation
    for q in range(s.shape[0]):
        labels.append(model_name)  # Model-level label
        uniqueLabels.append(f"{model_name}_{q}")  # Unique label for each simulation
    
    sims.append(s)

# Combine all simulations into a single array
sims = np.vstack(sims)

# Load the trained model
model = loadModel(model_path)

# Project the simulations into the reduced-dimensional space
points = model.predict(sims)

# Create a DataFrame for visualization
data = pd.DataFrame({
    'Dim1': points[:, 0],
    'Dim2': points[:, 1],
    'ModelLabel': labels,  # Fine-grained labels
    'UniqueLabel': uniqueLabels  # Unique labels for each simulation
})

# Save the DataFrame as a CSV file
output_csv_path = f'{fld}/projection_data_25.csv'
data.to_csv(output_csv_path, index=False)
print(f"Data saved to {output_csv_path}")

# Define a color palette for 14 groups
palette = sns.color_palette("tab20", len(model_names))

# Plot the results, coloring by fine-grained model labels
plt.figure(figsize=(10, 8))
sns.scatterplot(
    data=data,
    x='Dim1',
    y='Dim2',
    hue='ModelLabel',
    palette=palette,  # Use the palette for 14 groups
    s=50,
    alpha=0.8
)
plt.title('Neural Network Projection of Simulations (Fine-Grained Labels)')
plt.xlabel('Dimension 1')
plt.ylabel('Dimension 2')
plt.legend(title='Model Labels', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save the plot
output_plot_path = f'{fld}/projection_plot_25.png'
plt.savefig(output_plot_path)
print(f"Plot saved to {output_plot_path}")

# Show the plot
plt.show()