import pandas as pd
import os

data_dir = "/Users/chun/WTL Dropbox/YenChun Chen/TumorEvolutionLab/Project_Peritoneal_Metastasis/wgs/O36/sarek.cnvkit.grch38/variant_calling/cnvkit"
patient = "O36"

# Define thresholds for aberrant CNVs
aberrant_threshold = 0.5

folders = [folder for folder in os.listdir(
    data_dir) if os.path.isdir(os.path.join(data_dir, folder))]

index_dict = {}
for folder in folders:
    sample = folder.split('_')[0]
    folder_path = os.path.join(data_dir, folder)
    files = [file for file in os.listdir(folder_path) if file.endswith('.cnr')]
    data = os.path.join(folder_path, files[0])

    df = pd.read_csv(data, sep='\t')
    # Exclude chrX and chrY
    df = df[~df['chromosome'].isin(['chrX', 'chrY'])]

    # Calculate length of each segment
    df['length'] = df['end'] - df['start']

    # Identify aberrant segments
    df['is_aberrant'] = df['log2'].abs() > aberrant_threshold

    # Calculate the total length of aberrant segments per chromosome
    aberrant_lengths = df[df['is_aberrant']].groupby(
        'chromosome')['length'].sum().reset_index()
    aberrant_lengths.columns = ['chromosome', 'aberrant_length']

    # Calculate the total length of each chromosome
    total_lengths = df.groupby('chromosome')['length'].sum().reset_index()
    total_lengths.columns = ['chromosome', 'total_length']

    # Merge to get a complete dataframe
    merged_lengths = pd.merge(total_lengths, aberrant_lengths,
                              on='chromosome', how='left').fillna(0)

    # Calculate the percentage of aberrant length per chromosome
    merged_lengths['percentage_aberrant'] = merged_lengths['aberrant_length'] / \
        merged_lengths['total_length']

    # Calculate the mean percentage aberration (weighted GII)
    wGII = merged_lengths['percentage_aberrant'].mean()

    # Store the wGII value
    index_dict[sample] = wGII

# Convert the dictionary to a dataframe
wGII_df = pd.DataFrame.from_dict(index_dict, orient='index')
wGII_df.columns = ['wGII']

# Sort the dataframe by index in ascending order
wGII_df.sort_index(inplace=True)

# Turn index into a column called 'sample'
wGII_df.reset_index(inplace=True)
wGII_df.rename(columns={'index': 'sample'}, inplace=True)

# Save the dataframe to a CSV file
wGII_df.to_csv(os.path.join("results", f"{patient}_wGII_0.5.csv"), index=False)
