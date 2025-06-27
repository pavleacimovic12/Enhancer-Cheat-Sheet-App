# Enhancer-Cheat-Sheet-App
Enhancer Cheat Sheet App with metadata, imaging and peak accessibility visualization across brain cell types


# Hall of Fame Enhancers Analysis - Cheat Sheet

A Streamlit-based genomic data visualization platform for analyzing enhancer accessibility profiles across 55 Hall of Fame enhancers and 34 cell types.

## Features

- **Authentic Data Integration**: Uses real enhancer data from CSV and metadata from Feather files
- **Interactive Filtering**: Filter by enhancer, cargo, experiment type, proximal gene, and cell type
- **Embedded Visualizations**: Large contact sheets and Neuroglancer viewers embedded directly
- **pyGenomeTracks-style Plots**: Peak accessibility profiles with consistent Y-axis scaling
- **Comprehensive Metadata**: Real experimental information from Allen Institute data

## Data Sources

- **Peak Data**: `HOF_enhancers_peak_data_*.csv` - 783,158 accessibility records
- **Metadata**: `Enhancer_and_experiment_metadata_*.feather` - 411 experimental records

## Installation

```bash
# Clone the repository
git clone <your-repo-url>
cd hall-of-fame-enhancers

# Install dependencies
pip install streamlit pandas numpy plotly pyarrow

# Run the application
streamlit run app.py --server.port 5000
