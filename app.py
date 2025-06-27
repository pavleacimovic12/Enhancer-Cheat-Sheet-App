import streamlit as st
import pandas as pd
import numpy as np
import pyarrow.feather as feather
from data_processor import DataProcessor
from visualization import VisualizationGenerator
import os

# Configure page
st.set_page_config(
    page_title="Hall of Fame Enhancers Analysis",
    page_icon="🧬",
    layout="wide"
)

# Initialize data processor
@st.cache_data
def load_data():
    """Load and process all data files"""
    processor = DataProcessor()
    return processor.load_all_data()

# Load data
try:
    enhancer_metadata, peak_data, hof_enhancers = load_data()
    st.success("App is loaded!")
except Exception as e:
    st.error(f"❌ Error loading data: {str(e)}")
    st.stop()

# Main title
st.title("CHEAT SHEET FOR ENHANCER ANALYSIS")
st.markdown("---")

# Filter controls in main area
st.header("🔍 Enhancer Selection and Filters")
st.markdown("Choose an enhancer and apply additional filters to focus your analysis")

# Get unique values for filters from the appropriate sources
unique_enhancers = sorted(hof_enhancers['enhancer_id'].unique()) if 'enhancer_id' in hof_enhancers.columns else []

# Get filter values from the metadata (feather file) since that's where cargo, experiment, and gene info is stored
unique_cargos = []
unique_experiments = [] 
unique_genes = []

if not enhancer_metadata.empty:
    # Filter metadata to only include enhancers that are in our HOF enhancers list
    relevant_metadata = enhancer_metadata[enhancer_metadata['enhancer_id'].isin(unique_enhancers)]
    
    unique_cargos = sorted([x for x in relevant_metadata['cargo'].dropna().unique() if x != ''])
    unique_experiments = sorted([x for x in relevant_metadata['experiment'].dropna().unique() if x != ''])
    unique_genes = sorted([x for x in relevant_metadata['proximal_gene'].dropna().unique() if x != ''])
# Sort cell types numerically by their leading numbers (1-34)
import re
def extract_cell_type_number(cell_type):
    match = re.match(r'^(\d+)', str(cell_type))
    return int(match.group(1)) if match else 999

unique_cell_types = sorted(peak_data['cell_type'].unique(), key=extract_cell_type_number) if not peak_data.empty else []

# Main filter controls in columns
col1, col2, col3, col4, col5 = st.columns(5)

with col1:
    selected_enhancer = st.selectbox(
        "Select Enhancer",
        options=["All"] + unique_enhancers,
        index=0,
        help="Choose a specific enhancer to analyze"
    )

with col2:
    selected_cargo = st.selectbox(
        "Filter by Cargo",
        options=["All"] + unique_cargos,
        index=0,
        help="Filter by experimental cargo type"
    )

with col3:
    selected_experiment = st.selectbox(
        "Filter by Experiment", 
        options=["All"] + unique_experiments,
        index=0,
        help="Filter by experiment identifier"
    )

with col4:
    selected_gene = st.selectbox(
        "Filter by Proximal Gene",
        options=["All"] + unique_genes,
        index=0,
        help="Filter by nearest gene"
    )

with col5:
    selected_cell_type = st.selectbox(
        "Filter by Cell Type",
        options=["All"] + unique_cell_types,
        index=0,
        help="Focus on specific cell type for visualization"
    )

# Display summary metrics
st.markdown("---")
col1, col2, col3 = st.columns(3)
with col1:
    st.metric("Total Enhancers", len(unique_enhancers))
with col2:
    st.metric("Cell Types", len(unique_cell_types))
with col3:
    st.metric("Peak Records", f"{len(peak_data):,}")

# Apply filters using both HOF enhancers and metadata
if selected_enhancer != "All":
    # Filter to specific enhancer
    filtered_enhancers = hof_enhancers[hof_enhancers['enhancer_id'] == selected_enhancer].copy()
    relevant_metadata = enhancer_metadata[enhancer_metadata['enhancer_id'] == selected_enhancer] if not enhancer_metadata.empty else pd.DataFrame()
else:
    # Start with all HOF enhancers
    enhancer_ids_to_include = set(hof_enhancers['enhancer_id'].unique())
    
    # Apply metadata-based filters to narrow down the enhancer list
    if not enhancer_metadata.empty:
        filtered_metadata = enhancer_metadata.copy()
        
        if selected_cargo != "All":
            filtered_metadata = filtered_metadata[filtered_metadata['cargo'] == selected_cargo]
        
        if selected_experiment != "All":
            filtered_metadata = filtered_metadata[filtered_metadata['experiment'] == selected_experiment]
            
        if selected_gene != "All":
            filtered_metadata = filtered_metadata[filtered_metadata['proximal_gene'] == selected_gene]
        
        # Get the enhancer IDs that match the metadata filters
        matching_enhancer_ids = set(filtered_metadata['enhancer_id'].unique())
        enhancer_ids_to_include = enhancer_ids_to_include.intersection(matching_enhancer_ids)
    
    # Filter HOF enhancers to only include those that pass metadata filters
    filtered_enhancers = hof_enhancers[hof_enhancers['enhancer_id'].isin(list(enhancer_ids_to_include))].copy()
    relevant_metadata = enhancer_metadata[enhancer_metadata['enhancer_id'].isin(list(enhancer_ids_to_include))] if not enhancer_metadata.empty else pd.DataFrame()

# Display results
if filtered_enhancers.empty:
    st.warning("⚠️ No enhancers match the selected filters. Please adjust your filter criteria.")
elif selected_enhancer == "All":
    st.info(f"📊 Select a specific enhancer from the dropdown above to view detailed analysis")
    st.markdown("### Available Enhancers")
    
    # Show summary table of all enhancers
    enhancer_summary = []
    for idx, enhancer_row in filtered_enhancers.iterrows():
        enhancer_id = enhancer_row.get('enhancer_id', 'Unknown')
        enhancer_peaks = peak_data[peak_data['enhancer_id'] == enhancer_id]
        
        if not enhancer_peaks.empty:
            enhancer_summary.append({
                'Enhancer ID': enhancer_id,
                'Chromosome': enhancer_peaks.iloc[0]['chr'],
                'Start': f"{enhancer_peaks.iloc[0]['start']:,}",
                'End': f"{enhancer_peaks.iloc[0]['end']:,}",
                'Length (bp)': f"{enhancer_peaks.iloc[0]['end'] - enhancer_peaks.iloc[0]['start']:,}",
                'Cell Types': enhancer_peaks['cell_type'].nunique(),
                'Mean Accessibility': f"{enhancer_peaks['accessibility_score'].mean():.4f}",
                'Max Accessibility': f"{enhancer_peaks['accessibility_score'].max():.4f}"
            })
    
    if enhancer_summary:
        summary_df = pd.DataFrame(enhancer_summary)
        st.dataframe(summary_df, use_container_width=True)
else:
    # Process selected enhancer
    for idx, enhancer_row in filtered_enhancers.iterrows():
        enhancer_id = enhancer_row.get('enhancer_id', 'Unknown')
        
        # Create expandable section for each enhancer
        with st.expander(f"🎯 **{enhancer_id}**", expanded=(len(filtered_enhancers) == 1)):
            
            # 1. METADATA DISPLAY (First) - Use authentic metadata from feather file
            st.markdown("### 📋 Enhancer Metadata")
            
            # Get the actual metadata for this enhancer from the feather file
            enhancer_metadata = relevant_metadata[relevant_metadata['enhancer_id'] == enhancer_id] if not relevant_metadata.empty else pd.DataFrame()
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("**Experimental Information:**")
                if not enhancer_metadata.empty:
                    # Show actual metadata from your feather file
                    for _, meta_row in enhancer_metadata.iterrows():
                        if pd.notna(meta_row.get('cargo')) and meta_row.get('cargo') != '':
                            st.markdown(f"- **Cargo:** {meta_row['cargo']}")
                        if pd.notna(meta_row.get('experiment')) and meta_row.get('experiment') != '':
                            st.markdown(f"- **Experiment Type:** {meta_row['experiment']}")
                        if pd.notna(meta_row.get('proximal_gene')) and meta_row.get('proximal_gene') != '':
                            st.markdown(f"- **Proximal Gene:** {meta_row['proximal_gene']}")
                        break  # Just show the first match for now
                else:
                    st.markdown("- No metadata available for this enhancer")
            
            with col2:
                # Get genomic coordinates from peak data
                enhancer_peaks = peak_data[peak_data['enhancer_id'] == enhancer_id]
                if not enhancer_peaks.empty:
                    chr_info = enhancer_peaks.iloc[0]['chr']
                    start_pos = enhancer_peaks.iloc[0]['start']
                    end_pos = enhancer_peaks.iloc[0]['end']
                    length = end_pos - start_pos
                    
                    st.markdown("**Genomic Location:**")
                    st.markdown(f"- **Chromosome:** {chr_info}")
                    st.markdown(f"- **Start Position:** {start_pos:,}")
                    st.markdown(f"- **End Position:** {end_pos:,}")
                    st.markdown(f"- **Length:** {length:,} bp")
                    
                    # Additional statistics
                    total_cell_types = enhancer_peaks['cell_type'].nunique()
                    total_positions = enhancer_peaks['position_index'].nunique()
                    st.markdown(f"- **Cell Types Analyzed:** {total_cell_types}")
                    st.markdown(f"- **Genomic Positions:** {total_positions}")
            
            st.markdown("---")
            
            # 2. EMBEDDED IMAGING (Second) - Contact Sheets and Neuroglancer
            st.markdown("### 🖼️ Imaging Visualization")
            
            # Get all imaging data for this enhancer from the metadata dataframe
            if not relevant_metadata.empty:
                # Get the metadata row for this specific enhancer
                enhancer_meta_row = relevant_metadata[relevant_metadata['enhancer_id'] == enhancer_id]
                
                if not enhancer_meta_row.empty:
                    # Extract multiple visualization URLs from the authentic metadata
                    meta_row = enhancer_meta_row.iloc[0]
                    image_link = meta_row.get('Image_link', '') if pd.notna(meta_row.get('Image_link')) else ''
                    neuroglancer_1 = meta_row.get('Neuroglancer 1', '') if pd.notna(meta_row.get('Neuroglancer 1')) else ''
                    neuroglancer_3 = meta_row.get('Neuroglancer 3', '') if pd.notna(meta_row.get('Neuroglancer 3')) else ''
                    viewer_link = meta_row.get('neuroglancer_url', '') if pd.notna(meta_row.get('neuroglancer_url')) else ''
                    coronal_mip = meta_row.get('Coronal_MIP', '') if pd.notna(meta_row.get('Coronal_MIP')) else ''
                    sagittal_mip = meta_row.get('Sagittal_MIP', '') if pd.notna(meta_row.get('Sagittal_MIP')) else ''
                else:
                    # No metadata found for this enhancer
                    image_link = neuroglancer_1 = neuroglancer_3 = viewer_link = coronal_mip = sagittal_mip = ''
            else:
                # No metadata available
                image_link = neuroglancer_1 = neuroglancer_3 = viewer_link = coronal_mip = sagittal_mip = ''
            
            # Collect all valid URLs
            imaging_urls = []
            if image_link and image_link.startswith('http'):
                imaging_urls.append(('Contact Sheet', image_link))
            if neuroglancer_1 and neuroglancer_1.startswith('http'):
                imaging_urls.append(('Neuroglancer 1', neuroglancer_1))
            if neuroglancer_3 and neuroglancer_3.startswith('http'):
                imaging_urls.append(('Neuroglancer 3', neuroglancer_3))
            if viewer_link and viewer_link.startswith('http'):
                imaging_urls.append(('Viewer', viewer_link))
            if coronal_mip and coronal_mip.startswith('http'):
                imaging_urls.append(('Coronal MIP', coronal_mip))
            if sagittal_mip and sagittal_mip.startswith('http'):
                imaging_urls.append(('Sagittal MIP', sagittal_mip))
            
            if imaging_urls:
                # Display all available imaging modalities as large embedded viewers
                for i, (title, url) in enumerate(imaging_urls):
                    st.markdown(f"#### {title}")
                    
                    # Embed as large, wide iframe
                    st.markdown(
                        f'<iframe src="{url}" width="100%" height="700" frameborder="0" '
                        f'style="border: 2px solid #0066cc; border-radius: 8px; margin: 10px 0;"></iframe>', 
                        unsafe_allow_html=True
                    )
                    
                    if i < len(imaging_urls) - 1:  # Add separator between viewers
                        st.markdown("---")
            else:
                st.info("No imaging visualizations available for this enhancer")
            
            st.markdown("---")
            
            # 3. PYGENOMETRACKS VISUALIZATION (Third)
            st.markdown("### 📈 Peak Accessibility Profile Across Cell Types")
            
            # Filter peak data for this enhancer
            enhancer_peak_data = peak_data[peak_data['enhancer_id'] == enhancer_id].copy()
            
            if selected_cell_type != "All":
                enhancer_peak_data = enhancer_peak_data[enhancer_peak_data['cell_type'] == selected_cell_type]
                st.info(f"Showing data filtered for cell type: **{selected_cell_type}**")
            
            if not enhancer_peak_data.empty:
                try:
                    # Generate pyGenomeTracks-style visualization
                    viz_generator = VisualizationGenerator()
                    fig = viz_generator.create_peak_visualization(enhancer_peak_data, enhancer_id)
                    st.plotly_chart(fig, use_container_width=True)
                    
                    # Summary statistics
                    st.markdown("### 📊 Accessibility Statistics")
                    
                    cell_types_count = enhancer_peak_data['cell_type'].nunique()
                    max_accessibility = enhancer_peak_data['accessibility_score'].max()
                    mean_accessibility = enhancer_peak_data['accessibility_score'].mean()
                    std_accessibility = enhancer_peak_data['accessibility_score'].std()
                    
                    col1, col2, col3, col4 = st.columns(4)
                    with col1:
                        st.metric("Cell Types", cell_types_count)
                    with col2:
                        st.metric("Max Accessibility", f"{max_accessibility:.4f}")
                    with col3:
                        st.metric("Mean Accessibility", f"{mean_accessibility:.4f}")
                    with col4:
                        st.metric("Std Deviation", f"{std_accessibility:.4f}")
                    
                    # Top cell types by accessibility
                    if cell_types_count > 1:
                        st.markdown("**Top Cell Types by Mean Accessibility:**")
                        top_cell_types = (enhancer_peak_data.groupby('cell_type')['accessibility_score']
                                        .mean()
                                        .sort_values(ascending=False)
                                        .head(5))
                        
                        for i, (cell_type, score) in enumerate(top_cell_types.items(), 1):
                            st.markdown(f"{i}. **{cell_type}**: {score:.4f}")
                    
                except Exception as e:
                    st.error(f"Error generating visualization: {str(e)}")
                    
                    # Fallback: Show data table
                    st.markdown("**Raw Peak Data Preview:**")
                    display_cols = ['cell_type', 'position_index', 'accessibility_score']
                    available_cols = [col for col in display_cols if col in enhancer_peak_data.columns]
                    st.dataframe(
                        enhancer_peak_data[available_cols].head(20),
                        use_container_width=True
                    )
            else:
                st.warning("⚠️ No peak accessibility data available for this enhancer with current filters")
            
            st.markdown("---")

# Footer
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666; padding: 20px;'>
    <strong>Hall of Fame Enhancers Analysis Tool</strong><br>
    Genomic Data Visualization Platform for Enhancer Accessibility Analysis
</div>
""", unsafe_allow_html=True)
