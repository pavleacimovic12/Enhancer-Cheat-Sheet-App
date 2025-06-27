import pandas as pd
import pyarrow.feather as feather
import numpy as np
import os
import streamlit as st

class DataProcessor:
    def __init__(self):
        self.csv_path = "attached_assets/HOF_enhancers_peak_data_1751042112619.csv"
        self.feather_path = "attached_assets/Enhancer_and_experiment_metadata_1751042144891.feather"
    
    def load_all_data(self):
        """Load all data files and return processed datasets"""
        
        # Load peak data from CSV
        peak_data = self.load_peak_data()
        
        # Load metadata from feather file
        enhancer_metadata = self.load_metadata()
        
        # Extract Hall of Fame enhancers
        hof_enhancers = self.extract_hof_enhancers(enhancer_metadata, peak_data)
        
        return enhancer_metadata, peak_data, hof_enhancers
    
    def load_peak_data(self):
        """Load and process the peak accessibility data"""
        try:
            # Load CSV data
            df = pd.read_csv(self.csv_path)
            
            # Validate expected columns
            expected_cols = ['enhancer_id', 'chr', 'start', 'end', 'cell_type', 'position_index', 'accessibility_score']
            missing_cols = [col for col in expected_cols if col not in df.columns]
            
            if missing_cols:
                st.warning(f"Missing columns in peak data: {missing_cols}")
            
            # Data cleaning and validation
            initial_rows = len(df)
            df = df.dropna(subset=['enhancer_id', 'accessibility_score'])
            
            if len(df) < initial_rows:
                st.info(f"Removed {initial_rows - len(df)} rows with missing critical data")
            
            # Convert data types with error handling
            numeric_cols = ['start', 'end', 'position_index', 'accessibility_score']
            for col in numeric_cols:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
            
            # Remove rows with invalid numeric data
            df = df.dropna(subset=numeric_cols)
            
            # Validate genomic coordinates
            if 'start' in df.columns and 'end' in df.columns:
                invalid_coords = df['start'] >= df['end']
                if invalid_coords.any():
                    st.warning(f"Found {invalid_coords.sum()} records with invalid genomic coordinates")
                    df = df[~invalid_coords]
            
            return df
            
        except FileNotFoundError:
            st.error(f"Peak data file not found: {self.csv_path}")
            return pd.DataFrame()
        except Exception as e:
            st.error(f"Error loading peak data: {str(e)}")
            return pd.DataFrame()
    
    def load_metadata(self):
        """Load enhancer metadata from feather file"""
        try:
            # Try to read the actual feather file
            feather_path = "attached_assets/Enhancer_and_experiment_metadata_1751044039206.feather"
            if os.path.exists(feather_path):
                df = feather.read_feather(feather_path)
                
                # Rename columns to match expected names and filter for Hall of Fame enhancers
                if 'Enhancer_ID' in df.columns:
                    df = df.rename(columns={
                        'Enhancer_ID': 'enhancer_id',
                        'Cargo': 'cargo', 
                        'Experiment_Type': 'experiment',
                        'Proximal_Gene': 'proximal_gene',
                        'Viewer Link': 'neuroglancer_url'
                    })
                    
                    # Filter only Hall of Fame enhancers and remove empty enhancer IDs
                    df = df[df['enhancer_id'].notna() & (df['enhancer_id'] != '')]
                    
                    return df
                else:
                    st.warning("Metadata file missing Enhancer_ID column")
                    return pd.DataFrame()
            else:
                st.warning(f"Metadata file not found: {feather_path}")
                return pd.DataFrame()
            
        except Exception as e:
            st.warning(f"Could not load feather metadata file: {str(e)}")
            return pd.DataFrame()
    
    def extract_hof_enhancers(self, metadata_df, peak_data):
        """Extract the Hall of Fame enhancers and merge with metadata"""
        
        if peak_data.empty:
            return pd.DataFrame()
        
        # Get unique enhancer IDs from peak data (these are our Hall of Fame enhancers)
        hof_enhancer_ids = peak_data['enhancer_id'].unique()
        
        # If we have metadata, merge it
        if not metadata_df.empty and 'enhancer_id' in metadata_df.columns:
            # Filter metadata for HOF enhancers
            hof_metadata = metadata_df[metadata_df['enhancer_id'].isin(hof_enhancer_ids)].copy()
            
            # Ensure all HOF enhancers are represented
            missing_enhancers = set(hof_enhancer_ids) - set(hof_metadata['enhancer_id'])
            
            if missing_enhancers:
                # Create placeholder entries for missing enhancers
                missing_data = []
                for enhancer_id in missing_enhancers:
                    # Get basic info from peak data
                    enhancer_peaks = peak_data[peak_data['enhancer_id'] == enhancer_id]
                    cell_types = enhancer_peaks['cell_type'].unique()
                    
                    missing_data.append({
                        'enhancer_id': enhancer_id,
                        'cargo': 'Not specified',
                        'experiment': 'Not specified',
                        'cell_type': f"{len(cell_types)} cell types",
                        'proximal_gene': 'Not specified',
                        'neuroglancer_url': ''
                    })
                
                missing_df = pd.DataFrame(missing_data)
                hof_metadata = pd.concat([hof_metadata, missing_df], ignore_index=True)
        else:
            # Create comprehensive metadata from peak data analysis
            enhancer_info = []
            for enhancer_id in hof_enhancer_ids:
                enhancer_peaks = peak_data[peak_data['enhancer_id'] == enhancer_id]
                
                if not enhancer_peaks.empty:
                    cell_types = enhancer_peaks['cell_type'].unique()
                    mean_accessibility = enhancer_peaks['accessibility_score'].mean()
                    max_accessibility = enhancer_peaks['accessibility_score'].max()
                    
                    enhancer_info.append({
                        'enhancer_id': enhancer_id,
                        'cargo': 'Not specified',
                        'experiment': 'Hall of Fame Collection',
                        'cell_type': f"{len(cell_types)} cell types analyzed",
                        'proximal_gene': 'Not specified',
                        'neuroglancer_url': '',
                        'chromosome': enhancer_peaks.iloc[0]['chr'],
                        'start_position': int(enhancer_peaks.iloc[0]['start']),
                        'end_position': int(enhancer_peaks.iloc[0]['end']),
                        'length_bp': int(enhancer_peaks.iloc[0]['end'] - enhancer_peaks.iloc[0]['start']),
                        'mean_accessibility': round(mean_accessibility, 6),
                        'max_accessibility': round(max_accessibility, 6),
                        'num_cell_types': len(cell_types),
                        'total_data_points': len(enhancer_peaks)
                    })
            
            hof_metadata = pd.DataFrame(enhancer_info)
        
        return hof_metadata.sort_values('enhancer_id')
    
    def get_enhancer_summary(self, peak_data):
        """Generate comprehensive summary statistics for enhancers"""
        if peak_data.empty:
            return {}
        
        summary = {}
        
        # Basic counts
        summary['total_enhancers'] = peak_data['enhancer_id'].nunique()
        summary['total_cell_types'] = peak_data['cell_type'].nunique()
        summary['total_records'] = len(peak_data)
        
        # Accessibility score statistics
        accessibility_scores = peak_data['accessibility_score']
        summary['max_accessibility'] = accessibility_scores.max()
        summary['min_accessibility'] = accessibility_scores.min()
        summary['mean_accessibility'] = accessibility_scores.mean()
        summary['median_accessibility'] = accessibility_scores.median()
        summary['std_accessibility'] = accessibility_scores.std()
        
        # Genomic span analysis
        if 'start' in peak_data.columns and 'end' in peak_data.columns:
            enhancer_lengths = peak_data.groupby('enhancer_id').apply(
                lambda x: x.iloc[0]['end'] - x.iloc[0]['start']
            )
            summary['mean_enhancer_length'] = enhancer_lengths.mean()
            summary['median_enhancer_length'] = enhancer_lengths.median()
            summary['total_genomic_span'] = enhancer_lengths.sum()
        
        # Cell type distribution
        cell_type_counts = peak_data['cell_type'].value_counts()
        summary['most_common_cell_type'] = cell_type_counts.index[0]
        summary['least_common_cell_type'] = cell_type_counts.index[-1]
        
        # Chromosome distribution
        if 'chr' in peak_data.columns:
            chr_counts = peak_data['chr'].value_counts()
            summary['chromosomes_represented'] = len(chr_counts)
            summary['most_common_chromosome'] = chr_counts.index[0]
        
        return summary
    
    def validate_data_integrity(self, peak_data, metadata):
        """Validate data integrity and consistency"""
        validation_results = {
            'is_valid': True,
            'warnings': [],
            'errors': []
        }
        
        if peak_data.empty:
            validation_results['errors'].append("Peak data is empty")
            validation_results['is_valid'] = False
            return validation_results
        
        # Check for required columns
        required_cols = ['enhancer_id', 'chr', 'start', 'end', 'cell_type', 'accessibility_score']
        missing_cols = [col for col in required_cols if col not in peak_data.columns]
        if missing_cols:
            validation_results['errors'].append(f"Missing required columns: {missing_cols}")
            validation_results['is_valid'] = False
        
        # Check for data consistency
        enhancer_coords = peak_data.groupby('enhancer_id')[['chr', 'start', 'end']].nunique()
        inconsistent_coords = enhancer_coords[(enhancer_coords > 1).any(axis=1)]
        if not inconsistent_coords.empty:
            validation_results['warnings'].append(
                f"Inconsistent genomic coordinates for enhancers: {list(inconsistent_coords.index)}"
            )
        
        # Check accessibility score ranges
        acc_scores = peak_data['accessibility_score']
        if acc_scores.min() < 0:
            validation_results['warnings'].append("Found negative accessibility scores")
        
        if acc_scores.max() > 1:
            validation_results['warnings'].append("Found accessibility scores > 1")
        
        return validation_results
