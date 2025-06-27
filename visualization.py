import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import streamlit as st
from typing import Optional, List, Dict, Any

class VisualizationGenerator:
    def __init__(self):
        # Enhanced color palette inspired by pyGenomeTracks and genomic visualization tools
        self.colors = [
            '#8B0000', '#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', 
            '#FFEAA7', '#DDA0DD', '#98D8C8', '#F7DC6F', '#BB8FCE',
            '#85C1E9', '#F8C471', '#82E0AA', '#F1948A', '#AED6F1',
            '#A9DFBF', '#F9E79F', '#D7BDE2', '#A3E4D7', '#FAD7A0',
            '#CD5C5C', '#20B2AA', '#87CEEB', '#DDA0DD', '#F0E68C',
            '#FFB6C1', '#98FB98', '#87CEFA', '#F4A460', '#DA70D6',
            '#32CD32', '#FF69B4', '#00CED1', '#FF1493', '#00FF7F'
        ]
        
        # Cell type specific color mapping for consistency
        self.cell_type_colors = {}
    
    def get_cell_type_color(self, cell_type: str, index: int) -> str:
        """Get consistent color for cell type"""
        if cell_type not in self.cell_type_colors:
            self.cell_type_colors[cell_type] = self.colors[len(self.cell_type_colors) % len(self.colors)]
        return self.cell_type_colors[cell_type]
    
    def create_peak_visualization(self, peak_data: pd.DataFrame, enhancer_id: str) -> go.Figure:
        """Create a pyGenomeTracks-style visualization for peak accessibility"""
        
        if peak_data.empty:
            return self.create_empty_plot("No peak data available for visualization")
        
        # Get unique cell types and sort them numerically by their leading numbers
        cell_types = peak_data['cell_type'].unique()
        
        # Custom sorting function to extract and sort by leading numbers
        def extract_cell_type_number(cell_type):
            import re
            # Extract the leading number from cell type names like "11_CNU_HYa_GABA"
            match = re.match(r'^(\d+)', str(cell_type))
            return int(match.group(1)) if match else 999
        
        cell_types = sorted(cell_types, key=extract_cell_type_number)
        num_cell_types = len(cell_types)
        
        if num_cell_types == 0:
            return self.create_empty_plot("No cell types found in data")
        
        # Create subplots - one track per cell type (similar to pyGenomeTracks)
        fig = make_subplots(
            rows=num_cell_types,
            cols=1,
            shared_xaxes=True,
            vertical_spacing=0.015,  # Tight spacing like genomic browsers
            subplot_titles=[f"{ct}" for ct in cell_types],
            row_heights=[1] * num_cell_types
        )
        
        # Get genomic coordinates
        chr_name = peak_data.iloc[0]['chr']
        start_pos = int(peak_data.iloc[0]['start'])
        end_pos = int(peak_data.iloc[0]['end'])
        total_length = end_pos - start_pos
        
        # Process data for each cell type
        for i, cell_type in enumerate(cell_types, 1):
            cell_data = peak_data[peak_data['cell_type'] == cell_type].copy()
            
            if not cell_data.empty:
                # Sort by position index for proper genomic ordering
                cell_data = cell_data.sort_values('position_index')
                
                # Create genomic positions based on position_index
                # Map position indices to actual genomic coordinates
                max_pos_idx = cell_data['position_index'].max()
                min_pos_idx = cell_data['position_index'].min()
                
                if max_pos_idx > min_pos_idx:
                    cell_data['genomic_position'] = start_pos + (
                        (cell_data['position_index'] - min_pos_idx) / 
                        (max_pos_idx - min_pos_idx)
                    ) * total_length
                else:
                    # Handle case where all positions are the same
                    cell_data['genomic_position'] = start_pos + total_length / 2
                
                # Get consistent color for this cell type
                color = self.get_cell_type_color(cell_type, i-1)
                
                # Create the accessibility track as a filled area plot (more genomic browser-like)
                fig.add_trace(
                    go.Scatter(
                        x=cell_data['genomic_position'],
                        y=cell_data['accessibility_score'],
                        mode='lines',
                        fill='tozeroy',
                        name=cell_type,
                        line=dict(color=color, width=1),
                        fillcolor=color.replace('rgb', 'rgba').replace(')', ', 0.6)') if 'rgb' in color else color,
                        showlegend=False,
                        hovertemplate=(
                            f"<b>{cell_type}</b><br>" +
                            "Position: %{x:,}<br>" +
                            "Accessibility: %{y:.4f}<br>" +
                            "<extra></extra>"
                        )
                    ),
                    row=i, col=1
                )
                
                # Add peak markers for high accessibility regions
                high_accessibility = cell_data[cell_data['accessibility_score'] > cell_data['accessibility_score'].quantile(0.8)]
                if not high_accessibility.empty:
                    fig.add_trace(
                        go.Scatter(
                            x=high_accessibility['genomic_position'],
                            y=high_accessibility['accessibility_score'],
                            mode='markers',
                            marker=dict(
                                color=color,
                                size=4,
                                symbol='diamond',
                                line=dict(width=1, color='white')
                            ),
                            showlegend=False,
                            hoverinfo='skip'
                        ),
                        row=i, col=1
                    )
                
                # Add a baseline at y=0
                fig.add_hline(
                    y=0, 
                    line_dash="solid", 
                    line_color="lightgray", 
                    line_width=0.5,
                    row=i, col=1
                )
        
        # Calculate the global maximum accessibility score across ALL cell types for this enhancer
        global_max_score = peak_data['accessibility_score'].max()
        y_max = global_max_score * 1.1  # Add 10% padding for consistent scaling
        
        # Calculate appropriate height based on number of cell types
        plot_height = max(500, num_cell_types * 120)
        
        # Update layout to match pyGenomeTracks/genomic browser style
        fig.update_layout(
            title={
                'text': f"<b>Peak Accessibility Profile: {enhancer_id}</b><br><span style='font-size:12px'>{chr_name}:{start_pos:,}-{end_pos:,} ({total_length:,} bp)</span>",
                'x': 0.5,
                'xanchor': 'center',
                'font': {'size': 16, 'family': 'Arial, sans-serif'}
            },
            height=plot_height,
            margin=dict(l=150, r=50, t=100, b=80),  # More space for cell type labels
            plot_bgcolor='white',
            paper_bgcolor='white',
            font=dict(size=11, family='Arial, sans-serif'),
            hovermode='x unified'
        )
        
        # Update x-axis for the bottom subplot only
        fig.update_xaxes(
            title_text=f"<b>Genomic Position ({chr_name})</b>",
            showgrid=True,
            gridcolor='lightgray',
            gridwidth=0.5,
            tickformat=',',
            row=num_cell_types, col=1
        )
        
        # Update y-axes for all subplots with SAME SCALE for easy comparison
        for i in range(1, num_cell_types + 1):
            fig.update_yaxes(
                title_text="<b>Accessibility</b>",
                showgrid=True,
                gridcolor='lightgray',
                gridwidth=0.5,
                row=i, col=1,
                range=[0, y_max],  # Same scale for all cell types
                tickformat='.3f'
            )
        
        # Update x-axes for all but the bottom subplot (remove tick labels)
        for i in range(1, num_cell_types):
            fig.update_xaxes(
                showticklabels=False,
                showgrid=True,
                gridcolor='lightgray',
                gridwidth=0.5,
                row=i, col=1
            )
        
        # Update subplot titles styling
        if hasattr(fig.layout, 'annotations') and fig.layout.annotations:
            for i in range(min(num_cell_types, len(fig.layout.annotations))):
                fig.layout.annotations[i].update(
                    font=dict(size=12, color='#333'),
                    xanchor='left',
                    x=0.02
                )
        
        return fig
    
    def create_multi_enhancer_comparison(self, peak_data: pd.DataFrame, enhancer_ids: List[str]) -> go.Figure:
        """Create comparison visualization for multiple enhancers"""
        
        if peak_data.empty or not enhancer_ids:
            return self.create_empty_plot("No data available for comparison")
        
        # Calculate mean accessibility per enhancer per cell type
        comparison_data = (peak_data[peak_data['enhancer_id'].isin(enhancer_ids)]
                          .groupby(['enhancer_id', 'cell_type'])['accessibility_score']
                          .mean()
                          .reset_index())
        
        if comparison_data.empty:
            return self.create_empty_plot("No data found for selected enhancers")
        
        # Create heatmap
        pivot_data = comparison_data.pivot(index='cell_type', columns='enhancer_id', values='accessibility_score')
        
        fig = go.Figure(data=go.Heatmap(
            z=pivot_data.values,
            x=pivot_data.columns,
            y=pivot_data.index,
            colorscale='Viridis',
            hoveringinfo='z',
            hovertemplate="<b>%{y}</b><br>%{x}<br>Accessibility: %{z:.4f}<extra></extra>"
        ))
        
        fig.update_layout(
            title="Mean Accessibility Comparison Across Enhancers and Cell Types",
            xaxis_title="Enhancer ID",
            yaxis_title="Cell Type",
            height=max(400, len(pivot_data.index) * 25)
        )
        
        return fig
    
    def create_empty_plot(self, message: str = "No data to display") -> go.Figure:
        """Create an empty plot with a message"""
        fig = go.Figure()
        
        fig.add_annotation(
            text=message,
            xref="paper", yref="paper",
            x=0.5, y=0.5,
            showarrow=False,
            font=dict(size=16, color="gray", family="Arial, sans-serif"),
            align="center"
        )
        
        fig.update_layout(
            height=300,
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            plot_bgcolor='white',
            paper_bgcolor='white',
            margin=dict(l=50, r=50, t=50, b=50)
        )
        
        return fig
    
    def create_summary_dashboard(self, peak_data: pd.DataFrame) -> go.Figure:
        """Create comprehensive summary dashboard"""
        
        if peak_data.empty:
            return self.create_empty_plot("No peak data available for summary")
        
        # Create subplots for dashboard
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=[
                'Accessibility Distribution by Cell Type',
                'Enhancer Activity Overview', 
                'Top Enhancers by Mean Accessibility',
                'Genomic Coverage Analysis'
            ],
            specs=[[{"type": "box"}, {"type": "bar"}],
                   [{"type": "bar"}, {"type": "scatter"}]]
        )
        
        # 1. Box plot of accessibility by cell type (top 15 for readability)
        top_cell_types = peak_data['cell_type'].value_counts().head(15).index
        cell_type_data = peak_data[peak_data['cell_type'].isin(top_cell_types)]
        
        for i, ct in enumerate(top_cell_types):
            ct_data = cell_type_data[cell_type_data['cell_type'] == ct]['accessibility_score']
            color = self.get_cell_type_color(str(ct), i)
            fig.add_trace(
                go.Box(
                    y=ct_data, 
                    name=ct, 
                    showlegend=False,
                    marker_color=color,
                    boxmean=True
                ),
                row=1, col=1
            )
        
        # 2. Bar chart of enhancer count by cell type
        cell_type_counts = peak_data.groupby('cell_type')['enhancer_id'].nunique().head(15)
        fig.add_trace(
            go.Bar(
                x=cell_type_counts.index, 
                y=cell_type_counts.values, 
                showlegend=False,
                marker_color='lightblue',
                name="Enhancer Count"
            ),
            row=1, col=2
        )
        
        # 3. Top enhancers by mean accessibility
        mean_acc = peak_data.groupby('enhancer_id')['accessibility_score'].mean().nlargest(20)
        fig.add_trace(
            go.Bar(
                x=mean_acc.index, 
                y=mean_acc.values, 
                showlegend=False,
                marker_color='lightcoral',
                name="Mean Accessibility"
            ),
            row=2, col=1
        )
        
        # 4. Genomic span analysis
        genomic_info = peak_data.groupby('enhancer_id').agg({
            'start': 'first',
            'end': 'first',
            'accessibility_score': 'mean'
        }).reset_index()
        genomic_info['length'] = genomic_info['end'] - genomic_info['start']
        
        fig.add_trace(
            go.Scatter(
                x=genomic_info['length'],
                y=genomic_info['accessibility_score'],
                mode='markers',
                marker=dict(
                    size=8,
                    color=genomic_info['accessibility_score'],
                    colorscale='Viridis',
                    showscale=True,
                    colorbar=dict(title="Mean Accessibility")
                ),
                text=genomic_info['enhancer_id'],
                hovertemplate="<b>%{text}</b><br>Length: %{x:,} bp<br>Mean Accessibility: %{y:.4f}<extra></extra>",
                showlegend=False
            ),
            row=2, col=2
        )
        
        # Update layout
        fig.update_layout(
            height=800, 
            title_text="<b>Hall of Fame Enhancers - Summary Dashboard</b>",
            title_x=0.5,
            showlegend=False
        )
        
        # Update axis labels
        fig.update_xaxes(title_text="Cell Type", row=1, col=1, tickangle=45)
        fig.update_yaxes(title_text="Accessibility Score", row=1, col=1)
        
        fig.update_xaxes(title_text="Cell Type", row=1, col=2, tickangle=45)
        fig.update_yaxes(title_text="Number of Enhancers", row=1, col=2)
        
        fig.update_xaxes(title_text="Enhancer ID", row=2, col=1, tickangle=45)
        fig.update_yaxes(title_text="Mean Accessibility", row=2, col=1)
        
        fig.update_xaxes(title_text="Enhancer Length (bp)", row=2, col=2)
        fig.update_yaxes(title_text="Mean Accessibility", row=2, col=2)
        
        return fig
    
    def create_cell_type_specific_view(self, peak_data: pd.DataFrame, cell_type: str) -> go.Figure:
        """Create detailed view for a specific cell type across all enhancers"""
        
        cell_data = peak_data[peak_data['cell_type'] == cell_type]
        
        if cell_data.empty:
            return self.create_empty_plot(f"No data available for cell type: {cell_type}")
        
        # Group by enhancer and calculate statistics
        enhancer_stats = cell_data.groupby('enhancer_id').agg({
            'accessibility_score': ['mean', 'max', 'std', 'count'],
            'chr': 'first',
            'start': 'first',
            'end': 'first'
        }).reset_index()
        
        # Flatten column names
        enhancer_stats.columns = ['enhancer_id', 'mean_acc', 'max_acc', 'std_acc', 'count', 'chr', 'start', 'end']
        enhancer_stats['length'] = enhancer_stats['end'] - enhancer_stats['start']
        
        # Sort by mean accessibility
        enhancer_stats = enhancer_stats.sort_values('mean_acc', ascending=True)
        
        # Create horizontal bar chart
        fig = go.Figure()
        
        fig.add_trace(
            go.Bar(
                y=enhancer_stats['enhancer_id'],
                x=enhancer_stats['mean_acc'],
                orientation='h',
                marker=dict(
                    color=enhancer_stats['mean_acc'],
                    colorscale='Viridis',
                    showscale=True,
                    colorbar=dict(title="Mean Accessibility")
                ),
                text=[f"{acc:.4f}" for acc in enhancer_stats['mean_acc']],
                textposition='auto',
                hovertemplate=(
                    "<b>%{y}</b><br>" +
                    "Mean Accessibility: %{x:.4f}<br>" +
                    "Max Accessibility: %{customdata[0]:.4f}<br>" +
                    "Std Deviation: %{customdata[1]:.4f}<br>" +
                    "Data Points: %{customdata[2]}<br>" +
                    "Length: %{customdata[3]:,} bp<br>" +
                    "<extra></extra>"
                ),
                customdata=np.column_stack((
                    enhancer_stats['max_acc'],
                    enhancer_stats['std_acc'],
                    enhancer_stats['count'],
                    enhancer_stats['length']
                ))
            )
        )
        
        fig.update_layout(
            title=f"<b>Enhancer Accessibility Profile - {cell_type}</b>",
            xaxis_title="Mean Accessibility Score",
            yaxis_title="Enhancer ID",
            height=max(400, len(enhancer_stats) * 25),
            margin=dict(l=100, r=50, t=80, b=50)
        )
        
        return fig
