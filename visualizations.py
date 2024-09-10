from typing import Tuple, Optional

import numpy as np
import pandas as pd

from enzymeexplorer import chisquare_test_for_cut_position
import plotly.graph_objects as go


def plot_term_differences(df: pd.DataFrame) -> Tuple[go.Figure, go.Figure]:
    df = df.copy()

    # Keep only rows that have either N or C in the 'term' column
    df = df[df['peptide_term'].isin(['N', 'C'])]

    # Set up two empty figures
    fig_C = go.Figure()
    fig_N = go.Figure()

    # Define colors for each term
    colors = {'C-term': 'rgba(0, 0, 255, 0.5)',
              'N-term': 'rgba(255, 0, 0, 0.5)'}  # adjust RGBA values for your preference

    # Initialize lists to hold the data for each term
    data = {'C-term': {'C': {'x': [], 'y': []}, 'N': {'x': [], 'y': []}},
            'N-term': {'C': {'x': [], 'y': []}, 'N': {'x': [], 'y': []}}}

    # Calculate chi-squared test statistics
    chi_squared_C = chisquare_test_for_cut_position(df, 'c')
    chi_squared_N = chisquare_test_for_cut_position(df, 'n')

    # Iterate over the unique residues in your dataframe
    for residue in df['residue'].unique():
        # Filter the dataframe to include only rows with this residue
        filtered_df = df[df['residue'] == residue]

        # Further filter the dataframe for each combination of peptide_term and cut_position
        for term in filtered_df['peptide_term'].unique():
            for cut_pos in filtered_df['cut_position'].unique():
                # Compute the mean frequency of this residue for this peptide_term and cut_position
                mean_freq = \
                    filtered_df[(filtered_df['peptide_term'] == term) & (filtered_df['cut_position'] == cut_pos)][
                        'frequency'].mean()

                # Add the data to the appropriate list
                data[term + '-term'][cut_pos]['x'].append(residue)
                data[term + '-term'][cut_pos]['y'].append(mean_freq)

    # Add traces to the figures
    for term in ['C-term', 'N-term']:
        fig_C.add_trace(
            go.Bar(x=data[term]['C']['x'], y=data[term]['C']['y'], name=term, marker_color=colors[term], opacity=0.6))
        fig_N.add_trace(
            go.Bar(x=data[term]['N']['x'], y=data[term]['N']['y'], name=term, marker_color=colors[term], opacity=0.6))

    # Customize layout of the figures
    # Customize layout of the figures
    fig_C.update_layout(
        barmode='overlay',
        title_text='C cut position',
        xaxis_title='Residue',
        yaxis_title='Frequency',
        annotations=[
            dict(
                x=0.5,
                y=1.1,
                showarrow=False,
                text=f'Chi-squared: {chi_squared_C}',
                xref="paper",
                yref="paper"
            )
        ]
    )

    fig_N.update_layout(
        barmode='overlay',
        title_text='N cut position',
        xaxis_title='Residue',
        yaxis_title='Frequency',
        annotations=[
            dict(
                x=0.5,
                y=1.1,
                showarrow=False,
                text=f'Chi-squared: {chi_squared_N}',
                xref="paper",
                yref="paper"
            )
        ]
    )

    return fig_C, fig_N


def plot_volcano(df: pd.DataFrame, term: Optional[str], cut_position: Optional[str]) -> go.Figure:
    df = df.copy()
    # Filter the dataframe to include only rows with the specified term and cut_position
    if 'peptide_term' in df.columns and term is not None:
        df = df[df['peptide_term'] == term.upper()]

    if 'cut_position' in df.columns and cut_position is not None:
        df = df[df['cut_position'] == cut_position.upper()]


    # some rows have a 0.0 p_value, add a small random value to avoid taking the log of 0
    df['p_value'] = df['p_value'] + 1e-100

    # Generate a column for -log10(p-value)
    df['minus_log10_p_value'] = -np.log10(df['p_value'])

    # Create a new column 'significance' based on some threshold on the p_value
    df['significance'] = ['significant' if x < 0.05 else 'not significant' for x in df['p_value']]
    # Plotting
    fig = go.Figure()

    # Plot significant points
    significant_points = df[df['significance'] == 'significant']
    fig.add_trace(go.Scatter(
        x=significant_points['log2_change'],
        y=significant_points['minus_log10_p_value'],
        mode='markers+text',
        marker=dict(color='blue'),
        text=significant_points.residue,  # Use amino acid labels as text
        textposition='top center',
        name='Significant'
    ))

    # Plot non-significant points
    non_significant_points = df[df['significance'] == 'not significant']
    fig.add_trace(go.Scatter(
        x=non_significant_points['log2_change'],
        y=non_significant_points['minus_log10_p_value'],
        mode='markers+text',
        marker=dict(color='gray'),
        text=non_significant_points.residue,  # Use amino acid labels as text
        textposition='top center',
        name='Not Significant'
    ))

    # Plot horizontal line at p=0.05
    fig.add_shape(
        type="line",
        x0=df['log2_change'].min(),
        y0=-np.log10(0.05),
        x1=df['log2_change'].max(),
        y1=-np.log10(0.05),
        line=dict(color='red', dash='dash')
    )

    # Update layout
    fig.update_layout(
        title='Volcano plot',
        xaxis_title='Log2 fold change',
        yaxis_title='-Log10(p-value)',
        showlegend=True,
        legend=dict(title='Significance')
    )

    return fig


def plot_frequency_bar(df: pd.DataFrame, term: Optional[str], cut_position: Optional[str]) -> go.Figure:

    df = df.copy()
    # Filter the dataframe to include only rows with the specified term and cut_position
    if 'peptide_term' in df.columns and term is not None:
        df = df[df['peptide_term'] == term.upper()]

    if 'cut_position' in df.columns and cut_position is not None:
        df = df[df['cut_position'] == cut_position.upper()]

    fig = go.Figure()

    fig.add_trace(go.Bar(
        x=df.residue,
        y=df['frequency'],
        name='Observed Frequency',
        opacity=0.6  # added opacity here
    ))

    fig.add_trace(go.Bar(
        x=df.residue,
        y=df['expected_frequency'],
        name='Expected Frequency',
        opacity=0.6  # and here
    ))

    fig.update_layout(
        title='Comparison of Observed and Expected Frequencies',
        xaxis_title='Amino Acids',
        yaxis_title='Frequency',
        barmode='overlay',
        legend=dict(title='Frequency Type')
    )

    return fig


def plot_log2fold_change(df: pd.DataFrame, term: Optional[str], cut_position: Optional[str]) -> go.Figure:

    df = df.copy()

    # Filter the dataframe to include only rows with the specified term and cut_position
    if 'peptide_term' in df.columns and term is not None:
        df = df[df['peptide_term'] == term.upper()]

    if 'cut_position' in df.columns and cut_position is not None:
        df = df[df['cut_position'] == cut_position.upper()]

    fig = go.Figure()

    fig.add_trace(go.Bar(x=df.residue, y=df['log2_change']))

    fig.update_layout(
        title='Log2 Fold Change',
        xaxis_title='Amino Acids',
        yaxis_title='Log2 Fold Change',
        showlegend=False,
        height=500,
        width=800
    )

    return fig


def plot_radar(df: pd.DataFrame, term: str, cut_position: str) -> go.Figure:
    df = df.copy()

    df = df[(df['peptide_term'] == term.upper()) & (df['cut_position'] == cut_position.upper())]

    fig = go.Figure()

    fig.add_trace(go.Scatterpolar(
        r=df['frequency'],
        theta=df.residue,
        fill='toself',
        name='Observed Frequency'
    ))

    fig.add_trace(go.Scatterpolar(
        r=df['expected_frequency'],
        theta=df.residue,
        fill='toself',
        name='Expected Frequency'
    ))

    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0, max(df['frequency'].max(), df['expected_frequency'].max())]
            )),
        showlegend=True,
        title="Comparison of Observed and Expected Frequencies"
    )

    return fig
