from collections import Counter
from typing import List

import pandas as pd
from filterframes import from_dta_select_filter
from peptacular.peptide import strip_modifications
from peptacular.protein import identify_cleavage_sites
from scipy.stats import chi2_contingency
import numpy as np
from constants import BASELINE_VERTEBRATES_AA_FREQUENCY


def calculate_frequency_change(baseline_freq, observed_freq):
    changes = {}
    for amino_acid, base_freq in baseline_freq.items():
        if amino_acid in observed_freq:
            observed = observed_freq[amino_acid]
            change = ((observed - base_freq) / base_freq) * 100
            changes[amino_acid] = change
        else:
            changes[amino_acid] = 0
    return changes


def categorize_amino_acid(amino_acid):
    # Define the categorization dictionary with single characters for each category
    categorization = {
        'A': 'H',
        'C': 'H',
        'D': '-',
        'E': '-',
        'F': 'H',
        'G': 'H',
        'H': '+',
        'I': 'H',
        'K': '+',
        'L': 'H',
        'M': 'H',
        'N': 'P',
        'P': 'H',
        'Q': 'P',
        'R': '+',
        'S': 'P',
        'T': 'P',
        'V': 'H',
        'W': 'H',
        'Y': 'P'
    }

    # Return the category character for the given amino acid
    return categorization.get(amino_acid, 'U')

def convert_sequences_to_categories(sequence):
    categorized_sequence = [categorize_amino_acid(aa) for aa in sequence]
    return ''.join(categorized_sequence)


def explore_mc(peptides: List[str], c_term=False, n=1, enzyme='([KR])', aa_categories=False) -> pd.DataFrame:
    mc = []
    for peptide in peptides:
        peptide = peptide.split('.')[1]
        sites = identify_cleavage_sites(peptide, enzyme)

        if len(peptide) in sites:
            sites.remove(len(peptide))

        if n == 1:
            if c_term:
                for site in sites:
                    mc.append(peptide[site - 1])
            else:

                for site in sites:
                    mc.append(peptide[site])

        elif n == 2:
            for site in sites:
                mc.append(peptide[site - 1] + peptide[site])

    if n == 1:
        baseline_freqs = BASELINE_VERTEBRATES_AA_FREQUENCY

    elif n == 2:
        baseline_freqs = {}
        for aa1 in BASELINE_VERTEBRATES_AA_FREQUENCY:
            for aa2 in BASELINE_VERTEBRATES_AA_FREQUENCY:
                baseline_freqs[aa1 + aa2] = BASELINE_VERTEBRATES_AA_FREQUENCY[aa1] * \
                                            BASELINE_VERTEBRATES_AA_FREQUENCY[
                                                aa2]


    if aa_categories:
        mc = [convert_sequences_to_categories(aa) for aa in mc]
        tmp = {}
        for aa in baseline_freqs:
            converted_aa = convert_sequences_to_categories(aa)
            if converted_aa not in tmp:
                tmp[converted_aa] = 0
            tmp[converted_aa] += baseline_freqs[aa]
        baseline_freqs = tmp

    term_count = Counter(mc)

    # Normalize the frequencies
    total = sum(term_count.values())
    term_freq = {k: v / total for k, v in term_count.items()}
    changes = calculate_frequency_change(baseline_freqs, term_freq)

    # Convert percent changes to log2 fold changes
    log2_changes = {k: np.log2(v / 100 + 1) for k, v in changes.items()}

    df = pd.DataFrame.from_dict(log2_changes, orient='index', columns=['log2_change'])
    df['observed_frequency'] = [term_freq[aa] if aa in term_count else 0.0 for aa in df.index]
    df['observed_count'] = [int(term_count[aa]) if aa in term_count else 0 for aa in df.index]
    df['expected_frequency'] = [baseline_freqs[aa] for aa in df.index]
    df['expected_count'] = [baseline_freqs[aa] * total for aa in df.index]

    chi2_stats = []
    p_values = []

    # Iterate over each row in the DataFrame
    for _, row in df.iterrows():
        # Create a contingency table with observed and expected frequencies
        observed = [row['observed_count'], total - row['observed_count']]
        expected = [row['expected_count'], total - row['expected_count']]

        # Perform chi-square test of independence
        chi2, p, _, _ = chi2_contingency([observed, expected])

        # Append the test statistics and p-values to the respective lists
        chi2_stats.append(chi2)
        p_values.append(p)

    # Add the chi-square test statistics and p-values to the DataFrame
    df['chi2_statistic'] = chi2_stats
    df['p_value'] = p_values

    return df


def explore_sites(peptides: List[str], term='both', n=1, c_term=True, aa_categories=False) -> pd.DataFrame:
    if n == 1:
        if term == 'c':
            if c_term:
                terms = [peptide[-3] for peptide in peptides if peptide[-1] != '-']
            else:
                terms = [peptide[-1] for peptide in peptides if peptide[-1] != '-']

        elif term == 'n':
            if c_term:
                terms = [peptide[0] for peptide in peptides if peptide[0] != '-']
            else:
                terms = [peptide[2] for peptide in peptides if peptide[2] != '-']
        elif term == 'both':
            if c_term:
                n_terms = [peptide[0] for peptide in peptides if peptide[0] != '-']
                c_terms = [peptide[-3] for peptide in peptides if peptide[-1] != '-']
            else:
                n_terms = [peptide[2] for peptide in peptides if peptide[0] != '-']
                c_terms = [peptide[-1] for peptide in peptides if peptide[-1] != '-']
            terms = n_terms + c_terms
        else:
            raise ValueError('term must be one of "c", "n", or "both"')

        baseline_freqs = BASELINE_VERTEBRATES_AA_FREQUENCY
    elif n == 2:
        if term == 'c':
            terms = [peptide[-3] + peptide[-1] for peptide in peptides if peptide[-1] != '-']
        elif term == 'n':
            terms = [peptide[0] + peptide[2] for peptide in peptides if peptide[0] != '-']
        elif term == 'both':
            n_terms = [peptide[0] + peptide[2] for peptide in peptides if peptide[0] != '-']
            c_terms = [peptide[-3] + peptide[-1] for peptide in peptides if peptide[-1] != '-']
            terms = n_terms + c_terms
        else:
            raise ValueError('term must be one of "c", "n", or "both"')

        baseline_freqs = {}
        for aa1 in BASELINE_VERTEBRATES_AA_FREQUENCY:
            for aa2 in BASELINE_VERTEBRATES_AA_FREQUENCY:
                baseline_freqs[aa1 + aa2] = BASELINE_VERTEBRATES_AA_FREQUENCY[aa1] * BASELINE_VERTEBRATES_AA_FREQUENCY[
                    aa2]

    else:
        raise ValueError('n must be 1 or 2')

    if aa_categories:
        terms = [convert_sequences_to_categories(aa) for aa in terms]
        tmp = {}
        for aa in baseline_freqs:
            converted_aa = convert_sequences_to_categories(aa)
            if converted_aa not in tmp:
                tmp[converted_aa] = 0
            tmp[converted_aa] += baseline_freqs[aa]
        baseline_freqs = tmp

    term_count = Counter(terms)

    # Normalize the frequencies
    total = sum(term_count.values())
    term_freq = {k: v / total for k, v in term_count.items()}
    changes = calculate_frequency_change(baseline_freqs, term_freq)

    # Convert percent changes to log2 fold changes
    log2_changes = {k: np.log2(v / 100 + 1) for k, v in changes.items()}

    df = pd.DataFrame.from_dict(log2_changes, orient='index', columns=['log2_change'])
    df['observed_frequency'] = [term_freq[aa] if aa in term_count else 0.0 for aa in df.index]
    df['observed_count'] = [int(term_count[aa]) if aa in term_count else 0 for aa in df.index]
    df['expected_frequency'] = [baseline_freqs[aa] for aa in df.index]
    df['expected_count'] = [baseline_freqs[aa] * total for aa in df.index]

    chi2_stats = []
    p_values = []

    # Iterate over each row in the DataFrame
    for _, row in df.iterrows():
        # Create a contingency table with observed and expected frequencies
        observed = [row['observed_count'], total - row['observed_count']]
        expected = [row['expected_count'], total - row['expected_count']]

        # Perform chi-square test of independence
        chi2, p, _, _ = chi2_contingency([observed, expected])

        # Append the test statistics and p-values to the respective lists
        chi2_stats.append(chi2)
        p_values.append(p)

    # Add the chi-square test statistics and p-values to the DataFrame
    df['chi2_statistic'] = chi2_stats
    df['p_value'] = p_values

    return df


import plotly.graph_objects as go


def plot_volcano(df):
    # Generate a column for -log10(p-value)
    df['minus_log10_p_value'] = -np.log10(df['p_value'] + 1e-100)

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
        text=significant_points.index,  # Use amino acid labels as text
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
        text=non_significant_points.index,  # Use amino acid labels as text
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


def plot_bar(df):
    fig = go.Figure()

    fig.add_trace(go.Bar(
        x=df.index,
        y=df['observed_frequency'],
        name='Observed Frequency'
    ))

    fig.add_trace(go.Bar(
        x=df.index,
        y=df['expected_frequency'],
        name='Expected Frequency'
    ))

    fig.update_layout(
        title='Comparison of Observed and Expected Frequencies',
        xaxis_title='Amino Acids',
        yaxis_title='Frequency',
        barmode='group',
        legend=dict(title='Frequency Type')
    )

    return fig


def plot_log2fold_change(df):
    fig = go.Figure()

    fig.add_trace(go.Bar(x=df.index, y=df['log2_change']))

    fig.update_layout(
        title='Log2 Fold Change',
        xaxis_title='Amino Acids',
        yaxis_title='Log2 Fold Change',
        showlegend=False,
        height=500,
        width=800
    )

    return fig
