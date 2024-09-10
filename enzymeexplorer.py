from collections import Counter
from typing import List, Tuple, Set

import pandas as pd
from peptacular.protein import identify_cleavage_sites
from scipy.stats import chi2_contingency, chisquare, norm, stats
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


def get_enzyme_site_df(peptide_df: pd.DataFrame, n: int = 1) -> pd.DataFrame:

    tryptic_residues = {'K', 'R'}

    if n == 1:
        data = {'peptide_term':[], 'cut_position': [], 'peptide': [], 'residue': [], 'intensity': [], 'tryptic': []}
        for i, row in peptide_df.iterrows():

            peptide = row['StrippedSequence']
            intensity = row['TotalIntensity']

            c_term_c_cut_position = peptide[-3]
            data['peptide_term'].append('C')
            data['cut_position'].append('C')
            data['peptide'].append(peptide)
            data['residue'].append(c_term_c_cut_position)
            data['intensity'].append(intensity)
            data['tryptic'].append(c_term_c_cut_position in tryptic_residues)

            c_term_n_cut_position = peptide[-1]
            data['peptide_term'].append('C')
            data['cut_position'].append('N')
            data['peptide'].append(peptide)
            data['residue'].append(c_term_n_cut_position)
            data['intensity'].append(intensity)
            data['tryptic'].append(c_term_c_cut_position in tryptic_residues)

            n_term_c_cut_position = peptide[0]
            data['peptide_term'].append('N')
            data['cut_position'].append('C')
            data['peptide'].append(peptide)
            data['residue'].append(n_term_c_cut_position)
            data['intensity'].append(intensity)
            data['tryptic'].append(n_term_c_cut_position in tryptic_residues)


            data['peptide_term'].append('N')
            data['cut_position'].append('N')
            data['peptide'].append(peptide)
            data['residue'].append(peptide[2])
            data['intensity'].append(intensity)
            data['tryptic'].append(n_term_c_cut_position in tryptic_residues)

    elif n == 2:
        data = {'peptide_term': [], 'peptide': [], 'residue': [], 'intensity': [], 'tryptic': []}
        for i, row in peptide_df.iterrows():
            peptide = row['StrippedSequence']
            intensity = row['TotalIntensity']

            data['peptide_term'].append('C')
            data['peptide'].append(peptide)
            data['residue'].append(peptide[-3] + peptide[-1])
            data['intensity'].append(intensity)
            data['tryptic'].append(peptide[-3] in tryptic_residues)


            data['peptide_term'].append('N')
            data['peptide'].append(peptide)
            data['residue'].append(peptide[0] + peptide[2])
            data['intensity'].append(intensity)
            data['tryptic'].append(peptide[0] in tryptic_residues)

    else:
        raise ValueError('n must be 1 or 2')

    df = pd.DataFrame(data)
    return df


def get_enzyme_site_statistics(site_df: pd.DataFrame, residue_frequencies: dict = None) -> pd.DataFrame:

    residues_are_pairs = 'cut_position' not in site_df.columns

    if residue_frequencies is None:
        residue_frequencies = BASELINE_VERTEBRATES_AA_FREQUENCY

    if residues_are_pairs is True:
        tmp = {}
        for aa1 in residue_frequencies:
            for aa2 in residue_frequencies:
                tmp[aa1 + aa2] = residue_frequencies[aa1] * residue_frequencies[aa2]
        residue_frequencies = tmp

    group_columns = ['peptide_term']
    if residues_are_pairs is False:
        group_columns.append('cut_position')

    # Group by term, cut_position, and residue, and count the number of occurrences
    grouped = site_df.groupby(group_columns + ['residue']).size().reset_index(name='count')

    if residues_are_pairs is False:
        # Combine C and N terms
        combined = site_df.copy()
        combined['peptide_term'] = 'CN'
        combined_grouped = combined.groupby(group_columns + ['residue']).size().reset_index(name='count')
        # Concatenate grouped and combined_grouped dataframes
        grouped = pd.concat([grouped, combined_grouped], ignore_index=True)

    # Calculate the frequency of each residue per term and cut_position
    grouped['frequency'] = grouped.groupby(group_columns)['count'].transform(lambda x: x / x.sum())

    # Fill in missing residues with 0 frequency and count 0
    if residues_are_pairs is False:
        all_residues = pd.DataFrame([(term, cut_position, residue) for term in grouped['peptide_term'].unique()
                                     for cut_position in grouped['cut_position'].unique()
                                     for residue in residue_frequencies.keys()],
                                    columns=['peptide_term', 'cut_position', 'residue'])
    else:
        all_residues = pd.DataFrame([(term, residue) for term in grouped['peptide_term'].unique()
                                     for residue in residue_frequencies.keys()],
                                    columns=['peptide_term', 'residue'])
    grouped = pd.merge(all_residues, grouped, on=group_columns + ['residue'], how='left').fillna(0)
    grouped.sort_values(group_columns + ['residue'], inplace=True)

    # Calculate the expected count and frequency of each residue per term and cut_position (based on baseline_freqs)
    total_counts = grouped.groupby(group_columns)['count'].sum().reset_index(name='total_count')
    grouped = pd.merge(grouped, total_counts, on=group_columns)
    grouped['expected_frequency'] = grouped['residue'].apply(lambda x: residue_frequencies.get(x, 0))
    grouped['expected_count'] = grouped['expected_frequency'] * grouped['total_count']

    # Normalize observed and expected counts
    grouped['count_normalized'] = grouped.groupby(group_columns)['count'].transform(lambda x: x / x.sum())
    grouped['expected_count_normalized'] = grouped.groupby(group_columns)['expected_count'].transform(lambda x: x / x.sum())

    # Calculate Z-scores and p-values for each residue
    grouped['z_score'] = (grouped['frequency'] - grouped['expected_frequency']) / np.sqrt((grouped['expected_frequency'] * (1 - grouped['expected_frequency'])) / grouped['total_count'])
    #grouped['p_value'] = 2 * (1 - norm.cdf(abs(grouped['z_score'])))  # Two-tailed test
    grouped['p_value'] = norm.sf(abs(grouped['z_score'])) * 2  # Two-tailed test

    # Calculate the log2 fold change of the observed frequency over the expected frequency
    grouped['log2_change'] = np.log2(grouped['frequency'] / grouped['expected_frequency'])

    grouped.to_csv('enzyme_site_statistics.csv', index=False)
    return grouped


def get_enzyme_site_statistics2(site_df: pd.DataFrame, residue_frequencies: dict = None) -> pd.DataFrame:

    residues_are_pairs = len(site_df['residue'].iloc[0]) == 2

    if residue_frequencies is None:
        residue_frequencies = BASELINE_VERTEBRATES_AA_FREQUENCY

    if residues_are_pairs is True:
        tmp = {}
        for aa1 in residue_frequencies:
            for aa2 in residue_frequencies:
                tmp[aa1 + aa2] = residue_frequencies[aa1] * residue_frequencies[aa2]
        residue_frequencies = tmp

    group_columns = ['peptide_term', 'cut_position']
    # Group by term, cut_position, and residue, and count the number of occurrences
    grouped = site_df.groupby(group_columns + ['residue']).size().reset_index(name='count')

    if residues_are_pairs is False:
        # Combine C and N terms
        combined = site_df.copy()
        combined['peptide_term'] = 'CN'
        combined_grouped = combined.groupby(['cut_position', 'residue']).size().reset_index(name='count')
        # Concatenate grouped and combined_grouped dataframes
        grouped = pd.concat([grouped, combined_grouped], ignore_index=True)

    # Calculate the frequency of each residue per term and cut_position
    grouped['frequency'] = grouped.groupby(group_columns)['count'].transform(lambda x: x / x.sum())

    # Fill in missing residues with 0 frequency and count 0
    if residues_are_pairs is False:
        all_residues = pd.DataFrame([(term, cut_position, residue) for term in grouped['peptide_term'].unique()
                                     for cut_position in grouped['cut_position'].unique()
                                     for residue in residue_frequencies.keys()],
                                    columns=['peptide_term', 'cut_position', 'residue'])
    else:
        all_residues = pd.DataFrame([(term, residue) for term in grouped['peptide_term'].unique()
                                     for residue in residue_frequencies.keys()],
                                    columns=['peptide_term', 'residue'])
    grouped = pd.merge(all_residues, grouped, on=group_columns + ['residue'], how='left').fillna(0)
    grouped.sort_values(group_columns + ['residue'], inplace=True)

    # Calculate the expected count and frequency of each residue per term and cut_position (based on baseline_freqs)
    total_counts = grouped.groupby(group_columns)['count'].sum().reset_index(name='total_count')
    grouped = pd.merge(grouped, total_counts, on=group_columns)
    grouped['expected_frequency'] = grouped['residue'].apply(lambda x: residue_frequencies.get(x, 0))
    grouped['expected_count'] = grouped['expected_frequency'] * grouped['total_count']

    # Normalize observed and expected counts
    grouped['count_normalized'] = grouped.groupby(group_columns)['count'].transform(lambda x: x / x.sum())
    grouped['expected_count_normalized'] = grouped.groupby(group_columns)['expected_count'].transform(lambda x: x / x.sum())

    # Calculate Z-scores and p-values for each residue
    grouped['z_score'] = (grouped['frequency'] - grouped['expected_frequency']) / np.sqrt((grouped['expected_frequency'] * (1 - grouped['expected_frequency'])) / grouped['total_count'])
    #grouped['p_value'] = 2 * (1 - norm.cdf(abs(grouped['z_score'])))  # Two-tailed test
    grouped['p_value'] = norm.sf(abs(grouped['z_score'])) * 2  # Two-tailed test

    # Calculate the log2 fold change of the observed frequency over the expected frequency
    grouped['log2_change'] = np.log2(grouped['frequency'] / grouped['expected_frequency'])

    grouped.to_csv('enzyme_site_statistics.csv', index=False)
    return grouped


def chisquare_test_for_cut_position(site_df: pd.DataFrame, cut_position: str) -> Tuple[float, float]:
    cut_position = cut_position.upper()
    if cut_position == 'C':
        n_term_c_cut_pos_df = site_df[(site_df['peptide_term'] == 'N') & (site_df['cut_position'] == 'C')].sort_values(
            'residue')
        c_term_c_cut_pos_df = site_df[(site_df['peptide_term'] == 'C') & (site_df['cut_position'] == 'C')].sort_values('residue')
        res = stats.chisquare(n_term_c_cut_pos_df['frequency'], c_term_c_cut_pos_df['frequency'])
    elif cut_position == 'N':
        n_term_n_cut_pos_df = site_df[(site_df['peptide_term'] == 'N') & (site_df['cut_position'] == 'N')].sort_values(
            'residue')
        c_term_n_cut_pos_df = site_df[(site_df['peptide_term'] == 'C') & (site_df['cut_position'] == 'N')].sort_values('residue')
        res = stats.chisquare(n_term_n_cut_pos_df['frequency'], c_term_n_cut_pos_df['frequency'])
    else:
        raise ValueError(f'Invalid cut_position: {cut_position}')

    return res
