from typing import Set

import numpy as np
import pandas as pd
from peptacular.protein import identify_cleavage_sites
from scipy.stats import norm

from constants import BASELINE_VERTEBRATES_AA_FREQUENCY


def get_mc_df(peptides: Set[str], enzyme: str, n:int =1) -> pd.DataFrame:
    if n == 1:
        data = {'cut_position': [], 'peptide': [], 'residue': []}
        for peptide in peptides:

            peptide = peptide.split('.')[1]
            sites = identify_cleavage_sites(peptide, enzyme)

            if len(peptide) in sites:
                sites.remove(len(peptide))

                for site in sites:
                    data['cut_position'].append('C')
                    data['peptide'].append(peptide)
                    data['residue'].append(peptide[site - 1])

                    data['cut_position'].append('N')
                    data['peptide'].append(peptide)
                    data['residue'].append(peptide[site])
    elif n == 2:
        data = {'peptide': [], 'residue': []}
        for peptide in peptides:

            peptide = peptide.split('.')[1]
            sites = identify_cleavage_sites(peptide, enzyme)

            if len(peptide) in sites:
                sites.remove(len(peptide))

                for site in sites:
                    data['peptide'].append(peptide)
                    data['residue'].append(peptide[site - 1] + peptide[site])
    else:
        raise ValueError("n must be 1 or 2")
    df = pd.DataFrame(data)
    df['peptide_term'] = 'MC'
    return df


def get_mc_statistics(site_df: pd.DataFrame, residue_frequencies: dict = None) -> pd.DataFrame:
    if residue_frequencies is None:
        residue_frequencies = BASELINE_VERTEBRATES_AA_FREQUENCY

    # Group by term, cut_position, and residue, and count the number of occurrences
    grouped = site_df.groupby(['cut_position', 'residue']).size().reset_index(name='count')

    # Calculate the frequency of each residue per term and cut_position
    grouped['frequency'] = grouped.groupby(['cut_position'])['count'].transform(lambda x: x / x.sum())

    # Fill in missing residues with 0 frequency and count 0
    all_residues = pd.DataFrame([(cut_position, residue)
                                 for cut_position in grouped['cut_position'].unique()
                                 for residue in residue_frequencies.keys()],
                                columns=['cut_position', 'residue'])
    grouped = pd.merge(all_residues, grouped, on=['cut_position', 'residue'], how='left').fillna(0)
    grouped.sort_values(['cut_position', 'residue'], inplace=True)

    # Calculate the expected count and frequency of each residue per term and cut_position (based on baseline_freqs)
    total_counts = grouped.groupby(['cut_position'])['count'].sum().reset_index(name='total_count')
    grouped = pd.merge(grouped, total_counts, on=['cut_position'])
    grouped['expected_frequency'] = grouped['residue'].apply(lambda x: residue_frequencies.get(x, 0))
    grouped['expected_count'] = grouped['expected_frequency'] * grouped['total_count']

    # Normalize observed and expected counts
    grouped['count_normalized'] = grouped.groupby(['cut_position'])['count'].transform(lambda x: x / x.sum())
    grouped['expected_count_normalized'] = grouped.groupby(['cut_position'])['expected_count'].transform(lambda x: x / x.sum())

    # Calculate Z-scores and p-values for each residue
    grouped['z_score'] = (grouped['frequency'] - grouped['expected_frequency']) / np.sqrt((grouped['expected_frequency'] * (1 - grouped['expected_frequency'])) / grouped['total_count'])
    grouped['p_value'] = 2 * (1 - norm.cdf(abs(grouped['z_score'])))  # Two-tailed test

    # Calculate the log2 fold change of the observed frequency over the expected frequency
    grouped['log2_change'] = np.log2(grouped['frequency'] / grouped['expected_frequency'])

    grouped.to_csv('enzyme_site_statistics.csv', index=False)
    return grouped
