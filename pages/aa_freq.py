from collections import Counter
from io import StringIO

import numpy as np
import streamlit as st
import pandas as pd
from filterframes import from_dta_select_filter
from peptacular.peptide import strip_modifications
from constants import BASELINE_VERTEBRATES_AA_FREQUENCY

st.title("AA Freq Explorer")

filter_files = st.file_uploader("Upload DtaSelect-filter.txt files", accept_multiple_files=True, type='.txt')

run = st.button('Run')

if run:

    if not filter_files:
        st.write('No files uploaded')
        st.stop()

    for file in filter_files:

        st.subheader(file.name)

        file_io = StringIO(file.getvalue().decode("utf-8"))
        _, peptide_df, _, _ = from_dta_select_filter(file_io)
        peptides = []
        peptides.extend([strip_modifications(peptide[2:-2]) for peptide in peptide_df['Sequence'].tolist()])


        aa_counter = Counter()
        for peptide in peptides:
            aa_counter.update(peptide)  # Update with each peptide string
        total_aa = sum(aa_counter.values())

        aa_freqs = {aa: count / total_aa for aa, count in aa_counter.items()}

        # Prepare DataFrame
        aa_freq_df = pd.DataFrame.from_dict(aa_freqs, orient='index', columns=['Observed Frequency'])
        aa_freq_df['Baseline Frequency'] = aa_freq_df.index.map(BASELINE_VERTEBRATES_AA_FREQUENCY.get)


        # Calculate log2 fold change
        def calculate_log2fold_change(observed_freq, baseline_freq):
            if baseline_freq == 0:
                return np.nan  # Avoid division by zero
            ratio = observed_freq / baseline_freq
            return np.log2(ratio) if ratio > 0 else np.nan  # Handle non-positive ratios


        aa_freq_df['Log2Fold Change'] = aa_freq_df.apply(
            lambda row: calculate_log2fold_change(row['Observed Frequency'], row['Baseline Frequency']), axis=1)
        aa_freq_df = aa_freq_df.sort_values(by='Log2Fold Change', ascending=False)

        # Plotting and displaying data
        st.dataframe(aa_freq_df)
        st.scatter_chart(
            aa_freq_df[['Baseline Frequency', 'Observed Frequency']].sort_values(by='Baseline Frequency', ascending=False))

        st.caption('AA Log2Fold Change')
        st.bar_chart(data=aa_freq_df['Log2Fold Change'])
