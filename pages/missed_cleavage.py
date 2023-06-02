from io import StringIO

import streamlit as st
from filterframes import from_dta_select_filter
from peptacular.peptide import strip_modifications

from enzymeexplorer import explore_sites, plot_volcano, plot_bar, plot_log2fold_change, explore_mc

st.title("Enzyme Site Explorer")

filter_files = st.file_uploader("Upload DtaSelect-filter.txt files", accept_multiple_files=True, type='.txt')

single_site = st.checkbox('Single Residue', value=False)
single_site = 1 if single_site else 2
term = st.radio('Peptide Terminus', options=['n', 'c', 'both'])
use_categories = st.checkbox('Use Categories', value=False)
enzyme = st.text_input('Enzyme', value='([KR])')

residue_term = True
if single_site is 1:
    residue_term = st.radio('Cut position', options=['n', 'c'])
    residue_term = residue_term == 'c'

if st.button('Run'):
    peptides = set()

    if filter_files:
        for file in filter_files:
            file_io = StringIO(file.getvalue().decode("utf-8"))
            _, peptide_df, _, _ = from_dta_select_filter(file_io)
            peptides.update({strip_modifications(peptide) for peptide in peptide_df['Sequence'].tolist()})
    else:
        st.write('No files uploaded')
        st.stop()

    peptides.update({strip_modifications(peptide) for peptide in peptide_df['Sequence'].tolist()})

    df = explore_mc(peptides, n=single_site, c_term=residue_term, aa_categories=use_categories, enzyme=enzyme)
    df.sort_values(by='log2_change', inplace=True)

    st.plotly_chart(plot_volcano(df))
    st.plotly_chart(plot_bar(df))
    st.plotly_chart(plot_log2fold_change(df))