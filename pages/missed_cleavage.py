from io import StringIO

import streamlit as st
from filterframes import from_dta_select_filter
from peptacular.peptide import strip_modifications

from enzymeexplorer import get_enzyme_site_statistics
from mcexplorer import get_mc_df, get_mc_statistics
from visualizations import plot_volcano, plot_frequency_bar, plot_log2fold_change

st.title("Missed Cleavage Site Explorer")

filter_files = st.file_uploader("Upload DtaSelect-filter.txt files", accept_multiple_files=True, type='.txt')
enzyme = st.selectbox('Enzyme', ['([KR])', '([KR]|[^P])'])
use_categories = st.checkbox('Use Categories', value=False)

if st.button('Run'):

    peptides = set()

    if filter_files:
        for file in filter_files:
            file_io = StringIO(file.getvalue().decode("utf-8"))
            _, peptide_df, _, _ = from_dta_select_filter(file_io)
            peptides.update([strip_modifications(peptide) for peptide in peptide_df['Sequence'].tolist()])
    else:
        st.write('No files uploaded')
        st.stop()

    peptides = {peptide for peptide in peptides if '-' not in peptide}
    mc_df = get_mc_df(peptides, enzyme=enzyme)
    mc_df2 = get_mc_df(peptides, enzyme=enzyme, n=2)

    mc_site_df = get_enzyme_site_statistics(mc_df)
    mc_site_df2 = get_enzyme_site_statistics(mc_df2)
    mc_site_df2.sort_values(by='log2_change', inplace=True, ascending=False)

    st.metric(label='Number of peptides', value=len(peptides))

    st.subheader('MC First AA')
    st.plotly_chart(plot_volcano(mc_site_df, term='MC', cut_position='C'), use_container_width=True)
    st.plotly_chart(plot_frequency_bar(mc_site_df, term='MC', cut_position='C'), use_container_width=True)
    st.plotly_chart(plot_log2fold_change(mc_site_df, term='MC', cut_position='C'), use_container_width=True)

    st.subheader('MC Second AA')
    st.plotly_chart(plot_volcano(mc_site_df, term='MC', cut_position='N'), use_container_width=True)
    st.plotly_chart(plot_frequency_bar(mc_site_df, term='MC', cut_position='N'), use_container_width=True)
    st.plotly_chart(plot_log2fold_change(mc_site_df, term='MC', cut_position='N'), use_container_width=True)

    st.subheader('MC Combined')
    st.plotly_chart(plot_volcano(mc_site_df2, term='MC', cut_position=None), use_container_width=True)
    st.plotly_chart(plot_frequency_bar(mc_site_df2, term='MC', cut_position=None), use_container_width=True)
    st.plotly_chart(plot_log2fold_change(mc_site_df2, term='MC', cut_position=None), use_container_width=True)

