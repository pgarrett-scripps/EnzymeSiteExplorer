from io import StringIO

import pandas as pd
import streamlit as st
from filterframes import from_dta_select_filter
from peptacular.peptide import strip_modifications

from enzymeexplorer import get_enzyme_site_df, get_enzyme_site_statistics, chisquare_test_for_cut_position
from visualizations import plot_term_differences, plot_volcano, plot_frequency_bar, plot_log2fold_change, plot_radar

st.title("Enzyme Site Explorer")

filter_files = st.file_uploader("Upload DtaSelect-filter.txt files", accept_multiple_files=True, type='.txt')

use_categories = st.checkbox('Use Categories', value=False)
hide_trypsin = st.checkbox('Hide Trypsin', value=False)
use_intensity = st.checkbox('use intensity', value=False)


dfs = []
if st.button('Run'):
    peptides = []

    if filter_files:
        for file in filter_files:
            file_io = StringIO(file.getvalue().decode("utf-8"))
            _, peptide_df, _, _ = from_dta_select_filter(file_io)
            dfs.append(peptide_df)
            peptides.extend(list({strip_modifications(peptide) for peptide in peptide_df['Sequence'].tolist()}))
    else:
        st.write('No files uploaded')
        st.stop()

    peptide_df = pd.concat(dfs)

    # remove protein terminus peptides
    peptide_df = peptide_df[~peptide_df['Sequence'].str.contains('-')]
    peptide_df['StrippedSequence'] = peptide_df['Sequence'].apply(strip_modifications)

    df = get_enzyme_site_df(peptide_df, n=1)

    if hide_trypsin:
        df = df[~df['tryptic']]

    site_df = get_enzyme_site_statistics(df)


    st.metric(label='Number of peptides', value=len(peptide_df))

    df2 = get_enzyme_site_df(peptide_df, n=2)
    if hide_trypsin:
        df2 = df2[~df2['tryptic']]
    site_df2 = get_enzyme_site_statistics(df2)
    site_df2.sort_values(by='log2_change', inplace=True, ascending=False)

    with st.expander('Single Site Data'):
        st.dataframe(site_df)

    with st.expander('Double Site Data'):
        st.dataframe(site_df2)

    st.markdown("""
    
    ### Cleavage Frequencies
    
    #### C-Terminus
    - **C-Cut**: `x.xxX.x`
    - **N-Cut**: `x.xxx.X`
    
    #### N-Terminus
    - **C-Cut**: `X.xxx.x`
    - **N-Cut**: `x.Xxx.x`

    """)

    st.subheader('C Term')
    st.plotly_chart(plot_volcano(site_df2, term='C', cut_position=None), use_container_width=True)
    st.plotly_chart(plot_frequency_bar(site_df2, term='C', cut_position=None), use_container_width=True)
    st.plotly_chart(plot_log2fold_change(site_df2, term='C', cut_position=None), use_container_width=True)

    st.subheader('N Term')
    st.plotly_chart(plot_volcano(site_df2, term='N', cut_position=None), use_container_width=True)
    st.plotly_chart(plot_frequency_bar(site_df2, term='N', cut_position=None), use_container_width=True)
    st.plotly_chart(plot_log2fold_change(site_df2, term='N', cut_position=None), use_container_width=True)

    st.subheader('C-Term vs N-Term')
    fig_C, fig_N = plot_term_differences(site_df)
    st.plotly_chart(fig_C, use_container_width=True)
    st.plotly_chart(fig_N, use_container_width=True)

    st.subheader('C Cut Position')
    st.plotly_chart(plot_volcano(site_df, term='CN', cut_position='C'), use_container_width=True)
    st.plotly_chart(plot_frequency_bar(site_df, term='CN', cut_position='C'), use_container_width=True)
    st.plotly_chart(plot_log2fold_change(site_df, term='CN', cut_position='C'), use_container_width=True)
    st.plotly_chart(plot_radar(site_df, term='CN', cut_position='C'), use_container_width=True)

    st.subheader('N Cut Position')
    st.plotly_chart(plot_volcano(site_df, term='CN', cut_position='N'), use_container_width=True)
    st.plotly_chart(plot_frequency_bar(site_df, term='CN', cut_position='N'), use_container_width=True)
    st.plotly_chart(plot_log2fold_change(site_df, term='CN', cut_position='N'), use_container_width=True)
    st.plotly_chart(plot_radar(site_df, term='CN', cut_position='N'), use_container_width=True)

    with st.expander('C-Term'):
        st.subheader('C-Term N Cut Position')
        st.plotly_chart(plot_volcano(site_df, term='C', cut_position='N'), use_container_width=True)
        st.plotly_chart(plot_frequency_bar(site_df, term='C', cut_position='N'), use_container_width=True)
        st.plotly_chart(plot_log2fold_change(site_df, term='C', cut_position='N'), use_container_width=True)
        st.plotly_chart(plot_radar(site_df, term='C', cut_position='N'), use_container_width=True)

        st.subheader('C-Term C Cut Position')
        st.plotly_chart(plot_volcano(site_df, term='C', cut_position='C'), use_container_width=True)
        st.plotly_chart(plot_frequency_bar(site_df, term='C', cut_position='C'), use_container_width=True)
        st.plotly_chart(plot_log2fold_change(site_df, term='C', cut_position='C'), use_container_width=True)
        st.plotly_chart(plot_radar(site_df, term='C', cut_position='C'), use_container_width=True)

    with st.expander('N-Term'):
        st.subheader('N-Term N Cut Position')
        st.plotly_chart(plot_volcano(site_df, term='N', cut_position='N'), use_container_width=True)
        st.plotly_chart(plot_frequency_bar(site_df, term='N', cut_position='N'), use_container_width=True)
        st.plotly_chart(plot_log2fold_change(site_df, term='N', cut_position='N'), use_container_width=True)
        st.plotly_chart(plot_radar(site_df, term='N', cut_position='N'), use_container_width=True)

        st.subheader('N-Term C Cut Position')
        st.plotly_chart(plot_volcano(site_df, term='N', cut_position='C'), use_container_width=True)
        st.plotly_chart(plot_frequency_bar(site_df, term='N', cut_position='C'), use_container_width=True)
        st.plotly_chart(plot_log2fold_change(site_df, term='N', cut_position='C'), use_container_width=True)
        st.plotly_chart(plot_radar(site_df, term='N', cut_position='C'), use_container_width=True)
