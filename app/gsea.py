import streamlit as st
import decoupler as dc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO
import base64
import os
from .utils import parse_gene_input


plt.rcParams["font.family"] = "Arial"
plt.rcParams['svg.fonttype'] = 'none'


@st.cache_data(ttl='1d')
def perform_gsea(ranked_genes, gene_sets, threshold):   
    # Perform GSEA
    gsea_df = dc.get_gsea_df(
        df=ranked_genes,
        stat='stat',
        net=gene_sets,
        source='geneset',
        target='genesymbol',
    ).sort_values('NES', ascending=False)
    # Apply threshold on FDR p-value
    gsea_df = gsea_df.query('`FDR p-value` < @threshold')

    return gsea_df


@st.cache_data(ttl='1d')
def load_msigdb():
    if os.path.exists('data/msigdb.feather'):
        _msigdb = pd.read_feather('data/msigdb.feather')
    else:
        _msigdb = dc.get_resource('MSigDB')
        if not os.path.exists('data'):
            os.makedirs('data')
            _msigdb.to_feather('data/msigdb.feather')
    
    _msigdb = _msigdb[~_msigdb.duplicated(['geneset', 'genesymbol'])]
    unique_collections = _msigdb['collection'].unique()

    return _msigdb, unique_collections


def get_user_inputs(unique_collections):
    all_option = "Select All"
    options = [all_option] + list(unique_collections)
    default_collections = ['hallmark', 'kegg_pathways']
    selected_collections = st.multiselect('Select Gene Sets', options, default=default_collections)
    if all_option in selected_collections:
        selected_collections = list(unique_collections)

    # Parse user input
    user_genes = st.text_area('Enter Ranked Genes (separated by spaces)', height=200)
    user_genes = parse_gene_input(user_genes, remove_duplicates=False)
    st.write('Number of genes:', len(user_genes))

    user_genes_stat = st.text_area('Enter Ranked Genes Statistic (separated by spaces)', height=200)
    user_genes_stat = parse_gene_input(user_genes_stat, remove_duplicates=False)
    st.write('Number of statistics:', len(user_genes_stat))

    if len(user_genes) != len(user_genes_stat):
        st.error('The number of genes and statistics do not match.')
        return None, None, None, None, None, None
    else:
        ranked_genes = pd.DataFrame(
            index=user_genes,
            data={'stat': user_genes_stat}, 
            columns=['stat'])


        pvalue_threshold = st.slider('Set FDR p-value threshold', min_value=0.0, max_value=0.050, value=0.050, step=0.001)
        st.write("The current FDR p-value is ", pvalue_threshold)
        top_n = st.number_input('Number of top results to display', min_value=1, max_value=50, value=10)
        col1, col2 = st.columns(2, vertical_alignment="bottom")
        with col2:
            bar_color = st.color_picker('Color', '#ADD8E6')  # lightblue
        with col1:
            perform_gsea_button = st.button('Perform GSEA', use_container_width=True)
        return selected_collections, ranked_genes, pvalue_threshold, top_n, bar_color, perform_gsea_button


def plot_results(enr, top_n, bar_color):
    top_results = enr.head(top_n).sort_values('NES', ascending=True)
    plt.figure(figsize=(10, 0.6 * len(top_results)))
    plt.barh(top_results['Term'], top_results['NES'], color=bar_color)
    plt.xlabel('Normalized Enrichment Score (NES)', fontsize=12)
    plt.ylabel('Term', fontsize=12)
    plt.title(f'Top {top_n} Enriched Gene Sets')
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    buf = BytesIO()
    plt.savefig(buf, format='png', dpi=300)
    buf.seek(0)
    image_base64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close()
    st.image(f'data:image/png;base64,{image_base64}')


def main():
    _msigdb, unique_collections = load_msigdb()
    selected_collections, ranked_genes, pvalue_threshold, top_n, bar_color, perform_gsea_button = get_user_inputs(unique_collections)
    if perform_gsea_button:
        if selected_collections:
            gene_sets = _msigdb[_msigdb['collection'].isin(selected_collections)]
            enr = perform_gsea(ranked_genes, gene_sets, pvalue_threshold)
            if not enr.empty:
                st.subheader('GSEA Results')
                tab1, tab2 = st.tabs(["View as Table", "Plot Results"])
                with tab1:
                    st.dataframe(enr)
                with tab2:
                    plot_results(enr, top_n, bar_color)
            else:
                st.warning('No significant results found.')
        else:
            st.warning('Please enter some genes and select gene sets to analyze.')

if __name__ == "__main__":
    main()