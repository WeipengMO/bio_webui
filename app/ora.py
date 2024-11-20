import streamlit as st
import decoupler as dc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO
import base64
import os

plt.rcParams["font.family"] = "Arial"
plt.rcParams['svg.fonttype'] = 'none'


@st.cache_data(ttl='1d')
def perform_ora(genes, gene_sets, threshold):
    genes = genes.split()
    genes = [gene.strip() for gene in genes]
    enr_pvals = dc.get_ora_df(
        df=genes,
        net=gene_sets,
        source='geneset',
        target='genesymbol'
    ).sort_values('FDR p-value', ascending=True)
    enr_pvals = enr_pvals.query('`FDR p-value` < @threshold')
    return enr_pvals


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
    unique_genesets = _msigdb['collection'].unique()

    return _msigdb, unique_genesets


def get_user_inputs(unique_genesets):
    all_option = "Select All"
    options = [all_option] + list(unique_genesets)
    default_collections = ['hallmark', 'kegg_pathways']
    selected_collections = st.multiselect('Select Gene Sets', options, default=default_collections)
    if all_option in selected_collections:
        selected_collections = list(unique_genesets)
    user_genes = st.text_area('Enter Genes (separated by spaces)', height=200)
    pvalue_threshold = st.slider('Set FDR p-value threshold', min_value=0.0, max_value=0.050, value=0.050, step=0.001)
    st.write("The current FDR p-value is ", pvalue_threshold)
    top_n = st.number_input('Number of top results to display', min_value=1, max_value=50, value=10)
    col1, col2 = st.columns(2, vertical_alignment="bottom")
    with col2:
        bar_color = st.color_picker('Color', '#ADD8E6')  # lightblue
    with col1:
        perform_ora_button = st.button('Perform ORA', use_container_width=True)
    return selected_collections, user_genes, pvalue_threshold, top_n, bar_color, perform_ora_button


def plot_results(enr_pvals, top_n, bar_color):
    top_results = enr_pvals.head(top_n).sort_values('FDR p-value', ascending=False)
    plt.figure(figsize=(10, 0.6 * len(top_results)))
    plt.barh(top_results['Term'], -np.log10(top_results['FDR p-value']), color=bar_color)
    plt.xlabel('-log10(FDR p-value)', fontsize=12)
    plt.ylabel('Term', fontsize=12)
    plt.title(f'Top {top_n} Enriched Gene Sets')
    ax = plt.gca()
    for text in ax.get_yticklabels():
        plt.text(text.get_position()[0], text.get_position()[1], text.get_text(), ha='left', va='center')
    ax.set_yticks([])
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
    # setup_ui()
    _msigdb, unique_genesets = load_msigdb()
    selected_collections, user_genes, pvalue_threshold, top_n, bar_color, perform_ora_button = get_user_inputs(unique_genesets)
    if perform_ora_button or bar_color:
        if user_genes and selected_collections:
            gene_sets = _msigdb[_msigdb['collection'].isin(selected_collections)]
            enr_pvals = perform_ora(user_genes, gene_sets, pvalue_threshold)
            if not enr_pvals.empty:
                st.subheader('ORA Results')
                tab1, tab2 = st.tabs(["View as Table", "Plot Results"])
                with tab1:
                    st.dataframe(enr_pvals)
                with tab2:
                    plot_results(enr_pvals, top_n, bar_color)
            else:
                st.warning('No significant results found.')
        else:
            st.warning('Please enter some genes and select gene sets to analyze.')


if __name__ == "__main__":
    main()