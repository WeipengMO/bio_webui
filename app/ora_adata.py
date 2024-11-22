import streamlit as st
import os
import scanpy as sc
import pandas as pd
import numpy as np
import decoupler as dc
import matplotlib.pyplot as plt
from .utils import load_msigdb
from PyComplexHeatmap import DotClustermapPlotter


@st.cache_data(ttl=86400)  # Cache data for one day (in seconds)
def load_adata(adata_path):
    '''Load an AnnData object from the specified path.'''
    adata = sc.read_h5ad(adata_path)
    return adata


def get_rank_genes_from_groups(
    _adata, 
    groupby='leiden', 
    key='dea_leiden_filtered',
    n_genes=None,
    print_rank_genes=True,
    return_rank_genes=False):
    '''Get ranked genes from filtered rank genes groups.

    Parameters
    ----------
    _adata : AnnData
        The annotated data matrix.
    groupby : str
        The key of the observation grouping to consider.
    key : str
        The key of the ranking to consider.
    n_genes : int, optional
        Number of genes to retrieve.
    print_rank_genes : bool
        Whether to print the ranked genes.
    return_rank_genes : bool
        Whether to return the ranked genes.
    '''
    rank_genes = {}
    for group in _adata.obs[groupby].cat.categories:
        group_data = sc.get.rank_genes_groups_df(_adata, group=group, key=key).dropna()['names'].values
        if group_data.size > 0:
            rank_genes[group] = group_data[:n_genes] if n_genes else group_data
            if print_rank_genes:
                print(f"{group}: {', '.join(rank_genes[group])}")
    if return_rank_genes:
        return rank_genes
    

@st.cache_data(ttl=86400)  # Cache data for one day
def get_rank_genes(_adata, groupby, key, layer, n_genes=None):
    '''Compute and return ranked genes as a DataFrame.'''
    # Check if rank_genes_groups is already computed
    if key in _adata.uns and _adata.uns[key]['params']['groupby'] == groupby:
        # Skip recomputation if parameters match
        pass
    else:
        sc.tl.rank_genes_groups(_adata, groupby=groupby, method="wilcoxon", layer=layer)
        sc.tl.filter_rank_genes_groups(_adata)
    rank_genes = get_rank_genes_from_groups(
        _adata, groupby=groupby, key=key, n_genes=n_genes,
        print_rank_genes=False, return_rank_genes=True)
    rank_genes_df = pd.DataFrame({k: pd.Series(v) for k, v in rank_genes.items()})
    return rank_genes_df


@st.cache_data(ttl=86400)  # Cache data for one day
def run_ora(
        gene_df: pd.DataFrame, 
        gene_sets: pd.DataFrame, 
        source: str ='geneset', 
        target: str ='genesymbol',
        n_top: int = None):
    '''Run ORA analysis and cache the result.'''
    enrich_res = []
    for col in gene_df.columns:
        genes = list(set(gene_df[col]))
        enr_pvals = dc.get_ora_df(
                    df=genes,
                    net=gene_sets,
                    source=source,
                    target=target
                ).sort_values('FDR p-value', ascending=True)
        
        if n_top is not None:
            enr_pvals = enr_pvals.head(n_top)

        enr_pvals['cell_module'] = col
        enrich_res.append(enr_pvals)
    
    enrich_res = pd.concat(enrich_res)
    return enrich_res


def plot_results(rank_genes, gene_sets, top_n_terms):
    '''Plot ORA results using DotClustermapPlotter.'''
    enrich_res = run_ora(rank_genes, gene_sets, n_top=top_n_terms)
    enrich_res = enrich_res[enrich_res['FDR p-value'] < 0.05]
    # Calculate additional columns for plotting
    enrich_res['-log10(FDR p-value)'] = -np.log10(enrich_res['FDR p-value'])
    enrich_res['log10(Combined score)'] = np.log10(enrich_res['Combined score'])
    enrich_res['log10(Odds ratio)'] = np.log10(enrich_res['Odds ratio'])
    enrich_res['Term'] = enrich_res['Term'].apply(lambda x: ' '.join(x.split('_')[1:]))

    # Plot using DotClustermapPlotter
    fig, ax = plt.subplots()
    dm = DotClustermapPlotter(
        data=enrich_res, 
        x='cell_module', y='Term',
        value='log10(Odds ratio)', c='log10(Odds ratio)', s='-log10(FDR p-value)',
        row_cluster=True,
        col_cluster=False,
        cmap='Blues',
        show_rownames=True,
        show_colnames=True,
        row_names_side='left',
        yticklabels_kws={'labelsize': 8},
        verbose=0,
        legend_gap=15,
    )
    st.pyplot(fig)


@st.cache_data(ttl=86400)  # Cache data for one day
def cached_load_msigdb():
    '''Load MSigDB gene sets and cache the result.'''
    return load_msigdb()


def main():
    '''Main function to run the ORA analysis and display results.'''
    # Use cached load_msigdb
    _msigdb, unique_genesets = cached_load_msigdb()
    all_option = "Select All"
    options = [all_option] + list(unique_genesets)
    default_collections = ['hallmark', 'kegg_pathways']
    selected_collections = st.multiselect('Select Gene Sets', options, default=default_collections)
    if all_option in selected_collections:
        selected_collections = list(unique_genesets)

    gene_sets = _msigdb[_msigdb['collection'].isin(selected_collections)]
    
    adata_path = st.text_area('Input the path of AnnData file:', height=68)

    if os.path.exists(adata_path) and adata_path.endswith('.h5ad'):
        adata = load_adata(adata_path)
        st.info(str(adata).split('\n')[0])

        layer_keys = [None] + list(adata.layers.keys())
        layer_key = st.selectbox('Select a layer (optional):', layer_keys)
        if layer_key:
            adata.X = adata.layers[layer_key]
        
        obs_columns = adata.obs.columns
        group_label = st.selectbox('Select a group label:', obs_columns, index=None)

        top_n_genes = st.number_input('Number of top genes to display:', min_value=1, max_value=50, value=20)
        top_n_terms = st.number_input('Number of top terms to display:', min_value=1, max_value=50, value=10)

        run_buttom = st.button('Run ORA', use_container_width=False)
        

        if group_label and run_buttom:
            st.subheader('Cell type specific genes')
            with st.spinner('Getting cell type specific genes...'):
                rank_genes_df = get_rank_genes(
                    adata, 
                    groupby=group_label, 
                    key='rank_genes_groups_filtered', 
                    layer=layer_key, 
                    n_genes=top_n_genes)
            
            tab1, tab2 = st.tabs(["View as Table", "Plot Results"])
            with tab1:
                st.dataframe(rank_genes_df)
            with tab2:
                placeholder = st.info('Running ORA...')
                plot_results(rank_genes_df, gene_sets, top_n_terms)
                placeholder.empty()


    else:
        st.warning('Please input a valid path to an AnnData file.')