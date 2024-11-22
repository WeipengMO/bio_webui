import streamlit as st
import decoupler as dc
import os
import pandas as pd


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