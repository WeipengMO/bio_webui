from .utils import _survival
import streamlit as st
from io import BytesIO
import matplotlib.pyplot as plt
import base64
import yaml

    
plt.rcParams["font.family"] = "Arial"
plt.rcParams['svg.fonttype'] = 'none'

SURVIVAL_METRICS = ['OS', 'DSS', 'PFI']


@st.cache_data(ttl='1d')
def load_survival_data(data, survival_data):
    exp_data = data[survival_data]['exp']
    meta_data = data[survival_data]['meta']
    ad_tcga = _survival.Survival(exp_data, meta_data, meta_index_col='sample')

    return ad_tcga

def main():

    with open('data/survival_data.yaml', 'r') as file:
        data = yaml.safe_load(file)

    # Load survival data
    survival_data = st.selectbox(
        'Select Survival Data', 
        list(data.keys()),
        index=None
        )
    if survival_data is not None:
        ad_tcga = load_survival_data(data, survival_data)

        user_genes = st.text_area('Enter Genes (separated by spaces)', height=200)

        # Select survival metrics
        survival_metrics = st.selectbox(
            'Survival metrics', 
            [metric for metric in ad_tcga.obs.columns if metric in SURVIVAL_METRICS]
            )
        ad_tcga.obs['event'] = ad_tcga.obs[survival_metrics].copy()
        ad_tcga.obs['time'] = ad_tcga.obs[f'{survival_metrics}.time'].copy()

        # Select axis units
        axis_units = st.selectbox('Axis Units', ['Days', 'Months'])
        if axis_units == 'Months':
            ad_tcga.obs['time'] /= 30

        # Select the range of time
        max_time = st.number_input(
            f'Time Range ({axis_units})', 
            min_value=0.0, 
            max_value=ad_tcga.obs['time'].max(), 
            value=ad_tcga.obs['time'].max()
            )

    ci_show = st.checkbox('Show 95% Confidence Interval', value=False)

    run_button = st.button('Run', use_container_width=True)


    if run_button:
        genes = user_genes.split()
        genes = [gene.strip() for gene in genes]

        ad_tcga.group_meta(
            groupby=genes, 
            group_method='median', 
            event='event', 
            time='time',
            time_limit=max_time,
        )

        # Plot survival curves
        st.subheader('Survival Plots')
        buf = BytesIO()
        fig, ax = plt.subplots(figsize=(4, 4))
        ad_tcga.km_plot(
            ax=ax, 
            xlabel=axis_units,
            ylabel=survival_metrics,
            ci_show=ci_show
            )
        plt.tight_layout()
        plt.savefig(buf, format='png', dpi=150)
        buf.seek(0)
        image_base64 = base64.b64encode(buf.read()).decode('utf-8')
        st.markdown(
            f'<div style="display: flex; justify-content: center;">'
            f'<img src="data:image/png;base64,{image_base64}" alt="Kaplan-Meier Plot">'
            f'</div>',
            unsafe_allow_html=True
        )
