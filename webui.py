import streamlit as st


def setup_ui():
    st.set_page_config(page_title="Bio WebUI", page_icon="🤔")


def main():
    with st.sidebar:
        st.header("Bio WebUI")
        api_options = ('ORA', 'Survival')
        app_choice = st.selectbox(
            label="Choose an app to run", 
            index=None,
            options=api_options
        )

        if app_choice == "ORA":
            st.caption(
                "Over Representation Analysis (ORA)")
            from app import ora
            run = ora.main
            title = "ORA Analysis"

        elif app_choice == "Survival":
            st.caption(
                "Survival Analysis")
            from app import survival
            run = survival.main
            title = "Survival Analysis"
    
    if app_choice is not None:
        st.title(title)
        run()
            
            
if __name__ == "__main__":
    setup_ui()
    main()
    with st.sidebar:
        st.markdown("---")
        st.markdown(
            '<h5>Made by&nbsp<img src="https://weipengmo.github.io/media/icon_hubc9cb4cc0bdc761ade9c4cb60858dae8_70724_32x32_fill_lanczos_center_3.png" height="32">&nbsp <a href="https://github.com/WeipengMO">@weipeng</a></h5>',
            unsafe_allow_html=True,
        )