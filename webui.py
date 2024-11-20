import streamlit as st


def setup_ui():
    st.set_page_config(page_title="Bio WebUI", page_icon="📊")


def main():
    with st.sidebar:
        st.header("Bio WebUI")
        api_options = ("ORA")
        app_choice = st.selectbox(
            label="Choose an app to run", 
            options=api_options
        )

        if app_choice == "ORA":
            st.caption(
                "Over Representation Analysis (ORA)")

        run_button = st.button('Run')


    if run_button:
        if app_choice == "ORA":
            from app import ora
            run = ora.main
            title = "ORA Analysis"

        st.title(title)
        run()
            
            
if __name__ == "__main__":
    setup_ui()
    main()
    with st.sidebar:
        st.markdown("---")
        st.markdown(
            '<h5>Made by&nbsp<img src="https://raw.githubusercontent.com/WeipengMO/resume/refs/heads/main/assets/media/icon.png?token=GHSAT0AAAAAACFXPCRAE2CXFFXYPVCS2ZZGZZ5UUQQ" height="32">&nbsp <a href="https://github.com/WeipengMO">@weipeng</a></h5>',
            unsafe_allow_html=True,
        )