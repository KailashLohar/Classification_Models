import os
import base64
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO, StringIO

from scripts.utils import *
from scripts.custom_styles import *
from scripts.counterfactual_generator import CounterfactualGenerator
from scripts.chemical_space_explorer import ChemicalSpaceExplorer
    
def main():
    
    header_and_footer_color = "#282d56"
    logo_image_path = "scripts/my_logo.svg"
    header_background_path = "scripts/header_background.png"
    background_image = "scripts/back_ground.jpg"
    title = "Property Modelling"
    subtitle = "GenAI powered molecule optimization"
    more_info_url = "https://www.aganitha.ai/services/computational-chemistry/#:~:text=Solution%20area-,ADMET%20prediction,-Our%20solutions%20to"
    
    add_custom_header_and_footer(header_and_footer_color, logo_image_path, header_background_path, 
                                 background_image, title, subtitle, more_info_url)

    add_custom_css()

    initialize_session_state()
    input_option, uploaded_file, single_smiles = handle_input_option()

    model_choice = st.sidebar.selectbox("Select dataset",["in vitro (H-CLAT)",
                                                          "in vitro (KeratinoSens)",
                                                          "in vivo (LLNA)",
                                                          "in chemico (DPRA)",
                                                          "human"],index=0)

    submit_button = st.sidebar.button("Predict")

    if submit_button:
        df = handle_file_upload(input_option, uploaded_file, single_smiles)
        if df is not None:
            st.session_state.file_uploaded = True
            df = process_smiles_data(df)
            update_session_state_with_results(df, model_choice)
        else:
            st.session_state.file_uploaded = False

    if st.session_state.prediction_table:
        st.markdown(load_centered_table_styles(), unsafe_allow_html=True)
        prediction_title = model_choice.split("skin_sensitization_")[-1]
        st.markdown(f"<span style='color:#8e5572; font-weight:bold; font-size:22px; margin-top:-20px; margin-bottom: 0px;'>{prediction_title} predictions:</span>", unsafe_allow_html=True)
        st.components.v1.html(st.session_state.prediction_table, height=500, scrolling=True)
        st.markdown("<div style='margin-top: 10px; margin-bottom: 10px; color:#8e5572; font-weight: bold; font-size: 22px;'>What would you like to do next?</div>", unsafe_allow_html=True)

    if st.session_state.file_uploaded and st.session_state.gnn_results is not None:
        tabs = st.tabs(["Generate counterfactuals for sensitizers",
                        "Explore chemical space for non-sensitizers"])
    
        with tabs[0]:
            render_tab(
                tab_type="Sensitizer",
                smiles_list=st.session_state.sensitizer_smiles,
                model_choice=model_choice,
                button_label="Generate counterfactuals",
                button_key="counterfactual_button",
                action_function=lambda selected_smiles: CounterfactualGenerator(
                    smiles=selected_smiles, model_choice=model_choice
                ).generate_and_display_counterfactuals(),
                header_text="Select a sensitizer SMILES for counterfactual generation",
                spinner_text="Generating counterfactuals...",
            )
            st.write("")
            st.write("")
            st.write("")

        with tabs[1]:
            render_tab(
                tab_type="Non-sensitizer",
                smiles_list=st.session_state.non_sensitizer_smiles,
                model_choice=model_choice,
                button_label="Explore chemical space",
                button_key="similar_smiles_button",
                action_function=lambda selected_smiles: ChemicalSpaceExplorer(
                    smiles=selected_smiles, model_choice=model_choice
                ).generate_and_display_similar_smiles(),
                header_text="Select a non-sensitizer SMILES for generating similar smiles",
                spinner_text="Exploring chemical space...",
            )
            st.write("")
            st.write("")
            st.write("")


if __name__ == "__main__":
    main()
