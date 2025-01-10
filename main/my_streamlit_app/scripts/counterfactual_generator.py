import os
import base64
import pandas as pd
import streamlit as st
import exmol
from scripts.gnn import predict  

class CounterfactualGenerator:
    def __init__(self, smiles, model_choice):
        self.smiles = smiles
        self.model_choice = model_choice
    
    def generate_and_display_counterfactuals(self):
        predictions = predict(pd.DataFrame({'SMILES': [self.smiles]}), self.model_choice)

        if predictions['Prediction'].iloc[0] == 1:
            samples = exmol.sample_space(self.smiles, self.exmol_model, batched=True)
            cfs = exmol.cf_explain(samples)

            # Display counterfactual SMILES
            exmol.plot_cf(cfs)
            svg_cf = exmol.insert_svg(cfs, mol_fontsize=16)
            svg_cf_base64 = base64.b64encode(svg_cf.encode('utf-8')).decode('utf-8')
            
            # Add note
            st.markdown("""
                <div style="text-align: left; font-size:16px; color:black; margin-bottom: -5px; margin-left: 130px;">
                    Note: Here f(x) = 1 represents Sensitizer and f(x) = 0 represents Non-sensitizer.
                </div>
            """, unsafe_allow_html=True)

            # Display centered SVG

            st.markdown("""
                <style>
                    .svg-container {
                        border: 0px solid #6a8f6b; /* Border color */
                        border-radius: 10px; /* Rounded corners */
                        box-shadow: 0px 4px 6px rgba(0, 0, 0, 0.2) !important;
                        padding: 10px; /* Padding inside the container */
                        background-color: #ffffff; /* Background color inside border */
                        box-sizing: content-box;
                        text-align: center; /* Center align contents */
                        max-width: 800px; /* Add max-width to prevent full width */
                        margin-left: 130px; /* Shift images to the right by 100px */
                    }
                </style>
            """, unsafe_allow_html=True)
            
            st.markdown(f"""
                <div class="svg-container">
                    <div style='margin-top: 10px;'>
                        <img src="data:image/svg+xml;base64,{svg_cf_base64}" width="800px"/>
                    </div>
                </div>
            """, unsafe_allow_html=True)


            # Display table
            st.markdown("<h5 style='color:#8e5572; font-size:16px; margin-top:20px; text-align:center;'>Counterfactual SMILES with prediction</h5>", unsafe_allow_html=True)

            cf_smiles = [cf.smiles for cf in cfs]
            similarity_scores = [f"{cf.similarity:.2f}" for cf in cfs]
            cf_predictions = self.exmol_model(cf_smiles)
            cf_predictions_labels = ["Sensitizer" if pred == 1 else "Non-sensitizer" for pred in cf_predictions]
            
            cf_df = pd.DataFrame({'SMILES': cf_smiles, 'Similarity score': similarity_scores, 'Prediction': cf_predictions_labels})
            cf_df = cf_df.iloc[1:].reset_index(drop=True)
            cf_df.index = range(1, len(cf_df) + 1)

            styled_df = self.style_table(cf_df)
            styled_df_html = styled_df.to_html(index=True, escape=False)
            centered_html = f"""<div style='display: flex; justify-content: center; margin-top: 0px;'> {styled_df_html}"""
            st.markdown(centered_html, unsafe_allow_html=True)
        
            # Display download button
            csv_data = cf_df.to_csv(index=False).replace('\n', '%0A').replace(',', '%2C')
            st.markdown(f"""
                <div style="text-align: center; margin-top: 5px;">
                    <a href="data:text/csv;charset=utf-8,{csv_data}" download="counterfactuals_smiles.csv" class="download-button">
                        Download (.csv)
                    </a>
                </div>
            """, unsafe_allow_html=True)

        else:
            st.markdown("<span style='color:#0000FF; font-weight:bold;'>No counterfactuals generated as the input is not classified as a sensitizer.</span>", unsafe_allow_html=True)

        return predictions['Prediction'].iloc[0]

    def exmol_model(self, smiles_list):
        if isinstance(smiles_list, str):  
            smiles_list = [smiles_list]

        df = pd.DataFrame({'SMILES': smiles_list})
        predictions = predict(df, self.model_choice)
        return predictions['Prediction'].tolist()

    @staticmethod
    def style_table(cf_df):
        styled_df = cf_df.style \
                           .set_table_styles([{'selector': 'thead th', 
                                               'props': [('background-color', '#282d56'), 
                                                         ('color', 'white'), 
                                                         ('font-size', '14px'), 
                                                         ('text-align', 'center')]}]) \
                           .apply(lambda row: ['background-color: #EDF4FB;' 
                                                if row.name % 2 == 0 
                                                else 'background-color: #EDF4FB; color: black;' for _ in row], axis=1) \
                           .set_properties(**{'font-size': '14px', 'text-align': 'center'})
        
        return styled_df
