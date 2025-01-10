import os
import subprocess
import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import AllChem, DataStructs
from scripts.gnn import predict

class ChemicalSpaceExplorer:
    def __init__(self, smiles, model_choice):
        self.smiles = smiles
        self.model_choice = model_choice
    
    def generate_and_display_similar_smiles(self):
        examples_dir = "exahustive_search/examples"
        os.makedirs(examples_dir, exist_ok=True)
        
        smiles_path = os.path.join(examples_dir, "smiles.txt")
        with open(smiles_path, 'w') as f:
            f.write(self.smiles)
        
        command = [
            "python", "exahustive_search/generate_smiles.py",
            "--model", "exahustive_search/paper_checkpoints/ecfp4_with_counts_with_rank",
            "--input-smiles", smiles_path,
            "--samples", "100",
            "--result-path", "exahustive_search/samples.csv"
        ]
        
        subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

        samples_path = "exahustive_search/samples.csv"
        if not os.path.exists(samples_path):
            st.write("samples.csv file not found.")
            return

        df = pd.read_csv(samples_path)
        df['tanimoto'] = df['tanimoto'].round(2)
        df = df[(df['tanimoto'] != 1) & (df['is_valid'] == 1)]
        df = df.sort_values(by='tanimoto', ascending=False)

        generated_smiles_df = pd.DataFrame({'SMILES': df['generated_smiles']})
        
        # Use the predict function from gnn module
        predictions = predict(generated_smiles_df, self.model_choice)

        if 'SMILES' not in predictions.columns:
            st.error("'SMILES' column not found in predictions DataFrame!")
            return

        non_sensitizers_df = generated_smiles_df[predictions['Prediction'] == 0]
        df = df[df['generated_smiles'].isin(non_sensitizers_df['SMILES'])]
        
        selected_rows = []
        last_value = None
        for _, row in df.iterrows():
            if last_value is None or row['tanimoto'] <= last_value - 0.02:
                selected_rows.append(row)
                last_value = row['tanimoto']
        
        filtered_df = pd.DataFrame(selected_rows).head(6)
        filtered_df = filtered_df[['generated_smiles', 'tanimoto']]
        filtered_df.columns = ['SMILES', 'Similarity score']
        filtered_df['Similarity score'] = filtered_df['Similarity score'].map(lambda x: f"{x:.2f}")
        filtered_df = pd.merge(filtered_df, predictions[['SMILES', 'Prediction']], how='inner', on='SMILES')
        filtered_df['Prediction'] = filtered_df['Prediction'].apply(lambda x: "Non-sensitizer" if x == 0 else "Sensitizer")
        
        self.display_svg_images(filtered_df)

        st.markdown("<h5 style='color:#8e5572; font-size:16px; margin-top:50px; text-align:center;'>SMILES with similarity scores and predictions</h5>", unsafe_allow_html=True)
        filtered_df.index = range(1, len(filtered_df) + 1)
        
        styled_df = self.style_table(filtered_df)
        styled_df_html = styled_df.to_html(index=True, escape=False)
        centered_html = f"""<div style='display: flex; justify-content: center; margin-top: 0px;'> {styled_df_html}"""
        st.markdown(centered_html, unsafe_allow_html=True)

        csv_data = filtered_df.to_csv(index=False).replace('\n', '%0A').replace(',', '%2C')
        st.markdown(f"""
            <div style="text-align: center; margin-top: 5px;">
                <a href="data:text/csv;charset=utf-8,{csv_data}" download="similar_smiles.csv" class="download-button">
                    Download (.csv)
                </a>
            </div>
        """, unsafe_allow_html=True)


    def display_svg_images(self, filtered_df):
        smiles_list = filtered_df['SMILES'].tolist()
        similarity_scores = filtered_df['Similarity score'].tolist()
        
        svg_images = []
        for smile in smiles_list:
            mol = Chem.MolFromSmiles(smile)
            if mol is not None:
                # Reduce the image size by changing the dimensions
                drawer = rdMolDraw2D.MolDraw2DSVG(280, 280)  
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText().replace("svg:", "")
                svg_images.append(svg)
        
        # Add CSS to create borders around the SVG with proper padding
        st.markdown("""
            <style>
                .svg-columns-container {
                    display: flex;
                    justify-content: center; /* Center align columns */
                    margin-top: 10px; /* Add spacing above the section */
                    margin-bottom: 10px; /* Add spacing below */
                    gap: 10px; /* Reduce gap between columns */
                }
        
                .svg-container {
                    border: 0px solid #6a8f6b; /* Border color */
                    border-radius: 10px; /* Rounded corners */
                    box-shadow: 0px 4px 6px rgba(0, 0, 0, 0.2) !important;
                    padding: 0px; /* Padding inside the container */
                    background-color: #ffffff; /* Background color inside border */
                    box-sizing: content-box;
                    text-align: center; /* Center align contents */
                    max-width: 300px; /* Add max-width to prevent full width */
                    margin-left: 100px; /* Shift images to the right by 100px */
                    display: flex;
                    flex-direction: column;
                    justify-content: center;
                    align-items: center;
                }
            </style>
        """, unsafe_allow_html=True)
        
        # Wrap SVG images in a container with the custom class
        for i in range(0, len(svg_images), 2):
            st.markdown('<div class="svg-columns-container">', unsafe_allow_html=True)
            cols = st.columns(2)
            for col, svg, score in zip(cols, svg_images[i:i+2], similarity_scores[i:i+2]):
                with col:
                    st.markdown(f"""
                        <div class="svg-container">
                            <div style='margin: 0px;'>{svg}</div>
                            <p style='font-size:15px; margin: 0px 0px 10px;'>Similarity score: 
                                <strong><span style='color:green;'>{float(score):.2f}</span></strong>
                            </p>
                        </div>
                    """, unsafe_allow_html=True)
            st.markdown('</div>', unsafe_allow_html=True)


    @staticmethod
    def style_table(filtered_df):
        styled_df = filtered_df.style \
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
