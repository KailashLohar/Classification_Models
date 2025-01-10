import os
import base64
import subprocess
import pandas as pd
import streamlit as st

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Draw
from rdkit.Chem.Draw import rdMolDraw2D
from io import BytesIO, StringIO

from scripts.gnn import predict

def add_custom_header_and_footer(header_and_footer_color, logo_image_path, header_background_path, background_image, title, subtitle, more_info_url):
    """Adds custom header with animated subtitle, GitHub link button, and footer to the app."""
    with open(logo_image_path, "rb") as image_file:
        logo_image = base64.b64encode(image_file.read()).decode("utf-8")

    with open(header_background_path, "rb") as image_file:
        header_background = base64.b64encode(image_file.read()).decode()

    with open(background_image, "rb") as image_file:
        background_image = base64.b64encode(image_file.read()).decode()

    st.markdown(
        f"""
        <style>
            .stApp {{
                background: linear-gradient(rgba(255, 255, 255, 0.8), rgba(255, 255, 255, 0.8)), 
                            url(data:image/jpeg;base64,{background_image});
                background-size: cover;
                background-repeat: no-repeat;
                background-attachment: fixed;
                background-blend-mode: lighten; /* Optional: Adjust blend mode */
            }}
        </style>
        """,
        unsafe_allow_html=True
    )

    header_css = f'''
    <style>
        [data-testid="stAppViewContainer"] {{
            padding: 0;
        }}
        [data-testid="stHeader"] {{
            display: none; /* Remove default Streamlit header */
        }}
        @keyframes popIn {{
            0% {{ opacity: 0; transform: scale(0.5); }}
            100% {{ opacity: 1; transform: scale(1); }}
        }}
        
        .custom-header {{
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            background: linear-gradient(to right, rgba(255, 255, 255, 0.95), rgba(255, 255, 255, 0.9), rgba(255, 255, 255, 0.0)),
                        url(data:image/jpeg;base64,{header_background}) no-repeat center;
            background-size: cover;
            box-shadow: 0px 4px 15px rgba(0, 0, 0, 0.2);
            padding: 5px 10px;
            z-index: 1000;
            display: flex;
            align-items: center;
            height: 80px; /* Height for the header */
        }}

        .custom-header .center-content {{
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            flex: 1;
            text-align: center;
        }}

        .custom-header .title {{
            font-family: "Times New Roman", serif;
            font-size: 25px;
            font-weight: bold;
            color: {header_and_footer_color};
            margin: 0;
        }}
        
        .custom-header .subtitle {{
            font-family: "Times New Roman", serif;
            font-size: 16px;
            color: {header_and_footer_color};
            margin: 0;
            display: flex;
            gap: 5px;
        }}

        .custom-header img {{
            height: 40px;
            margin-left: 9px;
        }}

        .custom-footer {{
            position: fixed;
            bottom: 0;
            left: 0;
            width: 100%;
            background: {header_and_footer_color};
            color: white;
            text-align: center;
            padding: 8px 0;
            font-size: 14px;
            z-index: 1000;
        }}
    </style>
    <div class="custom-header">
        <img src="data:image/svg+xml;base64,{logo_image}" alt="Logo" />
        <div class="center-content">
            <div class="title">{title}</div>
            <div class="subtitle">{subtitle}</div>
        </div>
        <div class="moreinfo-button">
            <a href="{more_info_url}" target="_blank"
               style="background-color: {header_and_footer_color}; color: white; text-decoration: none; 
                      border: none; border-radius: 10px; padding: 8px 10px; font-size: 14px;">
                More information
            </a>
        </div>
    </div>
    <div class="custom-footer">
        Confidential and Proprietary. Copyright © 2017-25
    </div>
    '''
    st.markdown(header_css, unsafe_allow_html=True)



def add_custom_css():
    """Add consolidated custom CSS for the application."""
    st.markdown(
        """
        <style>
            /* Sidebar styling */
            [data-testid="stSidebar"] {
                width: 260px !important;
                position: fixed !important;
                left: 15px;
                top: 0 !important;
                height: 80vh !important;
                overflow: auto !important;
                background: #d4d7dc;
                z-index: 10;
                margin-top: 100px;
                color: #31473A;
                border: 0px solid #31473A;
                box-shadow: 0px 4px 15px rgba(0, 0, 0, 0.2);
                backdrop-filter: blur(5px);
                border-radius: 15px;
                padding: 10px;
            }

            [data-testid="stSidebar"] > div:first-child {
                margin-top: -80px !important;
            }

            [data-testid="stAppViewContainer"] {
                margin-left: 260px !important;
                padding: 1rem;
                width: calc(100% - 260px) !important;
            }

            [data-testid="stSidebar"] {
                overflow-x: hidden !important;
            }

            .block-container {
                padding: 2.5rem 2rem;
                max-width: 100% !important;
            }
    
            div[data-testid="stTabs"] hr {
                display: none !important; /* Hide horizontal line if present */
            }
    
            /* Default styling for all tabs */
            div[data-testid="stTabs"] button {
                background-color: #d4d7dc; /* Tab background color */
                color: black; /* Tab text color */
                font-size: 22px; /* Increase font size */
                padding: 10px 15px; /* Add padding for better appearance */
                border: none; /* Remove default border */
                border-radius: 15px 15px 0px 0px; /* Round top corners only */
                transition: background-color 0.3s ease; /* Smooth transition */
            }
    
            /* Styling for hovered tabs */
            div[data-testid="stTabs"] button:hover {
                background-color: #d4d7dc; /* Background color on hover */
                color: black; /* Hover text color */
            }
    
            /* Styling for the active (selected) tab */
            div[data-testid="stTabs"] button[aria-selected="true"] {
                background-color: #6a8f6b !important; /* Active tab background color */
                color: white !important; /* Active tab text color */
            }
    
            /* Ensure container styling remains consistent */
            div[data-testid="stTabs"] {
                margin-top: 12px !important;
                margin-left: 0 !important;
                margin-right: 0 !important;
                width: 100% !important;
            }

            [data-testid='stFileUploader'] {
                width: max-content;
            }
            [data-testid='stFileUploader'] section {
                padding: 0;
                float: left;
            }
            [data-testid='stFileUploader'] section > input + div {
                display: none;
            }
            [data-testid='stFileUploader'] section + div {
                float: right;
                padding-top: 0;
            }

            /* File uploader styling */
            div[data-testid="stFileUploader"] {
                div div {display: none !important;}
                label {color: blue !important;}
                margin-top: -0.5rem !important;
            }
            
            /* Radio button adjustments */
            div[data-baseweb="radio"] > div {
                gap: 0.5rem !important;
            }

            div[data-testid="stRadio"] > label {
                color: blue !important;
                font-size: 1rem !important;
                margin-top: 1rem !important;
                margin-bottom: 0.0rem !important;
            }



            input[type="text"] {
                width: 100% !important;
            }

            /* Button styling */
            div.stButton > button {
                background-color: #dde7dd !important;
                color: black !important;
                padding: 8px 20px !important;
                border-radius: 10px !important;
                font-size: 18px !important;
                border: 0px solid #282d56 !important;
                box-shadow: 0px 4px 6px rgba(0, 0, 0, 0.2) !important;
                cursor: pointer !important;
                transition: all 0.2s ease-in-out;
            }

            div.stButton > button:hover {
                background-color: #6a8f6b !important;
                color: white !important;
                box-shadow: 0px 6px 8px rgba(0, 0, 0, 0.3) !important;
            }

            /* Custom CSS for the download button */
            .download-button {
                background-color: #dde7dd; 
                color: black; 
                text-decoration: none; 
                padding: 10px 20px; 
                border-radius: 5px; 
                font-size: 14px; 
                transition: background-color 0.3s ease;
            }
            .download-button:visited {
                color: black; /* Ensure visited links stay black */
            }
            .download-button:hover {
                background-color: #6a8f6b; 
                color: white; 
                text-decoration: none; /* Prevent underline on hover */
            }
            
        </style>
        """,
        unsafe_allow_html=True,
    )


def initialize_session_state():
    session_keys = ['sensitizer_smiles', 'prediction_table', 'selected_smiles', 'file_uploaded', 
                    'gnn_results', 'counterfactual_generated', 'similar_smiles_generated']
    for key in session_keys:
        if key not in st.session_state:
            st.session_state[key] = [] if 'smiles' in key else None


def handle_input_option():
    input_option = st.sidebar.radio("User Input", ["Upload SMILES", "Enter SMILES", "Use sample"],
                                    help="Please upload a CSV file with SMILES in a column named 'SMILES'.")
    uploaded_file, single_smiles = None, ""
    upload_status = {"csv_file": False, "smiles_input": False, "sample_file": False,}

    if input_option == "Upload SMILES":
        uploaded_file = st.sidebar.file_uploader("Upload CSV File", type=["csv"], label_visibility="collapsed")
        if uploaded_file:
            upload_status["csv_file"] = True
            st.sidebar.markdown('<p style="font-size:14px; color:green;">CSV file uploaded ✅</p>', unsafe_allow_html=True)
    elif input_option == "Enter SMILES":
        single_smiles = st.sidebar.text_input("Enter a SMILES string", key="smiles_input")
        st.markdown('<div class="custom-text-input"></div>', unsafe_allow_html=True)
        if single_smiles:
            upload_status["smiles_input"] = True
            st.sidebar.markdown('<p style="font-size:14px; color:green;">SMILES entered ✅</p>', unsafe_allow_html=True)
    elif input_option == "Use sample":
        uploaded_file = StringIO(open("scripts/test_sensitizer.csv").read())
        upload_status["sample_file"] = True
        st.sidebar.markdown('<p style="font-size:14px; color:green;">Sample file loaded ✅</p>', unsafe_allow_html=True)

    return input_option, uploaded_file, single_smiles



def handle_file_upload(input_option, uploaded_file, single_smiles):
    
    if input_option == "Upload SMILES" and uploaded_file:
        return pd.read_csv(uploaded_file)
    elif input_option == "Enter SMILES" and single_smiles:
        return pd.DataFrame({'SMILES': [single_smiles]})
    elif input_option == "Use sample" and uploaded_file:
        return pd.read_csv(uploaded_file)
    else:
        st.error("Please upload a CSV file, enter a SMILES string, or load a sample.")
        return None


def process_smiles_data(df):
    smiles_col = find_smiles_column(df)
    if smiles_col:
        df = df.rename(columns={smiles_col: 'SMILES'})
        df = df.drop(columns=['Target'], errors='ignore')
    return df

def find_smiles_column(df):
    possible_smiles_columns = ['SMILES', 'SMILE', 'Smiles', 'Chemicals', 'Compounds']
    for column in df.columns:
        if column in possible_smiles_columns:
            return column
    return None


def load_training_data(model_choice):
    training_data_map = {
        "in vitro (H-CLAT)": "skin_hCLAT/training_data_hCLAT.csv",
        "in vitro (KeratinoSens)": "skin_KeratinoSens/training_data_KeratinoSens.csv",
        "in vivo (LLNA)": "skin_LLNA/training_data_LLNA.csv",
        "in chemico (DPRA)": "skin_DPRA/training_data_DPRA.csv",
        "human": "skin_Human/training_data_Human.csv",
    }
    
    training_data_file = training_data_map.get(model_choice)
    if training_data_file:
        training_data = pd.read_csv(training_data_file)
        return training_data['SMILES'].tolist()
    else:
        raise ValueError("Invalid model choice.")


def calculate_highest_similarity(smiles, model_choice):
    training_smiles = load_training_data(model_choice)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0.0

    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    highest_similarity = 0.0
    for train_smiles in training_smiles:
        train_mol = Chem.MolFromSmiles(train_smiles)
        if train_mol is not None:
            train_fp = AllChem.GetMorganFingerprintAsBitVect(train_mol, 2, nBits=1024)
            similarity = DataStructs.TanimotoSimilarity(fp, train_fp)
            highest_similarity = max(highest_similarity, similarity)

    return highest_similarity * 100


def generate_structure_image(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    img = Draw.MolToImage(mol, size=(200, 100))
    buffered = BytesIO()
    img.save(buffered, format="PNG")
    img_str = base64.b64encode(buffered.getvalue()).decode()
    return f'<img src="data:image/png;base64,{img_str}"/>'


def add_chemical_structure_column(results):
    results['Chemical structure'] = results['SMILES'].apply(generate_structure_image)
    return results

    
def style_table(filtered_df):
    styled_df = filtered_df.style \
                           .set_table_styles([
                               {'selector': 'thead th', 
                                'props': [('background-color', '#6a8f6b'),  
                                          ('color', '#FFFFFF'), 
                                          ('font-size', '14px'), 
                                          ('text-align', 'center'),
                                          ('font-weight', '100')]} 
                           ]) \
                           .apply(lambda row: ['background-color: #FFFFFF; color: black;'  # Light blue for all rows
                                                for _ in row], axis=1) \
                           .set_properties(**{'font-size': '14px', 'background-color': '#FFFFFF', 'color': 'black'}) 
        
    return styled_df

# def update_session_state_with_results(df, model_choice):
#     st.session_state.gnn_results = predict(df, model_choice)
#     results = st.session_state.gnn_results

#     st.session_state.sensitizer_smiles = results[results['Prediction'] == 1]['SMILES'].tolist()
#     st.session_state.non_sensitizer_smiles = results[results['Prediction'] == 0]['SMILES'].tolist()

#     results['Prediction'] = results['Prediction'].replace({0: 'Non-sensitizer', 1: 'Sensitizer'})
#     results['Prediction'] = results['Prediction'].apply(
#         lambda x: f"<span style='color: {'red' if x == 'Sensitizer' else 'green'};'>{x}</span>")
#     results['Confidence (%)'] = results['Confidence (%)'].apply(lambda x: f"{x:.0f}")
#     results['Applicability'] = results['SMILES'].apply(lambda x: calculate_highest_similarity(x, model_choice))
#     results['Applicability'] = results['Applicability'].apply(lambda x: f"{x:.0f}%")
#     results = add_chemical_structure_column(results)

#     columns_to_display = ['SMILES', 'Chemical structure', 'Prediction', 'Applicability']
#     if 'Ground Truth' in df.columns:
#         results['Ground Truth'] = df['Ground Truth']
#         columns_to_display.insert(2, 'Ground Truth')

#     results.index += 1

#     # Convert DataFrame to HTML
#     styled_df = style_table(results[columns_to_display])
#     if 'Ground Truth' in styled_df.columns:
#         styled_df_html = styled_df.to_html(index=False, escape=False)
#         scrollable_table_html = (
#             f"<div style='display: flex; justify-content: center; margin-top: 0px; margin-bottom: 0px;'>"
#             f"<div style='max-height: 500px; overflow-y: scroll; width: auto; border: 0px solid #282d56; margin: 0; padding: 0; box-sizing: border-box;'>"
#             f"<table style='border-collapse: collapse; width: 100%; table-layout: fixed; background-color: #FFFFFF; margin: 0; padding: 0;'>"
#             f"<style>"
#             f"  .col1 {{ text-align: center; }}"
#             f"  .col2 {{ text-align: center; }}"
#             f"  .col3 {{ text-align: center; }}"
#             f"  .col4 {{ text-align: center; }}"
#             f"  table {{ margin: 0; padding: 0; }}"
#             f"  tbody {{ margin: 0; padding: 0; box-sizing: border-box; }}"
#             f"  tr, td {{ margin: 0; padding: 0; border: none; }}"
#             f"</style>"
#             f"<colgroup>"
#             f"  <col style='width: 35px;'> <!-- Fix index column width -->"
#             f"  <col style='width: 300px;'> <!-- Fix SMILES column width -->"
#             f"  <col style='width: 200px;'> <!-- Fix Chemical structure column width -->"
#             f"  <col style='width: 100px;'> <!-- Fix Ground Truth column width -->"
#             f"  <col style='width: 100px;'> <!-- Fix Prediction column width -->"
#             f"  <col style='width: 80px;'> <!-- Fix Applicability column width -->"
#             f"</colgroup>"
#             f"<thead style='position: sticky; top: 0; background-color: #283456; color: white; z-index: 1;'>"
#             f"{styled_df_html.split('<thead>', 1)[1].split('</thead>', 1)[0]}"
#             f"</thead>"
#             f"<tbody style='background-color: #FFFFFF; word-wrap: break-word; margin: 0; padding: 0; box-sizing: border-box;'>"
#             f"{styled_df_html.split('<tbody>', 1)[1]}</tbody>"
#             f"</table>"
#             f"</div>"
#             f"</div>"
#         )
#         st.session_state.prediction_table = scrollable_table_html

#     else:
#         styled_df_html = styled_df.to_html(index=False, escape=False)
#         scrollable_table_html = (
#             f"<div style='display: flex; justify-content: center; margin-top: 0px; margin-bottom: 0px;'>"
#             f"<div style='max-height: 500px; overflow-y: scroll; width: auto; border: 0px solid #282d56; margin: 0; padding: 0; box-sizing: border-box;'>"
#             f"<table style='border-collapse: collapse; width: 100%; table-layout: fixed; background-color: #FFFFFF; margin: 0; padding: 0;'>"
#             f"<style>"
#             f"  .col1 {{ text-align: center; }}"
#             f"  .col2 {{ text-align: center; }}"
#             f"  .col3 {{ text-align: center; }}"
#             f"  table {{ margin: 0; padding: 0; }}"
#             f"  tbody {{ margin: 0; padding: 0; box-sizing: border-box; }}"
#             f"  tr, td {{ margin: 0; padding: 0; border: none; }}"
#             f"</style>"
#             f"<colgroup>"
#             f"  <col style='width: 35px;'> <!-- Fix index column width -->"
#             f"  <col style='width: 300px;'> <!-- Fix SMILES column width -->"
#             f"  <col style='width: 200px;'> <!-- Fix Chemical structure column width -->"
#             f"  <col style='width: 100px;'> <!-- Fix Prediction column width -->"
#             f"  <col style='width: 80px;'> <!-- Fix Applicability column width -->"
#             f"</colgroup>"
#             f"<thead style='position: sticky; top: 0; background-color: #283456; color: white; z-index: 1;'>"
#             f"{styled_df_html.split('<thead>', 1)[1].split('</thead>', 1)[0]}"
#             f"</thead>"
#             f"<tbody style='background-color: #FFFFFF; word-wrap: break-word; margin: 0; padding: 0; box-sizing: border-box;'>"
#             f"{styled_df_html.split('<tbody>', 1)[1]}</tbody>"
#             f"</table>"
#             f"</div>"
#             f"</div>"
#         )
#         st.session_state.prediction_table = scrollable_table_html

    
def generate_table_html(results, columns_to_display):
    """
    Generate HTML for a sortable table with customizable columns and specified widths.
    """
    column_widths = {
        "Index": "35px",
        "SMILES": "300px",
        "Chemical structure": "200px",
        "Ground Truth": "100px",
        "Prediction": "100px",
        "Applicability": "80px",
    }

    colgroup_html = (
        f'<col style="width: {column_widths.get("Index", "100px")};">'
        + ''.join(
            f'<col style="width: {column_widths.get(col, "100px")};">' for col in columns_to_display
        )
    )

    table_html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/tablesort/5.2.1/tablesort.min.js"></script>
        <style>
            th {{
                cursor: pointer;
                background-color: #283456;
                color: white;
                padding: 8px;
                position: sticky;
                top: 0;
                z-index: 2;
            }}
            td {{
                padding: 8px;
                text-align: center;
                word-wrap: break-word;
                overflow-wrap: break-word;
            }}
            table {{
                border-collapse: collapse;
                width: 100%;
                table-layout: fixed;
                background-color: #FFFFFF;
            }}
            th, td {{
                border: 1px solid #ddd;
            }}
            td:first-child {{
                position: sticky;
                left: 0;
                background-color: #f4f4f4;
                z-index: 1;
            }}
            .table-container {{
                max-height: 500px;
                overflow-y: auto;
            }}
        </style>
    </head>
    <body>
        <div class="table-container">
            <table id="sortable-table">
                <colgroup>
                    {colgroup_html}
                </colgroup>
                <thead>
                    <tr>
                        <th>#</th>
                        {''.join(f'<th>{col}</th>' for col in columns_to_display)}
                    </tr>
                </thead>
                <tbody>
                    {''.join(
                        f'<tr><td>{i}</td>' + ''.join(f'<td>{val}</td>' for val in row) + '</tr>'
                        for i, row in enumerate(results[columns_to_display].values, start=1)
                    )}
                </tbody>
            </table>
        </div>
        <script>
            new Tablesort(document.getElementById('sortable-table'));
        </script>
    </body>
    </html>
    """
    return table_html


def update_session_state_with_results(df, model_choice):
    st.session_state.gnn_results = predict(df, model_choice)
    results = st.session_state.gnn_results

    st.session_state.sensitizer_smiles = results[results['Prediction'] == 1]['SMILES'].tolist()
    st.session_state.non_sensitizer_smiles = results[results['Prediction'] == 0]['SMILES'].tolist()

    results['Prediction'] = results['Prediction'].replace({0: 'Non-sensitizer', 1: 'Sensitizer'})
    results['Prediction'] = results['Prediction'].apply(
        lambda x: f"<span style='color: {'red' if x == 'Sensitizer' else 'green'};'>{x}</span>"
    )
    results['Confidence (%)'] = results['Confidence (%)'].apply(lambda x: f"{x:.0f}")
    results['Applicability'] = results['SMILES'].apply(lambda x: calculate_highest_similarity(x, model_choice))
    results['Applicability'] = results['Applicability'].apply(lambda x: f"{x:.0f}%")
    results = add_chemical_structure_column(results)

    columns_to_display = ['SMILES', 'Chemical structure', 'Prediction', 'Applicability']
    if 'Ground Truth' in df.columns:
        results['Ground Truth'] = df['Ground Truth']
        columns_to_display.insert(2, 'Ground Truth')

    results.index += 1

    table_html = generate_table_html(results, columns_to_display)
    st.session_state.prediction_table = table_html


def render_tab(tab_type, smiles_list, model_choice, button_label, button_key, action_function, header_text, spinner_text):
    if smiles_list:
        st.markdown(f"<div style='font-size: 17px; color: #8e5572; font-weight: bold; margin-bottom: 0px;'>{header_text}</div>", unsafe_allow_html=True)
        selected_smiles = st.selectbox("", smiles_list, index=0, key=f"select_{tab_type.lower()}_smiles")
    
        if st.button(button_label, key=button_key):
            with st.spinner(spinner_text):
                action_function(selected_smiles)
            if tab_type == "Sensitizer":
                st.session_state.counterfactual_generated = True
            else:
                st.session_state.similar_smiles_generated = True
    else:
        st.write(f"No {tab_type.lower()} SMILES found.")
