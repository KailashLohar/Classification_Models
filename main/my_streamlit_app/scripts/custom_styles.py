def load_custom_styles():
    return """
    <style>
    .title-box {
        position: sticky;
        top: 0;
        background-color: #90a4ae;
        border-radius: 20px;
        padding: 5px 0;
        margin-top: -50px;
        margin-bottom: 10px;
        text-align: center;
        box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
        border: 0px solid #8b1312;
        width: 90%;
        margin-left: auto;
        margin-right: auto;
        z-index: 1000;
    }
    .title-box h1 {
        color: #000000;
        margin: 0;
        display: inline-block;
    }
    .model-info {
        position: absolute;
        right: 0px;
        top: 110px; /* Position adjusted */
        font-weight: bold;
        color: #000000; /* Black color for general text */
        text-align: right;
    }
    .model-info p {
        margin: 0;
        line-height: 1;
    }
    .model-info a {
        color: #005f8e; /* Blue color for hyperlink */
        text-decoration: none;
    }
    .model-info a:hover {
        text-decoration: underline;
    }
    </style>
    """


# def load_main_title():
#     return """
#     <div class="title-box">
#         <h1>Skin sensitivity prediction</h1>
#         <div class="model-info">
#             <p>The dataset used for this model is from</p>
#             <p><a href="https://pubs.acs.org/doi/10.1021/acs.chemrestox.3c00396" target="_blank">Wang et al. Chem. Res. Toxicol. 2024</a></p>
#         </div>
#     </div>
#     """

def load_main_title():
    return """
    <div class="title-box">
        <h1>Skin sensitivity prediction</h1>
    </div>
    """


def load_sidebar_styles():
    return """
    <style>
    [data-testid="stSidebar"][aria-expanded="true"] > div:first-child {
        background-color: #CBC3E3;
    }
    .sidebar .sidebar-content {
        margin-top: -50px;
    }
    </style>
    """

def load_sidebar_button_styles():
    return """
    <style>
    [data-testid="stSidebar"][aria-expanded="true"] > div:first-child {
        background-color: #CBC3E3;
    }
    .stButton>button {
        background-color: #CBCBCB;
        color: #333333;
        font-weight: bold;
        border: 2px solid #8B1000;
        border-radius: 10px;
        padding: 8px 20px;
    }
    [data-testid="stSidebar"] {
        overflow: hidden;
    }
    </style>
    """

def hide_streamlit_page_links():
    return """
    <style>
    /* Hide default Streamlit page links */
    [data-testid="stSidebarNav"] {
        display: none;
    }
    </style>
    """

def load_centered_table_styles():
    return """
    <style>
    .centered-table {
        width: 100%;
        border-collapse: collapse;
    }
    .centered-table th, .centered-table td {
        text-align: center;
        padding: 8px;
    }
    .centered-table td:first-child {
        text-align: left;
    }
    </style>
    """