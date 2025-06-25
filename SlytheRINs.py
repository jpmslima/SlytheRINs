# SlytheRINs: Residue Interaction Network Analysis Dashboard
# Developed by the EvoMol-Lab (github.com/evomol-lab) - BioME, UFRN, Brazil
# Use RING generated data from ring.biocomputingup.it

# Import necessary libraries
import streamlit as st
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
from pandas.errors import EmptyDataError
import os
import scipy.stats as st_sci # Renamed to avoid conflict with streamlit 'st'
import numpy as np
import io # Required for BytesIO for image download
from collections import Counter # For degree distribution
import requests # For downloading AlphaMissense data
import re # For regex operations (extracting UniProt ID from filename)
from fpdf import FPDF # For PDF generation
from datetime import datetime # For PDF timestamp
import zipfile # For creating zip archives

#Keep the original colors on png download
import plotly.io as pio
pio.templates.default = "plotly_white"

# --- Helper function for downloads ---
@st.cache_data # Use Streamlit's caching for efficiency
def convert_df_to_tsv(df):
    """Converts a DataFrame to a TSV string for download."""
    return df.to_csv(sep='\t', index=False).encode('utf-8')

@st.cache_data
def convert_fig_to_png(fig, scale=2): # Added scale parameter
    """Converts a Plotly figure to PNG bytes for download."""
    try:
        img_bytes = fig.to_image(format="png", scale=scale) 
        return img_bytes
    except Exception as e:
        #Error message without kaleido
        st.error(f"Could not generate PNG. An error ocurred during image conversion. Error: {e}")
        return io.BytesIO() #empty bytes to prevent crashing

# --- PDF Report Generation Class ---
class PDF(FPDF):
    def __init__(self, orientation='P', unit='mm', format='A4'):
        super().__init__(orientation, unit, format)
        self.set_auto_page_break(auto=True, margin=15)
        self.alias_nb_pages() 
        self.slytherins_logo_path = 'SlytheRINs-logo2.png' 
        self.font_name = 'Helvetica' # Default fallback font
        
        # IMPORTANT: Place DejaVuSans font files in a 'fonts' subdirectory
        # relative to your script, or update these paths.
        self.font_paths = {
            'Regular': 'fonts/DejaVuSans.ttf',
            'Bold': 'fonts/DejaVuSans-Bold.ttf',
            'Italic': 'fonts/DejaVuSans-Oblique.ttf', # Or DejaVuSans-Italic.ttf
            'BoldItalic': 'fonts/DejaVuSans-BoldOblique.ttf' # Or DejaVuSans-BoldItalic.ttf
        }
        
        font_loaded_successfully = True
        font_files_missing_messages = []

        try:
            if os.path.exists(self.font_paths['Regular']):
                self.add_font('DejaVu', '', self.font_paths['Regular'], uni=True)
            else:
                font_loaded_successfully = False
                font_files_missing_messages.append(f"Regular: '{self.font_paths['Regular']}' not found.")

            if os.path.exists(self.font_paths['Bold']):
                self.add_font('DejaVu', 'B', self.font_paths['Bold'], uni=True)
            else:
                font_files_missing_messages.append(f"Bold: '{self.font_paths['Bold']}' not found. Bold style may not work as expected.")

            if os.path.exists(self.font_paths['Italic']):
                self.add_font('DejaVu', 'I', self.font_paths['Italic'], uni=True)
            else:
                font_files_missing_messages.append(f"Italic: '{self.font_paths['Italic']}' not found. Italic style may not work as expected.")
            
            if os.path.exists(self.font_paths['BoldItalic']):
                self.add_font('DejaVu', 'BI', self.font_paths['BoldItalic'], uni=True)
            else:
                font_files_missing_messages.append(f"BoldItalic: '{self.font_paths['BoldItalic']}' not found. BoldItalic style may not work as expected.")

            if font_loaded_successfully and os.path.exists(self.font_paths['Regular']): 
                self.font_name = 'DejaVu' 
            elif font_files_missing_messages: 
                 missing_files_str = "\n".join(font_files_missing_messages)
                 st.warning(f"PDF Font Warning: One or more DejaVu font variant files were not found in the 'fonts/' directory. Using default font or available variants. Special characters or styles might not render correctly in the PDF. Missing:\n{missing_files_str}")
                 self.font_name = 'Helvetica' 
            
        except RuntimeError as e: 
            print(f"FPDF Runtime Error during font loading. Using '{self.font_name}'. Error: {e}")
            st.warning(f"PDF Font Warning: A runtime error occurred while loading custom fonts. Using default font. Special characters or styles might not render correctly. Error: {e}")
            self.font_name = 'Helvetica' 
        except Exception as e: 
            print(f"FPDF General Font Loading Error: {e}. Using '{self.font_name}'.")
            st.warning(f"PDF Font Warning: An unexpected error occurred while loading custom fonts. Using default font.")
            self.font_name = 'Helvetica' 
        
        self.set_font(self.font_name, '', 10) 

    def _standardize_text(self, text):
        text = str(text) 
        replacements = { "’": "'", "‘": "'", "“": '"', "”": '"', "–": "-", "—": "--"}
        for unicode_char, ascii_char in replacements.items():
            text = text.replace(unicode_char, ascii_char)
        return text

    def header(self):
        current_font_family = self.font_family; current_font_style = self.font_style; current_font_size = self.font_size_pt
        self.set_font(self.font_name, 'B', 10 if self.page_no() == 1 else 8)
        if self.page_no() == 1:
            try:
                if os.path.exists(self.slytherins_logo_path): self.image(self.slytherins_logo_path, x=10, y=8, w=30); self.ln(5) 
            except Exception as e: print(f"PDF Error: Could not embed logo. {e}") 
            self.set_font(self.font_name, 'B', 18)
            self.cell(0, 10, self._standardize_text('SlytheRINs Analysis Report'), 0, 1, 'C')
            self.set_font(self.font_name, '', 10)
            description_text = ('SlytheRINs: Residue Interaction Network Analysis Dashboard.\n'
                                'Developed by the EvoMol-Lab (github.com/evomol-lab).\n'
                                'BioME, UFRN, Brazil (bioinfo.imd.ufrn.br).\n'
                                'RING data from ring.biocomputingup.it')
            self.multi_cell(0, 5, self._standardize_text(description_text), 0, 'C'); self.ln(5)
            self.set_font(self.font_name, 'I', 9)
            self.cell(0, 8, self._standardize_text(f'Report Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'), 0, 1, 'C'); self.ln(10)
        else: 
            self.set_font(self.font_name, 'I', 8)
            self.cell(0, 10, self._standardize_text('SlytheRINs Report'), 0, 0, 'L')
            self.cell(0, 10, f'Page {self.page_no()}/{{nb}}', 0, 0, 'R'); self.ln(10) 
        self.set_font(current_font_family, style=current_font_style, size=current_font_size)

    def footer(self):
        current_font_family = self.font_family; current_font_style = self.font_style; current_font_size = self.font_size_pt
        self.set_y(-15); self.set_font(self.font_name, 'I', 8)
        self.cell(0, 10, f'Page {self.page_no()}/{{nb}}', 0, 0, 'C')
        self.set_font(current_font_family, style=current_font_style, size=current_font_size)

    def chapter_title(self, title):
        self.set_font(self.font_name, 'B', 14); self.ln(10) 
        self.cell(0, 10, self._standardize_text(title), 0, 1, 'L'); self.ln(4) 

    def section_title(self, title): 
        self.set_font(self.font_name, 'B', 12); self.ln(5)
        self.cell(0, 8, self._standardize_text(title), 0, 1, 'L')
        
    def body_text(self, text): 
        self.set_font(self.font_name, '', 9) 
        self.multi_cell(0, 5, self._standardize_text(text)); self.ln(2)

    def add_metric(self, label, value):
        self.set_font(self.font_name, 'B', 10)
        self.cell(60, 7, self._standardize_text(label) + ":", ln=0) 
        self.set_font(self.font_name, '', 10)
        self.multi_cell(0, 7, self._standardize_text(value), ln=1); self.ln(1)

    def add_plotly_figure_to_pdf(self, fig, title, caption=None, fig_width_mm=170):
        self.section_title(title) 
        if caption: self.body_text(caption) 
        try:
            img_bytes = fig.to_image(format="png", scale=1.5) 
            img_file = io.BytesIO(img_bytes)
            page_width = self.w - 2 * self.l_margin 
            img_render_width = min(fig_width_mm, page_width) 
            x_pos = (self.w - img_render_width) / 2
            self.image(img_file, x=x_pos, w=img_render_width); self.ln(5) 
        except Exception as e:
            #Removing kaleido
            self.set_font(self.font_name, 'I', 9); self.set_text_color(255, 0, 0) 
            self.multi_cell(0, 5, self._standardize_text(f"Error embedding plot '{title}' . Image conversion failed in PDF: {e}"))
            self.set_text_color(0, 0, 0); self.ln()
            
# --- Function to generate PDF ---
def generate_pdf_report(report_elements_list, am_source_msg):
    pdf = PDF()
    pdf.add_page()
    pdf.current_part = 0 

    for element in report_elements_list:
        part = element.get("part")
        if pdf.current_part != part: 
            if part == 1: pdf.chapter_title("Part 1: RIN Comparison & Network Complexity")
            elif part == 2: pdf.add_page(); pdf.chapter_title("Part 2: Chemical Interaction Analysis")
            elif part == 3: 
                pdf.add_page(); pdf.chapter_title("Part 3: AlphaMissense Integration")
                if am_source_msg: pdf.body_text(f"Data Source: {am_source_msg}")
            pdf.current_part = part
            
        if element["type"] == "metric":
            pdf.add_metric(element["label"], element["value"])
        elif element["type"] == "plot" or element["type"] == "am_plot" or \
             element["type"] == "interaction_plot" or element["type"] == "complexity_plot" or \
             element["type"] == "hub_table_plot": 
            pdf.add_plotly_figure_to_pdf(element["fig"], element["title"], element.get("caption"))
        elif element["type"] == "text_summary":
            pdf.section_title(element["title"])
            pdf.body_text(element["content"])
        elif element["type"] == "complexity_metrics":
            pdf.section_title(f"Complexity Metrics for: {element['graph_name']}")
            for k, v in element["metrics"].items():
                pdf.add_metric(k, v)
        elif element["type"] == "hub_table_text": 
            pdf.section_title(element["title"])
            pdf.set_font(pdf.font_name, '', 8) 
            pdf.multi_cell(0,5, element["df_string"])
            pdf.ln()
                
    pdf_output_object = pdf.output(dest='S')
    if isinstance(pdf_output_object, str): return pdf_output_object.encode('latin-1') 
    elif isinstance(pdf_output_object, bytearray): return bytes(pdf_output_object)
    elif isinstance(pdf_output_object, bytes): return pdf_output_object
    else: st.error("Unexpected PDF output type."); return b""


# --- Function to generate ZIP archive ---
def generate_zip_archive(report_elements_list, am_source_msg_zip):
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "a", zipfile.ZIP_DEFLATED, False) as zip_file:
        readme_content = f"SlytheRINs Analysis Report - Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        if am_source_msg_zip:
            readme_content += f"AlphaMissense Data Source: {am_source_msg_zip}\n\n"
        
        data_summary_for_readme = []

        for element in report_elements_list:
            part = element.get("part")
            el_type = element.get("type")
            raw_title = element.get("title", "Untitled_Element")
            sanitized_title = re.sub(r'[^\w\s-]', '', raw_title).strip().replace(' ', '_').replace('-', '_')
            
            if el_type == "df_text": 
                df_to_zip = element["df"]
                filename = f"data/part{part}_{sanitized_title}.tsv"
                try:
                    zip_file.writestr(filename, convert_df_to_tsv(df_to_zip))
                    data_summary_for_readme.append(f"- {filename}")
                except Exception as e:
                    print(f"Error writing df {sanitized_title} to zip: {e}") 
            elif el_type in ["plot", "am_plot", "interaction_plot", "complexity_plot", "hub_table_plot"]: 
                fig_to_zip = element["fig"]
                try:
                    png_bytes = convert_fig_to_png(fig_to_zip, scale=2) 
                    if png_bytes and len(png_bytes) > 0 : 
                        filename = f"plots/part{part}_{sanitized_title}.png"
                        zip_file.writestr(filename, png_bytes)
                        data_summary_for_readme.append(f"- {filename}")
                except Exception as e:
                     print(f"Error converting fig {sanitized_title} to png for zip: {e}") 
            elif el_type == "metric":
                 readme_content += f"Metric (Part {part}): {element['label']} = {element['value']}\n"
            elif el_type == "text_summary":
                 readme_content += f"Summary (Part {part} - {raw_title}):\n{element['content']}\n\n"
            elif el_type == "complexity_metrics":
                 readme_content += f"Complexity Metrics (Part {part} - {element['graph_name']}):\n"
                 for k,v in element['metrics'].items():
                     readme_content += f"  {k}: {v}\n"
                 readme_content += "\n"

        readme_content += "\nIncluded Data Files:\n" + "\n".join(data_summary_for_readme)
        zip_file.writestr("Report_Summary.txt", readme_content)
    return zip_buffer.getvalue()


# --- Page Configuration & Sidebar Setup ---
st.set_page_config(layout="wide") 
with st.sidebar:
    try:
        st.image('SlytheRINs-logo2.png', use_container_width=True)
    except FileNotFoundError: st.warning("SlytheRINs-logo2.png not found.")
    except Exception as e: st.warning(f"Could not load logo: {e}")

    st.title('**SlytheRINs**')
    st.write('*SlytheRINs is an application to analyze and compare Residue Interaction Networks (RINs) calculated from different protein structures and conformations.*')
    st.caption('Developed by the [EvoMol-Lab](https://github.com/evomol-lab) - BioME, UFRN, Brazil')

    st.header("Data Upload")
    st.subheader("1. RIN Data")
    uploaded_files = st.file_uploader("Upload [RING](https://ring.biocomputingup.it/) generated .edges.txt files", accept_multiple_files=True, type="txt")
    use_example_files = st.checkbox("...or use example RIN files")

    st.subheader("2. AlphaMissense Data")
    st.info("The download of AlphaMissense data from the [AlphaFold DB](https://alphafold.ebi.ac.uk/) only works for human proteins.")
    uniprot_id_input = st.text_input("Enter UniProt Accession (e.g., P03891)", key="uniprot_id_input")
    download_am_button = st.button("Download & Process from AlphaFold DB", key="download_am_button")
    
    uploaded_file_am = st.file_uploader("...or upload AlphaMissense .tsv/.csv file", type=["tsv", "csv"], key="am_uploader_sidebar")
    use_example_file_am = st.checkbox("...or use an example AlphaMissense file", key="am_example_checkbox_sidebar")

    st.header("Report Generation")
    if st.button("Prepare Data for Reports", key="prepare_report_button"):
        if 'report_elements' in st.session_state and st.session_state.report_elements:
            st.session_state.reports_ready = True 
            st.success("Report data prepared. Download options below.")
            st.session_state.trigger_pdf_regeneration = True 
            st.session_state.trigger_zip_regeneration = True
        else:
            st.warning("No analysis data available to generate reports. Please process some data first.")

    if st.session_state.get('reports_ready', False):
        if 'pdf_bytes' not in st.session_state or st.session_state.get('trigger_pdf_regeneration', False):
             with st.spinner("Generating PDF report..."):
                st.session_state.pdf_bytes = generate_pdf_report(st.session_state.get('report_elements', []), st.session_state.get('am_data_source_message_for_report', 'N/A'))
                st.session_state.trigger_pdf_regeneration = False 
        if 'pdf_bytes' in st.session_state and st.session_state.pdf_bytes:
            st.download_button(label="Download PDF Report", data=st.session_state.pdf_bytes, file_name="SlytheRINs_Report.pdf", mime="application/pdf", key="download_pdf_report_button")
        
        if 'zip_bytes' not in st.session_state or st.session_state.get('trigger_zip_regeneration', False):
            with st.spinner("Generating ZIP archive..."):
                st.session_state.zip_bytes = generate_zip_archive(st.session_state.get('report_elements', []), st.session_state.get('am_data_source_message_for_report', 'N/A'))
                st.session_state.trigger_zip_regeneration = False
        if 'zip_bytes' in st.session_state and st.session_state.zip_bytes:
            st.download_button(label="Download All Data & Plots (.zip)", data=st.session_state.zip_bytes, file_name="SlytheRINs_All_Data_Plots.zip", mime="application/zip", key="download_zip_button")


# --- Initialize session state for report elements ---
if 'report_elements' not in st.session_state: st.session_state.report_elements = []
if 'am_data_source_message_for_report' not in st.session_state: st.session_state.am_data_source_message_for_report = ""


# --- Data Loading Logic ---
edge_dataframes = []
edge_df_names = [] 
new_data_loaded_flag = False

if use_example_files:
    if not st.session_state.get('using_examples_for_report', False) or not edge_dataframes: 
        st.session_state.report_elements = []; st.session_state.reports_ready = False
        st.session_state.pop('pdf_bytes', None); st.session_state.pop('zip_bytes', None)
        st.session_state.trigger_pdf_regeneration = True; st.session_state.trigger_zip_regeneration = True
        new_data_loaded_flag = True
        st.session_state.using_examples_for_report = True; st.session_state.using_uploads_for_report = False
    example_dir = 'example_data'
    if os.path.isdir(example_dir):
        example_files_paths = [os.path.join(example_dir, f) for f in os.listdir(example_dir) if f.endswith('.txt')]
        if not example_files_paths: st.sidebar.warning(f"No '.txt' example files found in '{example_dir}'.")
        for fp in example_files_paths:
            try: df = pd.read_csv(fp, sep='\t'); edge_dataframes.append(df); edge_df_names.append(os.path.basename(fp))
            except Exception as e: st.sidebar.error(f"Error reading example {fp}: {e}")
    else: st.sidebar.error(f"Example directory '{example_dir}' not found.")
elif uploaded_files: 
    if not st.session_state.get('using_uploads_for_report', False) or not edge_dataframes:
        st.session_state.report_elements = []; st.session_state.reports_ready = False
        st.session_state.pop('pdf_bytes', None); st.session_state.pop('zip_bytes', None)
        st.session_state.trigger_pdf_regeneration = True; st.session_state.trigger_zip_regeneration = True
        new_data_loaded_flag = True
        st.session_state.using_uploads_for_report = True; st.session_state.using_examples_for_report = False
    for ufo in uploaded_files:
        try: df = pd.read_csv(ufo, sep='\t'); edge_dataframes.append(df); edge_df_names.append(ufo.name)
        except Exception as e: st.sidebar.error(f"Error reading uploaded {ufo.name}: {e}")
else: 
    if st.session_state.get('using_examples_for_report', False) or st.session_state.get('using_uploads_for_report', False):
        st.session_state.report_elements = []; st.session_state.reports_ready = False
        st.session_state.pop('pdf_bytes', None); st.session_state.pop('zip_bytes', None)
        st.session_state.trigger_pdf_regeneration = True; st.session_state.trigger_zip_regeneration = True
    st.session_state.using_examples_for_report = False; st.session_state.using_uploads_for_report = False

# AlphaMissense Data Loading
am_df_final = pd.DataFrame()
am_data_source_message = ""
if 'downloaded_am_data' not in st.session_state: st.session_state.downloaded_am_data = None
if 'downloaded_uniprot_id' not in st.session_state: st.session_state.downloaded_uniprot_id = None
am_class_mapping = {'LPath': 'pathogenic', 'LBen': 'benign', 'Amb': 'ambiguous'}

if download_am_button and uniprot_id_input:
    uniprot_id = uniprot_id_input.strip().upper()
    am_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-aa-substitutions.csv"
    try:
        with st.spinner(f"Downloading AlphaMissense data for {uniprot_id}..."):
            response = requests.get(am_url, timeout=30); response.raise_for_status()  
            temp_df = pd.read_csv(io.StringIO(response.text))
            if 'protein_variant' in temp_df.columns and 'am_pathogenicity' in temp_df.columns and 'am_class' in temp_df.columns:
                am_df_final = temp_df.rename(columns={'protein_variant': 'Mutation', 'am_pathogenicity': 'AM_Score', 'am_class': 'AM_Classification'})
                am_df_final['ProteinID'] = uniprot_id 
                am_df_final = am_df_final[['ProteinID', 'Mutation', 'AM_Score', 'AM_Classification']] 
                am_df_final['AM_Classification'] = am_df_final['AM_Classification'].map(am_class_mapping).fillna(am_df_final['AM_Classification'])
                st.session_state.downloaded_am_data = am_df_final.copy(); st.session_state.downloaded_uniprot_id = uniprot_id
                am_data_source_message = f"Using AlphaMissense data downloaded for {uniprot_id}."
                st.session_state.am_data_source_message_for_report = am_data_source_message
                st.sidebar.success(f"Successfully downloaded data for {uniprot_id}.")
                st.session_state.trigger_pdf_regeneration = True; st.session_state.trigger_zip_regeneration = True
            else:
                st.sidebar.error(f"Downloaded CSV for {uniprot_id} has unexpected columns."); st.session_state.downloaded_am_data = None; st.session_state.downloaded_uniprot_id = None
    except Exception as e: st.sidebar.error(f"Error with AM data for {uniprot_id}: {e}"); st.session_state.downloaded_am_data = None; st.session_state.downloaded_uniprot_id = None

if st.session_state.downloaded_am_data is not None and not st.session_state.downloaded_am_data.empty and \
   (not (download_am_button and uniprot_id_input) or (uniprot_id_input.strip().upper() == st.session_state.downloaded_uniprot_id)):
    am_df_final = st.session_state.downloaded_am_data
    if not am_data_source_message: am_data_source_message = f"Using previously downloaded AM data for {st.session_state.downloaded_uniprot_id}."; st.session_state.am_data_source_message_for_report = am_data_source_message
elif uploaded_file_am:
    try:
        temp_df = pd.read_csv(uploaded_file_am); uploaded_file_am.seek(0) 
        if 'protein_variant' in temp_df.columns and 'am_pathogenicity' in temp_df.columns and 'am_class' in temp_df.columns:
            am_df_final = temp_df.rename(columns={'protein_variant': 'Mutation', 'am_pathogenicity': 'AM_Score', 'am_class': 'AM_Classification'})
            match = re.search(r'AF-([A-Z0-9]+)-F1', uploaded_file_am.name, re.IGNORECASE); protein_id_from_filename = match.group(1) if match else "Unknown"
            am_df_final['ProteinID'] = protein_id_from_filename
            am_df_final = am_df_final[['ProteinID', 'Mutation', 'AM_Score', 'AM_Classification']]
            am_df_final['AM_Classification'] = am_df_final['AM_Classification'].map(am_class_mapping).fillna(am_df_final['AM_Classification'])
            am_data_source_message = f"Using uploaded file: {uploaded_file_am.name} (Detected AlphaFold DB CSV format)."
        else: 
            temp_df = pd.read_csv(uploaded_file_am, sep=None, header=None, engine='python', on_bad_lines='skip'); uploaded_file_am.seek(0)
            if temp_df.shape[1] == 4:
                 am_df_final = temp_df.copy(); am_df_final.columns = ['ProteinID', 'Mutation', 'AM_Score', 'AM_Classification']
                 am_data_source_message = f"Using uploaded file: {uploaded_file_am.name} (Assumed 4-column format)."
            else: st.sidebar.error(f"Uploaded file {uploaded_file_am.name} not recognized."); am_df_final = pd.DataFrame()
        if not am_df_final.empty: st.session_state.downloaded_am_data = None; st.session_state.downloaded_uniprot_id = None; st.session_state.am_data_source_message_for_report = am_data_source_message; st.session_state.trigger_pdf_regeneration = True; st.session_state.trigger_zip_regeneration = True
    except Exception as e: st.sidebar.error(f"Error reading uploaded AM file {uploaded_file_am.name}: {e}"); am_df_final = pd.DataFrame()
elif use_example_file_am:
    am_example_path = 'example_data/P03891-pred.tsv' 
    if os.path.exists(am_example_path):
        try:
            am_df_final = pd.read_csv(am_example_path, sep='\t', header=None) 
            if am_df_final.shape[1] == 4:
                am_df_final.columns = ['ProteinID', 'Mutation', 'AM_Score', 'AM_Classification']
                am_data_source_message = "Using example AlphaMissense file (P03891-pred.tsv)."
                st.session_state.downloaded_am_data = None; st.session_state.downloaded_uniprot_id = None
                st.session_state.am_data_source_message_for_report = am_data_source_message
                st.session_state.trigger_pdf_regeneration = True; st.session_state.trigger_zip_regeneration = True
            else: st.sidebar.error(f"Example AM file {am_example_path} does not have 4 columns."); am_df_final = pd.DataFrame()
        except Exception as e: st.sidebar.error(f"Error reading example AM file {am_example_path}: {e}")
    else: st.sidebar.warning(f"Example AM file '{am_example_path}' not found.")

# --- Main Content Area ---
st.header('**SlytheRINs Dashboard**') 
tab1, tab2, tab3, tab_docs = st.tabs(["RIN Comparison & Network Complexity", "Chemical Interactions", "AlphaMissense Integration", "Help & Documentation"]) 

all_nodes_set_tab1 = set() 
valid_dfs_for_processing_tab1 = True if edge_dataframes else False

with tab1:
    st.subheader('RIN Comparison & Network Complexity')
    if edge_dataframes:
        all_degrees_list = [] 
        all_graphs_nx = []    
        all_degrees_tab1, all_betweenness_tab1, all_clustering_tab1, all_triangles_tab1 = [], [], [], [] 
        all_eigenvector_centrality_tab1 = [] 
        assortativity_values_tab1 = []
        if new_data_loaded_flag: st.session_state.report_elements = []; st.session_state.trigger_pdf_regeneration = True; st.session_state.trigger_zip_regeneration = True
        for i, df_iter in enumerate(edge_dataframes): 
            if 'NodeId1' not in df_iter.columns or 'NodeId2' not in df_iter.columns:
                st.error(f"File '{edge_df_names[i]}' missing key columns."); valid_dfs_for_processing_tab1 = False; break
            try:
                G = nx.from_pandas_edgelist(df_iter, 'NodeId1', 'NodeId2'); all_graphs_nx.append(G); all_nodes_set_tab1.update(G.nodes())
                current_degrees = dict(G.degree()); all_degrees_list.append(current_degrees); all_degrees_tab1.append(current_degrees)      
                all_betweenness_tab1.append(nx.betweenness_centrality(G)); all_clustering_tab1.append(nx.clustering(G)); all_triangles_tab1.append(nx.triangles(G))
                try: 
                    if nx.is_connected(G): 
                        all_eigenvector_centrality_tab1.append(nx.eigenvector_centrality_numpy(G))
                    else: 
                        if G.number_of_nodes() > 0 : 
                            largest_cc_nodes = max(nx.connected_components(G), key=len, default=set())
                            if largest_cc_nodes: 
                                subgraph = G.subgraph(largest_cc_nodes)
                                if subgraph.number_of_nodes() > 1: 
                                    eigenvector_centrality_subgraph = nx.eigenvector_centrality_numpy(subgraph)
                                    full_graph_eigenvector = {node: eigenvector_centrality_subgraph.get(node, 0) for node in G.nodes()}
                                    all_eigenvector_centrality_tab1.append(full_graph_eigenvector)
                                else:
                                    all_eigenvector_centrality_tab1.append({node: 0 for node in G.nodes()}) 
                            else: 
                                all_eigenvector_centrality_tab1.append({node: 0 for node in G.nodes()})
                        else: 
                            all_eigenvector_centrality_tab1.append({})
                except (nx.NetworkXError, nx.NetworkXPointlessConcept, Exception) as e: 
                    st.warning(f"Could not calculate Eigenvector Centrality for graph '{edge_df_names[i]}'. Error: {e}")
                    all_eigenvector_centrality_tab1.append({node: 0 for node in G.nodes()} if G.nodes() else {}) 
                try: assortativity = nx.degree_assortativity_coefficient(G); assortativity_values_tab1.append(assortativity)
                except: st.warning(f"Assortativity error for '{edge_df_names[i]}'."); assortativity_values_tab1.append(np.nan)
            except Exception as e: st.error(f"Error processing '{edge_df_names[i]}': {e}"); valid_dfs_for_processing_tab1 = False; break
        if valid_dfs_for_processing_tab1 and all_nodes_set_tab1:
            sorted_residue_numbers = sorted(list(all_nodes_set_tab1), key=lambda x: int(x.split(':')[1]))
            formatted_labels = [f'{x.split(":")[1]}:{x.split(":")[3]}' for x in sorted_residue_numbers]
            node_to_formatted_label = dict(zip(sorted_residue_numbers, formatted_labels))
            all_degree_values = pd.DataFrame(all_degrees_tab1).reindex(columns=sorted_residue_numbers).fillna(0)
            all_betweenness_values = pd.DataFrame(all_betweenness_tab1).reindex(columns=sorted_residue_numbers).fillna(0)
            all_clustering_values = pd.DataFrame(all_clustering_tab1).reindex(columns=sorted_residue_numbers).fillna(0)
            all_eigenvector_values = pd.DataFrame(all_eigenvector_centrality_tab1).reindex(columns=sorted_residue_numbers).fillna(0) 
            triangle_counts_df = pd.DataFrame(all_triangles_tab1).reindex(columns=sorted_residue_numbers).fillna(0)
            
            degree_variance = all_degree_values.var(axis=0); 
            betweenness_variance = all_betweenness_values.var(axis=0); 
            clustering_variance = all_clustering_values.var(axis=0)
            eigenvector_variance = all_eigenvector_values.var(axis=0) 

            mean_degree = all_degree_values.mean(axis=0); std_degree = all_degree_values.std(axis=0); neg_std_degree = -std_degree
            mean_triangles = triangle_counts_df.mean(axis=0); std_triangles = triangle_counts_df.std(axis=0)
            mean_eigenvector = all_eigenvector_values.mean(axis=0) 
            std_eigenvector = all_eigenvector_values.std(axis=0)   

            N_graphs = len(all_degree_values) 
            st.session_state.report_elements.append({"type": "metric", "label": "Number of Networks Loaded", "value": N_graphs, "part": 1})
            st.session_state.report_elements.append({"type": "metric", "label": "Total Unique Residues", "value": len(all_nodes_set_tab1), "part": 1})
            st.metric(label="Number of Networks Loaded", value=N_graphs); st.metric(label="Total Unique Residues Identified", value=len(all_nodes_set_tab1))
            df_degree_variance_dl = pd.DataFrame({'Residue': sorted_residue_numbers, 'Formatted_Label': formatted_labels, 'Degree_Variance': degree_variance.loc[sorted_residue_numbers]})
            df_betweenness_variance_dl = pd.DataFrame({'Residue': sorted_residue_numbers, 'Formatted_Label': formatted_labels, 'Betweenness_Variance': betweenness_variance.loc[sorted_residue_numbers]})
            df_clustering_variance_dl = pd.DataFrame({'Residue': sorted_residue_numbers, 'Formatted_Label': formatted_labels, 'Clustering_Variance': clustering_variance.loc[sorted_residue_numbers]})
            df_eigenvector_variance_dl = pd.DataFrame({'Residue': sorted_residue_numbers, 'Formatted_Label': formatted_labels, 'Eigenvector_Variance': eigenvector_variance.loc[sorted_residue_numbers]})
            df_mean_std_degree_dl = pd.DataFrame({'Residue': sorted_residue_numbers, 'Formatted_Label': formatted_labels, 'Mean_Degree': mean_degree.loc[sorted_residue_numbers], 'STD_Degree': std_degree.loc[sorted_residue_numbers]})
            df_mean_std_eigenvector_dl = pd.DataFrame({'Residue': sorted_residue_numbers, 'Formatted_Label': formatted_labels, 'Mean_Eigenvector': mean_eigenvector.loc[sorted_residue_numbers], 'STD_Eigenvector': std_eigenvector.loc[sorted_residue_numbers]})
            df_triangle_counts_dl = pd.DataFrame({'Residue': sorted_residue_numbers, 'Formatted_Label': formatted_labels, 'Mean_Triangles': mean_triangles.loc[sorted_residue_numbers], 'STD_Triangles': std_triangles.loc[sorted_residue_numbers]})
            df_assortativity_dl = pd.DataFrame({'Graph_Index': list(range(1, len(assortativity_values_tab1) + 1)), 'Graph_Name': edge_df_names[:len(assortativity_values_tab1)], 'Assortativity': assortativity_values_tab1})
            
            st.markdown("#### Variance Plots (Across All Networks)")
            col_var1, col_var2 = st.columns(2) 
            with col_var1:
                caption_deg_var = "Variation in connectivity (degree) of residues across conformations. High variance highlights residues with constant change in the number of interactions."
                st.write("##### Degree Variance"); st.caption(caption_deg_var)
                fig_degree = go.Figure(go.Scatter(x=formatted_labels, y=degree_variance.loc[sorted_residue_numbers], mode='lines+markers'))
                fig_degree.update_layout(title='Degree Variance', height=400); st.plotly_chart(fig_degree, use_container_width=True)
                st.session_state.report_elements.append({"type": "plot", "fig": fig_degree, "title": "Degree Variance", "caption": caption_deg_var, "part": 1})
                st.session_state.report_elements.append({"type": "df_text", "df": df_degree_variance_dl, "title": "Degree_Variance_Data", "part": 1}) 
                dl_c1_1, dl_c2_1 = st.columns(2)
                with dl_c1_1: st.download_button("Data (TSV)", convert_df_to_tsv(df_degree_variance_dl), "degree_variance.tsv", key="dl_deg_var_tsv_tab1")
                with dl_c2_1: st.download_button("Plot (PNG)", convert_fig_to_png(fig_degree), "degree_variance_plot.png", key="dl_deg_var_png_tab1")
            with col_var2:
                caption_bet_var = "Variation in betweenness centrality of residues across conformations. High variance indicates residues that frequently change their role as bridges or bottlenecks."
                st.write("##### Betweenness-centrality Variance"); st.caption(caption_bet_var)
                fig_betweenness = go.Figure(go.Scatter(x=formatted_labels, y=betweenness_variance.loc[sorted_residue_numbers], mode='lines+markers'))
                fig_betweenness.update_layout(title='Betweenness Variance', height=400); st.plotly_chart(fig_betweenness, use_container_width=True)
                st.session_state.report_elements.append({"type": "plot", "fig": fig_betweenness, "title": "Betweenness Variance", "caption": caption_bet_var, "part": 1})
                st.session_state.report_elements.append({"type": "df_text", "df": df_betweenness_variance_dl, "title": "Betweenness_Variance_Data", "part": 1})
                dl_c1_2, dl_c2_2 = st.columns(2)
                with dl_c1_2: st.download_button("Data (TSV)", convert_df_to_tsv(df_betweenness_variance_dl), "betweenness_variance.tsv", key="dl_bet_var_tsv_tab1")
                with dl_c2_2: st.download_button("Plot (PNG)", convert_fig_to_png(fig_betweenness), "betweenness_variance_plot.png", key="dl_bet_var_png_tab1")
            
            col_var3, col_var4 = st.columns(2) 
            with col_var3:
                caption_clust_var = "Variation in clustering coefficient of residues across conformations. High variance indicates residues that frequently change their local connectivity patterns."
                st.write("##### Clustering Coefficient Variance"); st.caption(caption_clust_var)
                fig_clustering = go.Figure(go.Scatter(x=formatted_labels, y=clustering_variance.loc[sorted_residue_numbers], mode='lines+markers'))
                fig_clustering.update_layout(title='Clustering Coeff. Variance', height=400); st.plotly_chart(fig_clustering, use_container_width=True)
                st.session_state.report_elements.append({"type": "plot", "fig": fig_clustering, "title": "Clustering Coefficient Variance", "caption": caption_clust_var, "part": 1})
                st.session_state.report_elements.append({"type": "df_text", "df": df_clustering_variance_dl, "title": "Clustering_Coefficient_Variance_Data", "part": 1})
                dl_c1_3, dl_c2_3 = st.columns(2)
                with dl_c1_3: st.download_button("Data (TSV)", convert_df_to_tsv(df_clustering_variance_dl), "clustering_variance.tsv", key="dl_clust_var_tsv_tab1")
                with dl_c2_3: st.download_button("Plot (PNG)", convert_fig_to_png(fig_clustering), "clustering_variance_plot.png", key="dl_clust_var_png_tab1")
            with col_var4: 
                caption_eig_var = "Variation in eigenvector centrality. High variance indicates residues whose influence (considering connections to other influential residues) changes significantly."
                st.write("##### Eigenvector Centrality Variance"); st.caption(caption_eig_var)
                fig_eigenvector_var = go.Figure(go.Scatter(x=formatted_labels, y=eigenvector_variance.loc[sorted_residue_numbers], mode='lines+markers'))
                fig_eigenvector_var.update_layout(title='Eigenvector Centrality Variance', height=400); st.plotly_chart(fig_eigenvector_var, use_container_width=True)
                st.session_state.report_elements.append({"type": "plot", "fig": fig_eigenvector_var, "title": "Eigenvector Centrality Variance", "caption": caption_eig_var, "part": 1})
                st.session_state.report_elements.append({"type": "df_text", "df": df_eigenvector_variance_dl, "title": "Eigenvector_Variance_Data", "part": 1})
                dl_c1_eig, dl_c2_eig = st.columns(2)
                with dl_c1_eig: st.download_button("Data (TSV)", convert_df_to_tsv(df_eigenvector_variance_dl), "eigenvector_variance.tsv", key="dl_eig_var_tsv_tab1")
                with dl_c2_eig: st.download_button("Plot (PNG)", convert_fig_to_png(fig_eigenvector_var), "eigenvector_variance_plot.png", key="dl_eig_var_png_tab1")


            st.markdown("#### Mean & STD Plots (Across All Networks)") 
            col_mean_std1, col_mean_std2, col_mean_std3 = st.columns(3) 
            with col_mean_std1: 
                caption_tri_counts = "Display each residue’s average (mean) involvement in triangles (3 mutually connected residues) and its variability across conformations. Shaded areas display the standard deviation, indicating residue participation in local clusters."
                st.write("##### Triangle Counts"); st.caption(caption_tri_counts)
                fig_triangles = go.Figure([go.Scatter(x=formatted_labels, y=mean_triangles.loc[sorted_residue_numbers], mode='lines+markers', name='Mean'), go.Scatter(x=formatted_labels, y=(mean_triangles - std_triangles).loc[sorted_residue_numbers], fill=None, mode='lines', line_color='lightblue', name='-STD', showlegend=False), go.Scatter(x=formatted_labels, y=(mean_triangles + std_triangles).loc[sorted_residue_numbers], fill='tonexty', mode='lines', line_color='lightblue', name='+STD', showlegend=False)])
                fig_triangles.update_layout(title='Mean & STD of Triangle Counts', height=400); st.plotly_chart(fig_triangles, use_container_width=True)
                st.session_state.report_elements.append({"type": "plot", "fig": fig_triangles, "title": "Mean & STD of Triangle Counts", "caption": caption_tri_counts, "part": 1})
                st.session_state.report_elements.append({"type": "df_text", "df": df_triangle_counts_dl, "title": "Triangle_Counts_Data", "part": 1})
                dl_c1_4, dl_c2_4 = st.columns(2)
                with dl_c1_4: st.download_button("Data (TSV)", convert_df_to_tsv(df_triangle_counts_dl), "triangle_counts.tsv",key="dl_tri_counts_tsv_tab1")
                with dl_c2_4: st.download_button("Plot (PNG)", convert_fig_to_png(fig_triangles), "triangle_counts_plot.png", key="dl_tri_counts_png_tab1")
            with col_mean_std2:
                caption_deg_sd = "Display each residue’s average (mean) degree value and variability across distinct conformations. Shaded areas display the standard deviation, indicating residue connectivity stability."
                st.write("##### Degree Values"); st.caption(caption_deg_sd)
                fig_degree_sd = go.Figure([go.Scatter(x=formatted_labels, y=mean_degree.loc[sorted_residue_numbers], mode='lines+markers', name='Mean'), go.Scatter(x=formatted_labels, y=(mean_degree - std_degree).loc[sorted_residue_numbers], fill=None, mode='lines', line_color='lightblue', name='-STD', showlegend=False), go.Scatter(x=formatted_labels, y=(mean_degree + std_degree).loc[sorted_residue_numbers], fill='tonexty', mode='lines', line_color='lightblue', name='+STD', showlegend=False)])
                fig_degree_sd.update_layout(title='Mean & STD of Degree Values', height=400); st.plotly_chart(fig_degree_sd, use_container_width=True)
                st.session_state.report_elements.append({"type": "plot", "fig": fig_degree_sd, "title": "Mean & STD of Degree Values", "caption": caption_deg_sd, "part": 1})
                st.session_state.report_elements.append({"type": "df_text", "df": df_mean_std_degree_dl, "title": "Mean & STD Degree Data", "part": 1})
                dl_c1_5, dl_c2_5 = st.columns(2)
                with dl_c1_5: st.download_button("Data (TSV)", convert_df_to_tsv(df_mean_std_degree_dl), "mean_std_degree.tsv", key="dl_mean_std_deg_tsv_tab1")
                with dl_c2_5: st.download_button("Plot (PNG)", convert_fig_to_png(fig_degree_sd), "mean_std_degree_plot.png", key="dl_mean_std_deg_png_tab1")
            with col_mean_std3: 
                caption_eig_sd = "Displays each residue’s average (mean) eigenvector centrality and its variability. Eigenvector centrality measures a node's influence based on the centrality of its neighbors."
                st.write("##### Eigenvector Centrality"); st.caption(caption_eig_sd)
                fig_eigenvector_sd = go.Figure([
                    go.Scatter(x=formatted_labels, y=mean_eigenvector.loc[sorted_residue_numbers], mode='lines+markers', name='Mean Eigenvector'),
                    go.Scatter(x=formatted_labels, y=(mean_eigenvector - std_eigenvector).loc[sorted_residue_numbers], fill=None, mode='lines', line_color='lightcoral', name='-STD Eigenvector', showlegend=False), 
                    go.Scatter(x=formatted_labels, y=(mean_eigenvector + std_eigenvector).loc[sorted_residue_numbers], fill='tonexty', mode='lines', line_color='lightcoral', name='+STD Eigenvector', showlegend=False)
                ])
                fig_eigenvector_sd.update_layout(title='Mean & STD of Eigenvector Centrality', height=400, yaxis_title="Eigenvector Centrality"); st.plotly_chart(fig_eigenvector_sd, use_container_width=True)
                st.session_state.report_elements.append({"type": "plot", "fig": fig_eigenvector_sd, "title": "Mean & STD of Eigenvector Centrality", "caption": caption_eig_sd, "part": 1})
                st.session_state.report_elements.append({"type": "df_text", "df": df_mean_std_eigenvector_dl, "title": "Mean_STD_Eigenvector_Data", "part": 1})
                dl_c1_eig_sd, dl_c2_eig_sd = st.columns(2)
                with dl_c1_eig_sd: st.download_button("Data (TSV)", convert_df_to_tsv(df_mean_std_eigenvector_dl), "mean_std_eigenvector.tsv", key="dl_mean_std_eig_tsv_tab1")
                with dl_c2_eig_sd: st.download_button("Plot (PNG)", convert_fig_to_png(fig_eigenvector_sd), "mean_std_eigenvector_plot.png", key="dl_mean_std_eig_png_tab1")

            st.markdown("#### Top 10 Hubs (Mean Eigenvector Centrality)")
            mean_eigenvector_series = all_eigenvector_values.mean(axis=0)
            top_10_hubs_eigenvector = mean_eigenvector_series.nlargest(10)
            
            hub_data = []
            for residue_raw, centrality_score in top_10_hubs_eigenvector.items():
                hub_label = "Unknown" 
                mean_deg_for_hub = mean_degree.get(residue_raw, 0) 
                try:
                    parts = residue_raw.split(':')
                    if len(parts) >= 3: 
                        res_id = parts[2] 
                        pos = parts[1]    
                        hub_label = f"{res_id}{pos}" 
                except IndexError: pass 
                except Exception as e: print(f"Error formatting hub label for {residue_raw}: {e}")
                hub_data.append({
                    'Raw NodeId': residue_raw, # Changed order, Raw NodeId first
                    'Hub Residue (AAAPOS)': hub_label, 
                    'Mean Degree': mean_deg_for_hub, 
                    'Mean Eigenvector Centrality': centrality_score
                })
            
            df_top_hubs_eigenvector_dl = pd.DataFrame(hub_data)
            st.dataframe(df_top_hubs_eigenvector_dl[['Raw NodeId', 'Hub Residue (AAAPOS)', 'Mean Degree', 'Mean Eigenvector Centrality']].style.format({'Mean Eigenvector Centrality': "{:.4f}", 'Mean Degree': "{:.2f}"}))
            st.session_state.report_elements.append({"type": "df_text", "df": df_top_hubs_eigenvector_dl, "title": "Top_10_Hubs_by_Mean_Eigenvector_Centrality", "part": 1}) 
            hub_table_string = df_top_hubs_eigenvector_dl[['Raw NodeId','Hub Residue (AAAPOS)', 'Mean Degree', 'Mean Eigenvector Centrality']].to_string(index=False)
            st.session_state.report_elements.append({"type": "hub_table_text", "df_string": hub_table_string, "title": "Top 10 Hubs by Mean Eigenvector Centrality", "part": 1}) 
            st.download_button("Download Top 10 Hubs Data (TSV)", convert_df_to_tsv(df_top_hubs_eigenvector_dl), "top_10_hubs_eigenvector.tsv", key="dl_top_hubs_eigenvector_tsv")


            st.markdown("#### Degree Distribution per Residue (Box Plot - Across All Networks)")
            caption_deg_boxplot = "Shows the distribution (median, quartiles, outliers) of degree values for each residue across all analyzed networks."
            st.caption(caption_deg_boxplot)
            df_degree_boxplot_data = all_degree_values.melt(var_name='Residue_Raw', value_name='Degree')
            df_degree_boxplot_data['Residue'] = df_degree_boxplot_data['Residue_Raw'].map(node_to_formatted_label)
            df_degree_boxplot_data['Residue_Raw_Cat'] = pd.Categorical(df_degree_boxplot_data['Residue_Raw'], categories=sorted_residue_numbers, ordered=True)
            df_degree_boxplot_data = df_degree_boxplot_data.sort_values('Residue_Raw_Cat')
            if not df_degree_boxplot_data.empty:
                fig_degree_boxplot = px.box(df_degree_boxplot_data, x='Residue', y='Degree', title='Distribution of Degree Values per Residue')
                fig_degree_boxplot.update_layout(height=500); st.plotly_chart(fig_degree_boxplot, use_container_width=True)
                st.session_state.report_elements.append({"type": "plot", "fig": fig_degree_boxplot, "title": "Degree Distribution per Residue (Box Plot)", "caption": caption_deg_boxplot, "part": 1})
                st.session_state.report_elements.append({"type": "df_text", "df": df_degree_boxplot_data[['Residue', 'Degree', 'Residue_Raw']], "title": "Degree Boxplot Data", "part": 1})
                col1_bp, col2_bp = st.columns(2) 
                with col1_bp: st.download_button("Data (TSV)", convert_df_to_tsv(df_degree_boxplot_data[['Residue', 'Degree', 'Residue_Raw']]), "degree_boxplot_data.tsv", key="dl_deg_boxplot_tsv_tab1")
                with col2_bp: st.download_button("Plot (PNG)", convert_fig_to_png(fig_degree_boxplot), "degree_boxplot.png", key="dl_deg_boxplot_png_tab1")
            else: st.warning("Could not generate data for the degree distribution box plot.")
            
            st.markdown("#### Degree Standard Deviation (DSD) Plot")
            caption_dsd = "A SlytheRINs' specific metric. It displays each residue's Degree Standard Deviation (DSD), highlighting how variable its number of interactions is across conformations. Zero or low DSD values (shaded regions of the plot) indicate stable connections, while high DSD values indicates dynamic behavior."
            st.caption(caption_dsd)
            fig_dsd = go.Figure([ go.Scatter(x=formatted_labels, y=std_degree.loc[sorted_residue_numbers], mode='lines+markers', name='STD'), go.Scatter(x=formatted_labels, y=neg_std_degree.loc[sorted_residue_numbers], mode='lines+markers', name='-STD')])
            if formatted_labels: fig_dsd.add_shape(type='rect', x0=formatted_labels[0], x1=formatted_labels[-1], y0=-0.5, y1=0.5, fillcolor='lightgrey', opacity=0.5, layer='below', line_width=0); fig_dsd.add_shape(type='rect', x0=formatted_labels[0], x1=formatted_labels[-1], y0=-1.0, y1=1.0, fillcolor='lightgrey', opacity=0.3, layer='below', line_width=0)
            fig_dsd.update_layout(title='Degree Standard Deviation (DSD)', height=500); st.plotly_chart(fig_dsd, use_container_width=True)
            st.session_state.report_elements.append({"type": "plot", "fig": fig_dsd, "title": "Degree Standard Deviation (DSD)", "caption": caption_dsd, "part": 1})
            st.download_button("Download DSD Plot (PNG)", convert_fig_to_png(fig_dsd), "dsd_plot.png", key="dl_dsd_png_tab1")

            st.markdown("#### Assortativity & Significance (Per Network)")
            assort_valid = [v for v in assortativity_values_tab1 if not np.isnan(v)] 
            if len(assort_valid) > 1:
                stat_assort, p_assort = st_sci.ttest_1samp(assort_valid, popmean=0); signif_assort = p_assort < 0.05; texto_assort = f"{'* ' if signif_assort else ''}p = {p_assort:.3e}"
                assort_summary = f"Overall Assortativity T-Test (vs 0): Mean = {np.mean(assort_valid):.4f} | t = {stat_assort:.3f} | p = {p_assort:.3e} {'(*p < 0.05*)' if signif_assort else ''}"
                st.write(f"**{assort_summary}**"); st.session_state.report_elements.append({"type": "text_summary", "title": "Overall Assortativity T-Test", "content": assort_summary, "part": 1})
            else: st.warning("Not enough valid assortativity values for overall t-test."); texto_assort = "N/A"
            fig_assortativity = go.Figure(go.Scatter(x=list(range(1, len(assortativity_values_tab1) + 1)), y=assortativity_values_tab1, mode='lines+markers'))
            fig_assortativity.add_hline(y=0, line_dash="dash", annotation_text="zero"); fig_assortativity.add_annotation(xref="paper", yref="paper", x=0.98, y=0.98, text=texto_assort, showarrow=False, font=dict(size=10))
            fig_assortativity.update_layout(title="Assortativity Coefficient per Network", height=400); st.plotly_chart(fig_assortativity, use_container_width=True)
            st.session_state.report_elements.append({"type": "plot", "fig": fig_assortativity, "title": "Assortativity Coefficient per Network", "caption": "Degree assortativity coefficient for each network. Positive values indicate hubs connect to hubs; negative values indicate hubs connect to low-degree nodes.", "part": 1})
            st.session_state.report_elements.append({"type": "df_text", "df": df_assortativity_dl, "title": "Assortativity Data", "part": 1})
            col1_6, col2_6 = st.columns(2) 
            with col1_6: st.download_button("Data (TSV)", convert_df_to_tsv(df_assortativity_dl), "assortativity.tsv", key="dl_assort_tsv_tab1")
            with col2_6: st.download_button("Plot (PNG)", convert_fig_to_png(fig_assortativity), "assortativity_plot.png", key="dl_assort_png_tab1")
            
            st.markdown("#### Significance Plots (Across All Networks)")
            col_sig1, col_sig2 = st.columns(2)
            with col_sig1:
                caption_paired_t = "Tests if a residue's degree is significantly different from the average degree of all *other* residues across the networks. Low p-values (e.g. < 0.05, below red line) suggest significance."
                st.markdown("###### Degree Significance (Paired T-test vs. Others)"); st.caption(caption_paired_t)
                paired_t_p_values = []
                if N_graphs > 1:
                    for residue in sorted_residue_numbers:
                        x_degrees = all_degree_values[residue].values; overall_means_excluding = []
                        for idx, row in all_degree_values.iterrows():
                            if residue in row.index: total = row.sum() - row[residue]; count = len(row) -1 if residue in row.index else len(row) 
                            else: total = row.sum(); count = len(row)
                            overall_means_excluding.append(total / count if count > 0 else 0)
                        if len(x_degrees) == len(overall_means_excluding) and len(x_degrees) > 1 and (np.std(x_degrees) > 0 or np.std(overall_means_excluding) > 0):
                            try: t_stat, p_val_t = st_sci.ttest_rel(x_degrees, overall_means_excluding); paired_t_p_values.append(p_val_t)
                            except ValueError: paired_t_p_values.append(np.nan)
                        else: paired_t_p_values.append(np.nan) 
                    df_significance_paired_t = pd.DataFrame({"Residue": sorted_residue_numbers, "Formatted_Label": formatted_labels, "Degree_Variance": degree_variance.loc[sorted_residue_numbers], "Paired_T_p_value": paired_t_p_values})
                    fig_pvalues_paired = go.Figure(go.Scatter(x=formatted_labels, y=paired_t_p_values, mode='lines+markers'))
                    if formatted_labels: fig_pvalues_paired.add_shape(type="line", x0=formatted_labels[0], x1=formatted_labels[-1], y0=0.05, y1=0.05, line=dict(color="red", dash="dash"))
                    fig_pvalues_paired.update_layout(title="P-Values (Paired T-Test)", height=400); st.plotly_chart(fig_pvalues_paired, use_container_width=True)
                    st.session_state.report_elements.append({"type": "plot", "fig": fig_pvalues_paired, "title": "Degree Significance (Paired T-Test)", "caption": caption_paired_t, "part": 1})
                    st.session_state.report_elements.append({"type": "df_text", "df": df_significance_paired_t, "title": "Paired T-Test Data", "part": 1})
                    dl_c1_sub_7, dl_c2_sub_7 = st.columns(2) 
                    with dl_c1_sub_7: st.download_button("Data (TSV)", convert_df_to_tsv(df_significance_paired_t), "residue_significance_paired_t.tsv", key="dl_paired_t_tsv_tab1")
                    with dl_c2_sub_7: st.download_button("Plot (PNG)", convert_fig_to_png(fig_pvalues_paired), "paired_t_test_plot.png", key="dl_paired_t_png_tab1")
                else: st.warning("Need >1 network.")
            with col_sig2:
                caption_chi2_var = "Tests if a residue's degree variation *within itself* across networks is significantly different from a random (Poisson) process. Low p-values suggest non-random variation."
                st.markdown("###### Degree Variation Significance (vs. Poisson)"); st.caption(caption_chi2_var)
                variation_p_values = []
                if N_graphs > 1:
                    for residue in sorted_residue_numbers:
                        degrees_for_residue = all_degree_values[residue].values; mean_deg_residue = np.mean(degrees_for_residue); var_deg_residue = np.var(degrees_for_residue, ddof=1) 
                        if var_deg_residue == 0 or mean_deg_residue == 0: p_value_chi2 = 1.0
                        else: chi2_stat = (N_graphs - 1) * var_deg_residue / mean_deg_residue; cdf_val = st_sci.chi2.cdf(chi2_stat, N_graphs - 1); p_value_chi2 = 2 * min(cdf_val, 1 - cdf_val) 
                        variation_p_values.append(p_value_chi2)
                    df_variation_pvals_dl = pd.DataFrame({"Residue": sorted_residue_numbers, "Formatted_Label": formatted_labels, "Degree_Variation_P_Value": variation_p_values})
                    fig_var_pvals = go.Figure(go.Scatter(x=formatted_labels, y=variation_p_values, mode='lines+markers'))
                    if formatted_labels: fig_var_pvals.add_shape(type="line", x0=formatted_labels[0], x1=formatted_labels[-1], y0=0.05, y1=0.05, line=dict(color="red", dash="dash"))
                    fig_var_pvals.update_layout(title="P-Values (Degree Variation)", height=400); st.plotly_chart(fig_var_pvals, use_container_width=True)
                    st.session_state.report_elements.append({"type": "plot", "fig": fig_var_pvals, "title": "Degree Variation Significance (vs. Poisson)", "caption": caption_chi2_var, "part": 1})
                    st.session_state.report_elements.append({"type": "df_text", "df": df_variation_pvals_dl, "title": "Degree Variation P-Value Data", "part": 1})
                    dl_c1_sub_8, dl_c2_sub_8 = st.columns(2) 
                    with dl_c1_sub_8: st.download_button("Data (TSV)", convert_df_to_tsv(df_variation_pvals_dl), "degree_variation_significance.tsv", key="dl_chi2_var_tsv_tab1")
                    with dl_c2_sub_8: st.download_button("Plot (PNG)", convert_fig_to_png(fig_var_pvals), "degree_variation_plot.png", key="dl_chi2_var_png_tab1")
                else: st.warning("Need >1 network.")

            st.markdown("#### Network Complexity Analysis (Per Selected Network)")
            caption_complexity = "This section provides metrics for a selected individual network to help assess its complexity (e.g., random vs. scale-free characteristics)."
            st.caption(caption_complexity)
            if not all_graphs_nx: st.warning("No graph objects for complexity analysis.")
            else:
                graph_name_options = edge_df_names
                if graph_name_options:
                    selected_graph_name_complex = st.selectbox("Select a network for complexity analysis:", graph_name_options, key="complexity_graph_selector")
                    if selected_graph_name_complex:
                        try:
                            selected_index = edge_df_names.index(selected_graph_name_complex)
                            G_complex = all_graphs_nx[selected_index] 
                            st.markdown(f"##### Metrics for Network: `{selected_graph_name_complex}`")
                            col_metric1, col_metric2, col_metric3 = st.columns(3)
                            num_nodes = G_complex.number_of_nodes(); num_edges = G_complex.number_of_edges(); density = nx.density(G_complex); avg_clustering_graph = nx.average_clustering(G_complex)
                            complexity_metrics_data = {"Number of Nodes": num_nodes, "Number of Edges": num_edges, "Network Density": f"{density:.4f}", "Avg. Clustering Coeff.": f"{avg_clustering_graph:.4f}"}
                            col_metric1.metric(label="Nodes", value=num_nodes); col_metric2.metric(label="Edges", value=num_edges); col_metric3.metric(label="Density", value=f"{density:.4f}"); col_metric1.metric(label="Avg. Clustering", value=f"{avg_clustering_graph:.4f}")
                            path_len_text = "N/A"
                            if nx.is_connected(G_complex): avg_shortest_path = nx.average_shortest_path_length(G_complex); path_len_text = f"{avg_shortest_path:.4f}"; col_metric2.metric(label="Avg. Shortest Path", value=path_len_text)
                            else:
                                if num_nodes > 0 and G_complex.number_of_edges() > 0: 
                                    largest_cc_nodes = max(nx.connected_components(G_complex), key=len); subgraph = G_complex.subgraph(largest_cc_nodes)
                                    if subgraph.number_of_nodes() > 1 : avg_shortest_path = nx.average_shortest_path_length(subgraph); path_len_text = f"{avg_shortest_path:.4f} (LCC)"; col_metric2.metric(label="Avg. Shortest Path (LCC)", value=f"{avg_shortest_path:.4f}"); st.caption("*Avg. Shortest Path on LCC.*")
                                    else: col_metric2.metric(label="Avg. Shortest Path (LCC)", value="N/A")
                                else: col_metric2.metric(label="Avg. Shortest Path", value="N/A")
                            complexity_metrics_data["Avg. Shortest Path"] = path_len_text
                            st.session_state.report_elements.append({"type": "complexity_metrics", "graph_name": selected_graph_name_complex, "metrics": complexity_metrics_data, "part": 1})
                            
                            caption_deg_dist = f"Degree distribution for network '{selected_graph_name_complex}'. A straight line on log-log scale can indicate scale-free properties."
                            st.markdown("###### Degree Distribution"); st.caption(caption_deg_dist)
                            if num_nodes > 0:
                                degrees_current_graph = [d for n, d in G_complex.degree()]
                                if not degrees_current_graph: st.write("No degrees to plot.")
                                else:
                                    degree_counts = Counter(degrees_current_graph); deg, cnt = zip(*degree_counts.items())
                                    df_degree_dist_dl = pd.DataFrame({'Degree': deg, 'Count': cnt}).sort_values(by='Degree')
                                    
                                    fig_deg_dist_bar = go.Figure(); fig_deg_dist_bar.add_trace(go.Bar(x=deg, y=cnt))
                                    log_x_dd_bar = st.checkbox(f"Log X-axis (Bar)", key=f"logx_dd_bar_{selected_graph_name_complex.replace('.', '_')}")
                                    log_y_dd_bar = st.checkbox(f"Log Y-axis (Bar)", key=f"logy_dd_bar_{selected_graph_name_complex.replace('.', '_')}")
                                    xaxis_type_dd_bar = "log" if log_x_dd_bar else "linear"; yaxis_type_dd_bar = "log" if log_y_dd_bar else "linear"
                                    fig_deg_dist_bar.update_layout(title=f"Degree Dist. (Bar) for '{selected_graph_name_complex}'", height=400, xaxis_type=xaxis_type_dd_bar, yaxis_type=yaxis_type_dd_bar); 
                                    st.plotly_chart(fig_deg_dist_bar, use_container_width=True)
                                    st.session_state.report_elements.append({"type": "complexity_plot", "fig": fig_deg_dist_bar, "title": f"Degree Distribution (Bar) for '{selected_graph_name_complex}'", "caption": caption_deg_dist, "part": 1})
                                    
                                    st.markdown("###### Degree Distribution (Scatter with Trendline)");
                                    fig_deg_dist_scatter = px.scatter(df_degree_dist_dl, x='Degree', y='Count', trendline="ols", title=f"Degree Dist. (Scatter) for '{selected_graph_name_complex}'")
                                    log_x_dd_scatter = st.checkbox(f"Log X-axis (Scatter)", key=f"logx_dd_scatter_{selected_graph_name_complex.replace('.', '_')}")
                                    log_y_dd_scatter = st.checkbox(f"Log Y-axis (Scatter)", key=f"logy_dd_scatter_{selected_graph_name_complex.replace('.', '_')}")
                                    xaxis_type_dd_scatter = "log" if log_x_dd_scatter else "linear"; yaxis_type_dd_scatter = "log" if log_y_dd_scatter else "linear"
                                    fig_deg_dist_scatter.update_layout(height=400, xaxis_type=xaxis_type_dd_scatter, yaxis_type=yaxis_type_dd_scatter)
                                    st.plotly_chart(fig_deg_dist_scatter, use_container_width=True)
                                    st.session_state.report_elements.append({"type": "complexity_plot", "fig": fig_deg_dist_scatter, "title": f"Degree Distribution (Scatter) for '{selected_graph_name_complex}'", "caption": caption_deg_dist + " OLS trendline shown.", "part": 1})
                                    
                                    st.session_state.report_elements.append({"type": "df_text", "df": df_degree_dist_dl, "title": f"Degree_Distribution_Data_{selected_graph_name_complex}", "part": 1})
                                    dl_c1_sub_9, dl_c2_sub_9 = st.columns(2)
                                    with dl_c1_sub_9: st.download_button(f"Data (TSV)", convert_df_to_tsv(df_degree_dist_dl), f"degree_dist_{selected_graph_name_complex}.tsv", key=f"dl_deg_dist_tsv_sel_{selected_graph_name_complex.replace('.', '_')}")
                                    with dl_c2_sub_9: st.download_button(f"Scatter Plot (PNG)", convert_fig_to_png(fig_deg_dist_scatter), f"degree_dist_scatter_{selected_graph_name_complex}_plot.png", key=f"dl_deg_dist_png_sel_scatter_{selected_graph_name_complex.replace('.', '_')}")
                            else: st.write("Graph empty.")
                        except ValueError: st.error(f"Could not find selected graph '{selected_graph_name_complex}'.")
                        except IndexError: st.error(f"Index error for graph '{selected_graph_name_complex}'.")
                else: st.info("No networks for complexity analysis.")
            
            st.markdown("#### Aggregated Degree Distribution (All Networks)")
            caption_agg_deg_dist = "Shows the mean count of nodes (P(k)) for each degree (k) across all uploaded networks. The shaded area represents +/- one standard deviation."
            st.caption(caption_agg_deg_dist)
            all_degree_distributions_data = []
            for idx, G_iter in enumerate(all_graphs_nx):
                graph_name_iter = edge_df_names[idx]
                degrees_iter = [d for n, d in G_iter.degree()]
                if degrees_iter:
                    degree_counts_iter = Counter(degrees_iter)
                    for k_val, pk_val in degree_counts_iter.items():
                        all_degree_distributions_data.append({'k': k_val, 'P(k)': pk_val, 'graph_name': graph_name_iter})
            if all_degree_distributions_data:
                df_all_deg_dist = pd.DataFrame(all_degree_distributions_data)
                df_agg_deg_dist = df_all_deg_dist.groupby('k')['P(k)'].agg(['mean', 'std']).reset_index()
                df_agg_deg_dist['std'] = df_agg_deg_dist['std'].fillna(0) 
                fig_agg_deg_dist = go.Figure([
                    go.Scatter(name='Mean P(k)', x=df_agg_deg_dist['k'], y=df_agg_deg_dist['mean'], mode='lines+markers', line=dict(color='rgb(31, 119, 180)')),
                    go.Scatter(name='Upper Bound (Mean+STD)', x=df_agg_deg_dist['k'], y=df_agg_deg_dist['mean'] + df_agg_deg_dist['std'], mode='lines', marker=dict(color="#444"), line=dict(width=0), showlegend=False),
                    go.Scatter(name='Lower Bound (Mean-STD)', x=df_agg_deg_dist['k'], y=df_agg_deg_dist['mean'] - df_agg_deg_dist['std'], marker=dict(color="#444"), line=dict(width=0), mode='lines', fillcolor='rgba(68, 68, 68, 0.3)', fill='tonexty', showlegend=False)
                ])
                log_x_agg_dd = st.checkbox(f"Log X-axis (Aggregated)", key="logx_agg_dd")
                log_y_agg_dd = st.checkbox(f"Log Y-axis (Aggregated)", key="logy_agg_dd")
                xaxis_type_agg_dd = "log" if log_x_agg_dd else "linear"; yaxis_type_agg_dd = "log" if log_y_agg_dd else "linear"
                fig_agg_deg_dist.update_layout(title='Mean Degree Distribution (+/- STD) Across All Networks', xaxis_title='Degree (k)', yaxis_title='Mean Number of Nodes P(k)', height=500, xaxis_type=xaxis_type_agg_dd, yaxis_type=yaxis_type_agg_dd)
                st.plotly_chart(fig_agg_deg_dist, use_container_width=True)
                st.session_state.report_elements.append({"type": "plot", "fig": fig_agg_deg_dist, "title": "Aggregated Degree Distribution (All Networks)", "caption": caption_agg_deg_dist, "part": 1})
                st.session_state.report_elements.append({"type": "df_text", "df": df_agg_deg_dist, "title": "Aggregated_Degree_Distribution_Data", "part": 1})
                col1_agg_dd, col2_agg_dd = st.columns(2)
                with col1_agg_dd: st.download_button("Data (TSV)", convert_df_to_tsv(df_agg_deg_dist), "aggregated_degree_distribution.tsv", key="dl_agg_deg_dist_tsv")
                with col2_agg_dd: st.download_button("Plot (PNG)", convert_fig_to_png(fig_agg_deg_dist), "aggregated_degree_distribution_plot.png", key="dl_agg_deg_dist_png")
            else: st.warning("No data available for aggregated degree distribution plot.")
        elif valid_dfs_for_processing_tab1 and not all_nodes_set_tab1: st.warning("Dataframes loaded, but no nodes extracted.")
        elif not edge_dataframes: st.info("Upload .edges.txt files or use examples (sidebar) for Part 1.")
    else: st.info("Upload .edges.txt files or use examples (sidebar) for Part 1 analysis.")

with tab2: 
    st.subheader('Chemical Interaction Analysis')
    if edge_dataframes and valid_dfs_for_processing_tab1 and all_nodes_set_tab1 : 
        st.write('**Legend examples:** HBOND - Hydrogen bonds; SSBOND - Disulphide bridges; IONIC - Ionic bond; VDW - van der Waals; etc.')
        def format_node_id_for_interactions(node_id): parts = node_id.split(':'); return f'{parts[1]}:{parts[3]}' if len(parts) >= 4 else node_id 
        interaction_counts_dict = {}; has_interaction_data = False; all_interaction_nodes_set = set()
        for df_iter in edge_dataframes:
            if 'Interaction' in df_iter.columns and 'NodeId1' in df_iter.columns:
                has_interaction_data = True
                valid_node_ids = df_iter[df_iter['NodeId1'].astype(str).str.count(':') == 3]['NodeId1']
                all_interaction_nodes_set.update(valid_node_ids.unique())
                counts = df_iter[df_iter['NodeId1'].astype(str).str.count(':') == 3].groupby(['Interaction', 'NodeId1']).size().reset_index(name='Count')
                for interaction_type in counts['Interaction'].unique():
                    if interaction_type not in interaction_counts_dict: interaction_counts_dict[interaction_type] = []
                    interaction_counts_dict[interaction_type].append(counts[counts['Interaction'] == interaction_type])
        if not has_interaction_data: st.warning("No suitable files for Part 2.")
        elif not all_interaction_nodes_set: st.warning("Interaction data found, but no valid NodeIds extracted.")
        else:
            interaction_stats_counts_dict = {}
            if not all_interaction_nodes_set: st.warning("No interaction nodes for Part 2.")
            else:
                sorted_interaction_nodes = sorted(list(all_interaction_nodes_set), key=lambda x: int(x.split(':')[1]))
                for interaction_type, data_list in interaction_counts_dict.items():
                    if not data_list: continue
                    combined_data = pd.concat(data_list)
                    if 'NodeId1' not in combined_data.columns: continue
                    stats = combined_data.groupby('NodeId1')['Count'].agg(['mean', 'std']).reindex(sorted_interaction_nodes).fillna(0).reset_index()
                    stats['FormattedNodeId1'] = stats['NodeId1'].apply(format_node_id_for_interactions)
                    stats = stats.set_index('NodeId1').reindex(sorted_interaction_nodes).reset_index()
                    stats['FormattedNodeId1'] = stats['NodeId1'].apply(format_node_id_for_interactions) 
                    interaction_stats_counts_dict[interaction_type] = stats.fillna(0)
                def create_interaction_count_plot(category_name, stats_df):
                    figInt = go.Figure(go.Scatter(x=stats_df['FormattedNodeId1'], y=stats_df['mean'], error_y=dict(type='data', array=stats_df['std'].fillna(0)), mode='markers', marker=dict(size=10), name=category_name))
                    figInt.update_layout(title=f'{category_name} Counts', height=400); return figInt
                if not interaction_stats_counts_dict: st.warning("No interaction statistics calculated.")
                else:
                    num_interaction_cols = 2; interaction_categories = list(interaction_stats_counts_dict.keys())
                    for i in range(0, len(interaction_categories), num_interaction_cols):
                        cols_interaction = st.columns(num_interaction_cols)
                        for j in range(num_interaction_cols):
                            if i + j < len(interaction_categories):
                                category_name = interaction_categories[i+j]
                                stats_df_for_plot = interaction_stats_counts_dict[category_name]
                                with cols_interaction[j]:
                                    caption_interaction = f"Mean count of '{category_name}' interactions per residue, with standard deviation shown as error bars."
                                    st.markdown(f"###### {category_name} Interactions"); st.caption(caption_interaction)
                                    if stats_df_for_plot.empty or 'FormattedNodeId1' not in stats_df_for_plot.columns or 'mean' not in stats_df_for_plot.columns:
                                        st.write(f"No data for {category_name}."); continue
                                    figInt = create_interaction_count_plot(category_name, stats_df_for_plot)
                                    st.plotly_chart(figInt, use_container_width=True)
                                    st.session_state.report_elements.append({"type": "interaction_plot", "fig": figInt, "title": f"{category_name} Interaction Counts", "caption": caption_interaction, "part": 2})
                                    st.session_state.report_elements.append({"type": "df_text", "df": stats_df_for_plot, "title": f"{category_name}_Interaction_Data", "part": 2})
                                    dl_c1_sub, dl_c2_sub = st.columns(2); 
                                    with dl_c1_sub: st.download_button(f"Data (TSV)", convert_df_to_tsv(stats_df_for_plot), f"interaction_{category_name}_counts.tsv",key=f"dl_int_tsv_{category_name}_tab2")
                                    with dl_c2_sub: st.download_button(f"Plot (PNG)", convert_fig_to_png(figInt), f"interaction_{category_name}_plot.png", key=f"dl_int_png_{category_name}_tab2")
    else: st.info("Upload .edges.txt files (sidebar) with 'Interaction'/'NodeId1' columns for Part 2.")

with tab3: 
    st.subheader('AlphaMissense Integration')
    if am_data_source_message: st.success(am_data_source_message)
    if not am_df_final.empty:
        try:
            am_df_final['Mutation'] = am_df_final['Mutation'].astype(str)
            am_df_final['Position_Extract'] = am_df_final['Mutation'].str.extract(r'[A-Za-z](\d+)[A-Za-z]')[0]
            am_df_final.dropna(subset=['Position_Extract'], inplace=True)
            am_df_final['Position'] = am_df_final['Position_Extract'].astype(int)
            am_df_final['AM_Score'] = pd.to_numeric(am_df_final['AM_Score'], errors='coerce')
            am_df_final.dropna(subset=['Position', 'AM_Score'], inplace=True)
            if am_df_final.empty: st.warning("No valid AlphaMissense data after cleaning.")
            else:
                caption_am_class = "Counts of AlphaMissense classifications (pathogenic, benign, ambiguous) for mutations at each residue position."
                st.markdown("#### AlphaMissense Classification Counts"); st.caption(caption_am_class)
                classification_counts = am_df_final.groupby(['Position', 'AM_Classification']).size().unstack(fill_value=0)
                color_map_am = {'pathogenic': 'orange', 'benign': 'green', 'ambiguous': 'blue'} 
                figAM1 = go.Figure()
                plot_order = [col for col in color_map_am.keys() if col in classification_counts.columns]; [plot_order.append(col) for col in classification_counts.columns if col not in plot_order]
                for ct in plot_order:
                    if ct in classification_counts.columns: figAM1.add_trace(go.Bar(name=ct, x=classification_counts.index, y=classification_counts[ct], marker_color=color_map_am.get(ct, 'grey')))
                figAM1.update_layout(barmode='stack', title='AlphaMissense Classification Counts', height=500); st.plotly_chart(figAM1, use_container_width=True)
                st.session_state.report_elements.append({"type": "am_plot", "fig": figAM1, "title": "AlphaMissense Classification Counts", "caption": caption_am_class, "part": 3})
                st.session_state.report_elements.append({"type": "df_text", "df": classification_counts.reset_index(), "title": "AlphaMissense_Classification_Counts_Data", "part": 3})
                col1_am1, col2_am1 = st.columns(2); 
                with col1_am1:st.download_button("Data (TSV)", convert_df_to_tsv(classification_counts.reset_index()),"AM_classification_counts.tsv",key="dl_am_class_tsv_tab3")
                with col2_am1:st.download_button("Plot (PNG)",convert_fig_to_png(figAM1),"AM_classification_plot.png",key="dl_am_class_png_tab3")

                caption_am_score = "Distribution of AlphaMissense pathogenicity scores for mutations at each residue position. Higher scores indicate higher likelihood of pathogenicity."
                st.markdown("#### AlphaMissense Score Distribution"); st.caption(caption_am_score)
                last_pos_am = am_df_final['Position'].max()
                figAM2 = px.box(am_df_final, x='Position', y='AM_Score', title='Predicted AlphaMissense Scores per Position')
                figAM2.add_shape(type='rect',x0=0,x1=last_pos_am,y0=0,y1=0.34,fillcolor='green',opacity=0.1,line_width=0); figAM2.add_shape(type='rect',x0=0,x1=last_pos_am,y0=0.34,y1=0.564,fillcolor='grey',opacity=0.1,line_width=0); figAM2.add_shape(type='rect',x0=0,x1=last_pos_am,y0=0.564,y1=1.0,fillcolor='red',opacity=0.1,line_width=0)
                figAM2.update_layout(height=500); st.plotly_chart(figAM2, use_container_width=True)
                st.session_state.report_elements.append({"type": "am_plot", "fig": figAM2, "title": "AlphaMissense Score Distribution", "caption": caption_am_score, "part": 3})
                st.session_state.report_elements.append({"type": "df_text", "df": am_df_final[['ProteinID', 'Mutation', 'Position', 'AM_Score', 'AM_Classification']], "title": "AlphaMissense_Full_Data", "part": 3})
                col1_am2, col2_am2 = st.columns(2); 
                with col1_am2:st.download_button("Data (TSV)",convert_df_to_tsv(am_df_final[['ProteinID','Mutation','Position','AM_Score','AM_Classification']]),"AlphaMissense_data.tsv",key="dl_am_full_tsv_tab3")
                with col2_am2:st.download_button("Plot (PNG)",convert_fig_to_png(figAM2),"AlphaMissense_score_plot.png",key="dl_am_full_png_tab3")
        except Exception as e: st.error(f"Error processing AM file: {e}. Check format.")
    elif uniprot_id_input or uploaded_file_am or use_example_file_am : st.warning("Could not load/process AM data. Check input (sidebar).")
    else: st.info("Enter UniProt ID, upload AM file, or use example (sidebar) for Part 3.")

with tab_docs: # New Tab for Documentation
    st.subheader("Help & Documentation")
    try:
        with open("documentation.md", "r", encoding="utf-8") as f:
            doc_markdown = f.read()
        st.markdown(doc_markdown, unsafe_allow_html=True)
    except FileNotFoundError:
        st.error("`documentation.md` file not found. Please ensure it is in the same directory as the script.")
    except Exception as e:
        st.error(f"Could not load documentation: {e}")


if not edge_dataframes and not (uploaded_files or use_example_files):
     st.info("Welcome to SlytheRINs! Upload .edges.txt files or use examples (sidebar) to begin analysis.")

