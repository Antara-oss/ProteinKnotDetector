import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import topoly
import os
import urllib3
import numpy as np

# --- Configuration ---
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
API_ENDPOINT = "https://api.esmatlas.com/foldSequence/v1/pdb/"
MINIMUM_PLDDT = 70.0
SAMPLING_ITERATIONS = 1000
KNOT_THRESHOLD = 0.50
WINDOW_SIZE = 400
OVERLAP = 200

# --- Page Setup ---
st.set_page_config(page_title="Protein Knot Detector", layout="wide")

st.title("Protein Knot Detector")
st.markdown("Automated Topological Analysis Pipeline")

# --- Sidebar (Technical Specs Only) ---
with st.sidebar:
    st.header("System Specifications")
    st.markdown("""
    **Structure Engine:** ESMFold API
    **Topology Algorithm:** Alexander Polynomial
    **Sampling:** 1000 Stochastic Closures
    """)
    st.divider()
    st.caption("Powered by Topoly & ESMAtlas")

# --- Helper Functions ---

def sanitize_sequence(sequence):
    sequence = sequence.upper()
    sequence = "".join(sequence.split())
    standard_codes = set("ACDEFGHIKLMNPQRSTVWY")
    clean_seq = "".join([char for char in sequence if char in standard_codes])
    return clean_seq

def calculate_plddt(pdb_path):
    plddt_values = []
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    try:
                        val = float(line[60:66].strip())
                        plddt_values.append(val)
                    except ValueError: continue
        if not plddt_values: return 0.0
        return np.mean(plddt_values)
    except Exception:
        return 0.0

def run_topology_analysis(pdb_path):
    try:
        alex_dist = topoly.alexander(pdb_path, closure=2, tries=SAMPLING_ITERATIONS)
        if not alex_dist: return '0_1', 0.0
        prob_unknotted = alex_dist.get('0_1', 0.0)
        prob_knotted = 1.0 - prob_unknotted
        dominant_knot = max(alex_dist, key=alex_dist.get)
        return dominant_knot, prob_knotted
    except Exception:
        return None, 0

def render_mol(pdb_file):
    with open(pdb_file, "r") as f:
        pdb_content = f.read()
    view = py3Dmol.view(width=800, height=500)
    view.addModel(pdb_content, "pdb")
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.zoomTo()
    showmol(view, height=500, width=800)

# --- Core Logic ---

sequence_input = st.text_area("Enter Sequence (FASTA format):", height=150)

if st.button("Run Analysis", type="primary"):
    if not sequence_input:
        st.error("Please enter a valid sequence.")
    else:
        clean_seq = sanitize_sequence(sequence_input)
        seq_len = len(clean_seq)
        
        # --- LARGE SEQUENCE LOGIC ---
        if seq_len > 400:
            st.info(f"Sequence Length: {seq_len} residues. Engaging Sliding Window Mode.")
            
            progress_bar = st.progress(0)
            status_text = st.empty()
            results = []
            
            step = WINDOW_SIZE - OVERLAP
            chunks = []
            for start in range(0, seq_len, step):
                end = min(start + WINDOW_SIZE, seq_len)
                if (end - start) < 50: break
                chunks.append((start, end))

            total_chunks = len(chunks)
            
            for i, (start, end) in enumerate(chunks):
                chunk_id = f"fragment_{start+1}_{end}"
                filename = f"{chunk_id}.pdb"
                sub_seq = clean_seq[start:end]
                
                status_text.text(f"Processing Fragment {i+1}/{total_chunks}...")
                
                try:
                    response = requests.post(API_ENDPOINT, data=sub_seq, verify=False)
                    if response.status_code == 200:
                        with open(filename, "w") as f:
                            f.write(response.text)
                        
                        plddt = calculate_plddt(filename)
                        knot_type, prob = run_topology_analysis(filename)
                        
                        if prob > KNOT_THRESHOLD:
                            results.append({
                                "Fragment": f"{start+1}-{end}",
                                "Type": knot_type,
                                "Probability": f"{prob*100:.1f}%",
                                "pLDDT": f"{plddt:.1f}",
                                "File": filename
                            })
                except Exception:
                    pass
                
                progress_bar.progress((i + 1) / total_chunks)

            status_text.text("Analysis Complete")
            
            if results:
                st.error(f"Knots Detected in {len(results)} fragments")
                st.table(results)
                
                selected_frag = st.selectbox("Visualize Fragment:", [r["File"] for r in results])
                if selected_frag:
                    render_mol(selected_frag)
            else:
                st.success("No knots detected.")

        # --- SMALL SEQUENCE LOGIC ---
        else:
            filename = "analysis_result.pdb"
            with st.spinner("Analyzing structure..."):
                try:
                    response = requests.post(API_ENDPOINT, data=clean_seq, verify=False)
                    if response.status_code == 200:
                        with open(filename, "w") as f:
                            f.write(response.text)
                        
                        # Layout: Report on Left, Structure on Right
                        col1, col2 = st.columns([1, 2])
                        
                        with col1:
                            st.subheader("Analysis Report")
                            
                            plddt = calculate_plddt(filename)
                            knot_type, prob = run_topology_analysis(filename)
                            
                            st.metric("Structure Quality (pLDDT)", f"{plddt:.1f}")
                            st.metric("Knot Probability", f"{prob*100:.1f}%")
                            
                            st.divider()
                            
                            if plddt < MINIMUM_PLDDT:
                                st.warning("Low Confidence Structure")
                            
                            if prob > KNOT_THRESHOLD:
                                st.error(f"Knot Detected: {knot_type}")
                            elif prob < 0.2:
                                st.success("Unknotted (0_1)")
                            else:
                                st.info("Ambiguous Topology")
                                
                        with col2:
                            st.subheader("3D Structure")
                            render_mol(filename)
                                
                    else:
                        st.error("API Connection Failed")
                except Exception as e:
                    st.error(f"Error: {e}")
