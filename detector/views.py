import tempfile
import uuid
from pathlib import Path

import numpy as np
import requests
import topoly
import urllib3
from django.shortcuts import render

# Suppress SSL warnings
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# --- CONFIGURATION PARAMETERS ---
API_ENDPOINT = "https://api.esmatlas.com/foldSequence/v1/pdb/"
MINIMUM_PLDDT = 70.0
SAMPLING_ITERATIONS = 1000
KNOT_THRESHOLD = 0.50
REQUEST_TIMEOUT = 120
STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

# --- SLIDING WINDOW SETTINGS (TUNED FOR FREE API) ---
WINDOW_SIZE = 400    # Lowered to 400 to avoid Error 413
OVERLAP = 200        # Keeps context between chunks


def sanitize_sequence(raw_sequence):
    lines = []
    for line in raw_sequence.splitlines():
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        lines.append(line)

    sequence = "".join(lines).upper()
    return "".join(char for char in sequence if char in STANDARD_AMINO_ACIDS)

def retrieve_structure(sequence, filename):
    sequence = sequence.replace("\n", "").replace(" ", "").strip()
    try:
        response = requests.post(
            API_ENDPOINT,
            data=sequence,
            verify=False,
            timeout=REQUEST_TIMEOUT,
        )
        if response.status_code == 200:
            with open(filename, "w") as f:
                f.write(response.text)
            return filename
        elif response.status_code == 413:
            print(f"[ERROR] Sequence too long for API (Status 413).")
            return None
        else:
            print(f"[ERROR] API Status: {response.status_code}")
            return None
    except Exception as e:
        print(f"[ERROR] Connection failed: {e}")
        return None

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
        mean = np.mean(plddt_values)
        return mean * 100 if mean <= 1.0 else mean
    except:
        return 0.0

def run_topology_analysis(pdb_path):
    try:
        alex_dist = topoly.alexander(pdb_path, closure=2, tries=SAMPLING_ITERATIONS)
        if not alex_dist: return '0_1', 0.0, 0.0
        
        prob_unknotted = alex_dist.get('0_1', 0.0)
        prob_knotted = 1.0 - prob_unknotted
        dominant_knot = max(alex_dist, key=alex_dist.get)
        return dominant_knot, 0, prob_knotted
    except:
        return None, 0, 0

def analyze_chunk(sequence, start_index, end_index, work_dir):
    chunk_id = f"fragment_{start_index}_{end_index}"
    filename = Path(work_dir) / f"{chunk_id}.pdb"
    
    print(f"\n[INFO] Processing Fragment: Residues {start_index} to {end_index}")
    
    pdb_file = retrieve_structure(sequence, filename)
    
    if pdb_file:
        plddt = calculate_plddt(pdb_file)
        knot_type, _, prob = run_topology_analysis(pdb_file)
        
        print(f"   -> Quality (pLDDT): {plddt:.1f}")
        print(f"   -> Knot Probability: {prob*100:.1f}%")
        
        # Return True only if knot is detected AND reliable
        if prob > KNOT_THRESHOLD:
            return True, knot_type, plddt
        
    return False, None, 0

def run_sliding_window(full_sequence, work_dir):
    seq_len = len(full_sequence)
    print(f"[INFO] Large sequence detected ({seq_len} > 400).")
    print(f"[INFO] Engaging Sliding Window Mode (Window={WINDOW_SIZE}, Overlap={OVERLAP})")
    
    knots_found = []
    
    # Loop through sequence with step size
    step = WINDOW_SIZE - OVERLAP
    for start in range(0, seq_len, step):
        end = min(start + WINDOW_SIZE, seq_len)
        
        # Don't process tiny chunks (<50 residues)
        if (end - start) < 50: break
            
        sub_seq = full_sequence[start:end]
        is_knotted, k_type, qual = analyze_chunk(sub_seq, start+1, end, work_dir)
        
        if is_knotted:
            knots_found.append({
                "range": f"{start+1}-{end}",
                "type": k_type,
                "quality": qual
            })
            
    return knots_found

if __name__ == "__main__":
    pass # Bypassed the terminal input so the web server doesn't crash


# ==========================================
# DJANGO WEB BRIDGE
# ==========================================
def index(request):
    context = {}
    if request.method == "POST":
        raw_sequence = request.POST.get("sequence", "")
        sequence = sanitize_sequence(raw_sequence)

        context['sequence'] = raw_sequence
        context['sequence_length'] = len(sequence)

        if raw_sequence and not sequence:
            context['error'] = "No valid amino acid sequence was found in the FASTA input."
        elif sequence:
            with tempfile.TemporaryDirectory(prefix="protein-knot-") as temp_dir:
                if len(sequence) > 400:
                    results = run_sliding_window(sequence, temp_dir)
                    
                    if results:
                        best = results[0]
                        context['quality'] = round(best['quality'], 1)
                        context['knot_prob_display'] = "Detected in fragment"
                        context['knot_prob_is_percent'] = False
                        context['knot_type'] = best['type']
                        context['is_knotted'] = True
                        context['success'] = True
                    else:
                        context['error'] = "No knots detected in any fragment."
                else:
                    run_id = str(uuid.uuid4())[:8]
                    fname = Path(temp_dir) / f"web_result_{run_id}.pdb"
                    pdb_file = retrieve_structure(sequence, fname)
                    
                    if pdb_file:
                        qual = calculate_plddt(pdb_file)
                        k_type, _, prob = run_topology_analysis(pdb_file)
                        
                        context['quality'] = round(qual, 1)
                        context['knot_prob_display'] = round(prob * 100, 1)
                        context['knot_prob_is_percent'] = True
                        context['knot_type'] = k_type if prob > KNOT_THRESHOLD else "Unknotted"
                        context['is_knotted'] = prob > KNOT_THRESHOLD
                        context['success'] = True
                        
                        try:
                            with open(pdb_file, "r") as f:
                                context['pdb_data'] = f.read()
                        except OSError:
                            pass
                    else:
                        context['error'] = "ESMFold API Error or Connection Failed."
                    
    return render(request, 'detector/index.html', context)
