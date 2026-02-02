import requests
import topoly
import os
import sys
import urllib3
import numpy as np

# Suppress SSL warnings
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# --- CONFIGURATION PARAMETERS ---
API_ENDPOINT = "https://api.esmatlas.com/foldSequence/v1/pdb/"
MINIMUM_PLDDT = 70.0
SAMPLING_ITERATIONS = 1000
KNOT_THRESHOLD = 0.50

# --- SLIDING WINDOW SETTINGS (TUNED FOR FREE API) ---
WINDOW_SIZE = 400    # Lowered to 400 to avoid Error 413
OVERLAP = 200        # Keeps context between chunks

def retrieve_structure(sequence, filename):
    sequence = sequence.replace("\n", "").replace(" ", "").strip()
    try:
        response = requests.post(API_ENDPOINT, data=sequence, verify=False)
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
    except: return 0.0

def run_topology_analysis(pdb_path):
    try:
        alex_dist = topoly.alexander(pdb_path, closure=2, tries=SAMPLING_ITERATIONS)
        if not alex_dist: return '0_1', 0.0, 0.0
        
        prob_unknotted = alex_dist.get('0_1', 0.0)
        prob_knotted = 1.0 - prob_unknotted
        dominant_knot = max(alex_dist, key=alex_dist.get)
        
        return dominant_knot, 0, prob_knotted
    except: return None, 0, 0

def analyze_chunk(sequence, start_index, end_index):
    chunk_id = f"fragment_{start_index}_{end_index}"
    filename = f"{chunk_id}.pdb"
    
    print(f"\n[INFO] Processing Fragment: Residues {start_index} to {end_index}")
    
    pdb_file = retrieve_structure(sequence, filename)
    
    if pdb_file:
        plddt = calculate_plddt(pdb_file)
        knot_type, _, prob = run_topology_analysis(pdb_file)
        
        print(f"   -> Quality (pLDDT): {plddt:.1f}")
        print(f"   -> Knot Probability: {prob*100:.1f}%")
        
        # Return True only if knot is detected AND reliable
        if prob > KNOT_THRESHOLD:
            return True, knot_type, plddt, filename
        
    return False, None, 0, None

def run_sliding_window(full_sequence):
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
        is_knotted, k_type, qual, fname = analyze_chunk(sub_seq, start+1, end)
        
        if is_knotted:
            knots_found.append({
                "range": f"{start+1}-{end}",
                "type": k_type,
                "quality": qual
            })
            
    return knots_found

if __name__ == "__main__":
    print("\n--- PROTEIN TOPOLOGY DETECTION PIPELINE v2.1 (Free API Tuned) ---")
    raw_seq = input("Enter Sequence (FASTA): ").strip()
    
    # --- UPDATED DECISION LOGIC (Trigger at 400) ---
    if len(raw_seq) > 400:
        results = run_sliding_window(raw_seq)
        
        print("\n" + "="*50)
        print("FINAL ANALYSIS REPORT")
        print("-" * 50)
        if results:
            print(f"CONCLUSION: KNOTS DETECTED in {len(results)} fragments.")
            for res in results:
                print(f"Location: Res {res['range']} | Type: {res['type']} | Qual: {res['quality']:.1f}")
        else:
            print("CONCLUSION: No knots detected in any fragment.")
            
    else:
        fname = "analysis_result.pdb"
        print(f"[INFO] Sequence length ({len(raw_seq)}) within API limits.")
        if retrieve_structure(raw_seq, fname):
            qual = calculate_plddt(fname)
            k_type, _, prob = run_topology_analysis(fname)
            
            print("\n" + "="*50)
            print("FINAL ANALYSIS REPORT")
            print("-" * 50)
            print(f"QUALITY STATUS: {'PASS' if qual >= 70 else 'UNRELIABLE'} ({qual:.1f})")
            print(f"Knotted Probability: {prob*100:.1f}%")
            
            if prob > KNOT_THRESHOLD:
                print(f"CONCLUSION: KNOT DETECTED ({k_type})")
            else:
                print(f"CONCLUSION: UNKNOTTED")
    print("="*50 + "\n")
