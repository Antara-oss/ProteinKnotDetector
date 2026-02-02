import requests
import topoly
import os
import sys
import urllib3
import numpy as np

# Suppress SSL warnings for cleaner terminal output
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# --- CONFIGURATION PARAMETERS ---
API_ENDPOINT = "https://api.esmatlas.com/foldSequence/v1/pdb/"
MINIMUM_PLDDT = 70.0        # Quality cutoff for reliable structures
SAMPLING_ITERATIONS = 1000  # Number of stochastic closures for Topoly
KNOT_THRESHOLD = 0.50       # Probability threshold to classify as knotted

def retrieve_structure(sequence, job_id):
    # NOTE: Add comment here about fetching data from ESMFold
    
    filename = f"{job_id}.pdb"
    # sanitize sequence input
    sequence = sequence.replace("\n", "").replace(" ", "").strip()
    
    print(f"[INFO] Initiating API request for ID: {job_id}")
    
    # Limit set to 2000 residues
    if len(sequence) > 2000:
        print("[ERROR] Sequence length exceeds API limit (2000 residues).")
        return None

    try:
        response = requests.post(API_ENDPOINT, data=sequence, verify=False)
        
        if response.status_code == 200:
            with open(filename, "w") as f:
                f.write(response.text)
            print(f"[INFO] Structure successfully saved to {filename}")
            return filename
        else:
            print(f"[ERROR] API returned status code: {response.status_code}")
            return None
            
    except Exception as e:
        print(f"[ERROR] Connection failure: {e}")
        return None

def calculate_plddt(pdb_path):
    # NOTE: Add comment here about extracting B-factors for quality check
    
    print(f"[INFO] Assessing structural quality (pLDDT)...")
    plddt_values = []
    
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    try:
                        # Extract B-factor columns (60-66)
                        val = float(line[60:66].strip())
                        plddt_values.append(val)
                    except ValueError:
                        continue
        
        if not plddt_values:
            return 100.0
        
        mean_plddt = np.mean(plddt_values)
        
        # Normalization to 0-100 scale if necessary
        if mean_plddt <= 1.0:
            mean_plddt *= 100
            
        print(f"[INFO] Mean pLDDT Score: {mean_plddt:.2f}")
        return mean_plddt
        
    except Exception as e:
        print(f"[WARNING] Could not calculate pLDDT: {e}")
        return 0.0

def run_topology_analysis(pdb_path):
    # NOTE: Add comment here about the Alexander/Jones polynomial logic
    
    print(f"[INFO] Running stochastic topology analysis (n={SAMPLING_ITERATIONS})...")
    
    try:
        # Primary Analysis: Alexander Polynomial
        alex_dist = topoly.alexander(pdb_path, closure=2, tries=SAMPLING_ITERATIONS)
        
        # Secondary Validation: Jones Polynomial
        jones_dist = topoly.jones(pdb_path, closure=2, tries=SAMPLING_ITERATIONS)
        
        if not alex_dist:
            return '0_1', 0.0, 0.0

        # Calculate probabilities
        prob_unknotted = alex_dist.get('0_1', 0.0)
        prob_knotted = 1.0 - prob_unknotted
        
        # Identify dominant topology
        dominant_knot = max(alex_dist, key=alex_dist.get)
        dominant_conf = alex_dist[dominant_knot]
        
        return dominant_knot, dominant_conf, prob_knotted

    except Exception as e:
        print(f"[ERROR] Topology analysis failed: {e}")
        return None, 0, 0

if __name__ == "__main__":
    print("\n--- PROTEIN TOPOLOGY DETECTION PIPELINE ---")
    
    # Direct Input
    raw_seq = input("Enter Sequence (FASTA): ").strip()
    job_name = "analysis_result"
    
    target_file = retrieve_structure(raw_seq, job_name)
    
    # Execution Block
    if target_file:
        quality_score = calculate_plddt(target_file)
        knot_type, conf_score, knot_prob = run_topology_analysis(target_file)
        
        print("\n" + "="*50)
        print("FINAL ANALYSIS REPORT")
        print("-" * 50)
        
        # Quality Verdict
        if quality_score < MINIMUM_PLDDT:
            print(f"QUALITY STATUS: UNRELIABLE ({quality_score:.1f})")
            print("Warning: Structure low confidence may yield false positives.")
        else:
            print(f"QUALITY STATUS: PASS ({quality_score:.1f})")

        print("-" * 50)

        # Topology Verdict
        print(f"Iterations: {SAMPLING_ITERATIONS}")
        print(f"Unknotted Probability: {(1.0-knot_prob)*100:.1f}%")
        print(f"Knotted Probability:   {knot_prob*100:.1f}%")
        
        print("-" * 50)
        
        if knot_prob > KNOT_THRESHOLD:
            print(f"CONCLUSION: KNOT DETECTED")
            print(f"Type: {knot_type}")
        elif knot_prob < 0.2:
            print(f"CONCLUSION: UNKNOTTED")
        else:
            print(f"CONCLUSION: AMBIGUOUS / SHALLOW KNOT")
            
        print("="*50 + "\n")
