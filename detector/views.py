from django.shortcuts import render
import requests
import topoly
import os
import uuid
import urllib3
import numpy as np

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

API_ENDPOINT = "https://api.esmatlas.com/foldSequence/v1/pdb/"
KNOT_THRESHOLD = 0.50
WINDOW_SIZE = 400
OVERLAP = 200

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
        alex_dist = topoly.alexander(pdb_path, closure=2, tries=1000)
        if not alex_dist: return '0_1', 0.0
        prob_unknotted = alex_dist.get('0_1', 0.0)
        prob_knotted = 1.0 - prob_unknotted
        dominant_knot = max(alex_dist, key=alex_dist.get)
        return dominant_knot, prob_knotted
    except: return None, 0.0

def index(request):
    context = {}
    if request.method == "POST":
        sequence = request.POST.get("sequence", "").replace("\n", "").replace(" ", "").strip()
        context['sequence'] = sequence
        context['sequence_length'] = len(sequence)

        if sequence:
            if len(sequence) > 400:
                # --- SLIDING WINDOW MODE ---
                step = WINDOW_SIZE - OVERLAP
                best_prob, best_qual = 0.0, 0.0
                best_k_type = "Unknotted"
                best_pdb = ""
                is_knotted = False
                processed_any = False

                for start in range(0, len(sequence), step):
                    end = min(start + WINDOW_SIZE, len(sequence))
                    if (end - start) < 50: break
                    
                    sub_seq = sequence[start:end]
                    run_id = str(uuid.uuid4())[:8]
                    temp_file = f"temp_{run_id}.pdb"
                    
                    try:
                        resp = requests.post(API_ENDPOINT, data=sub_seq, verify=False)
                        if resp.status_code == 200:
                            processed_any = True
                            with open(temp_file, "w") as f:
                                f.write(resp.text)
                            
                            qual = calculate_plddt(temp_file)
                            k_type, prob = run_topology_analysis(temp_file)
                            
                            # Keep the chunk with the highest knot probability to show on the UI
                            if prob >= best_prob:
                                best_prob = prob
                                best_qual = qual
                                best_k_type = k_type
                                best_pdb = resp.text
                                
                            if prob > KNOT_THRESHOLD:
                                is_knotted = True
                    except Exception:
                        pass
                    finally:
                        if os.path.exists(temp_file):
                            os.remove(temp_file)
                
                if processed_any:
                    context['quality'] = round(best_qual, 1)
                    context['knot_prob'] = round(best_prob * 100, 1)
                    context['knot_type'] = best_k_type if best_prob > KNOT_THRESHOLD else "0_1"
                    context['is_knotted'] = is_knotted
                    context['success'] = True
                    context['pdb_data'] = best_pdb
                else:
                    context['error'] = "Sliding window failed to process any chunks."

            else:
                # --- STANDARD MODE (<= 400) ---
                run_id = str(uuid.uuid4())[:8]
                temp_file = f"temp_{run_id}.pdb"
                try:
                    response = requests.post(API_ENDPOINT, data=sequence, verify=False)
                    if response.status_code == 200:
                        pdb_content = response.text
                        with open(temp_file, "w") as f:
                            f.write(pdb_content)
                        
                        qual = calculate_plddt(temp_file)
                        k_type, prob = run_topology_analysis(temp_file)
                        
                        context['quality'] = round(qual, 1)
                        context['knot_prob'] = round(prob * 100, 1)
                        context['knot_type'] = k_type if prob > KNOT_THRESHOLD else "0_1"
                        context['is_knotted'] = prob > KNOT_THRESHOLD
                        context['success'] = True
                        context['pdb_data'] = pdb_content
                    else:
                        context['error'] = f"ESMFold API Error: {response.status_code}"
                except Exception as e:
                    context['error'] = f"Processing Error: {str(e)}"
                finally:
                    if os.path.exists(temp_file):
                        os.remove(temp_file)
                        
    return render(request, 'detector/index.html', context)
