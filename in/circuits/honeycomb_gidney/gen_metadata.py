import json
import os
import re
import stim

dir_path = os.path.dirname(os.path.abspath(__file__))

for f in sorted(os.listdir(dir_path)):
    if not f.endswith('.stim'):
        continue
    
    # Parse filename params
    params = dict(re.findall(r'(\w+)=([^,]+)', f.replace('.stim', '')))
    
    c = stim.Circuit.from_file(os.path.join(dir_path, f))
    
    w = int(params['w'])
    h = int(params['h'])
    r = int(params['r'])
    
    metadata = {
        "code_family": "honeycomb",
        "code_name": "honeycomb_gidney",
        "source": "https://quantum-journal.org/papers/q-2022-09-21-813/",
        "noise_model": params["model"],
        "p": float(params["p"]),
        "w": w,
        "h": h,
        "rounds": r,
        "observable": params["obs"],
        "sheared": params["sheared"] == "True",
        "expected_distance": min(w//2, h//3),
        "k": c.num_observables,
        "num_qubits": c.num_qubits,
        "num_measurements": c.num_measurements,
        "num_detectors": c.num_detectors,
        "num_observables": c.num_observables,
    }
    
    json_path = os.path.join(dir_path, f.replace('.stim', '.json'))
    with open(json_path, 'w') as jf:
        json.dump(metadata, jf, indent=2)
    print(f"Wrote {os.path.basename(json_path)}")

