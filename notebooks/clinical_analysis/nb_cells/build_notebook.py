import json, os, glob

CELLS_DIR = "/home/alon/menow_home_ass/nb_cells"

def md(s):
    return {"cell_type":"markdown","metadata":{},"source":[s]}

def code(lines):
    return {"cell_type":"code","metadata":{},"source":[l+"\n" for l in lines],"outputs":[],"execution_count":None}

def load_code(filename):
    path = os.path.join(CELLS_DIR, filename)
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f.readlines()]
    return code(lines)

cells = []

# Title
cells.append(md("# PBTA_RNA Clinical Data Analysis\n\n**Exploratory analysis of pediatric brain tumor clinical data.**\n\nThis notebook covers:\n- Patient demographics and survival outcomes\n- Cancer predispositions\n- Sample-level tumor annotations\n- Cross-dataset integration\n- Survival analysis by cancer group and molecular subtype\n- Tumor purity and ploidy analysis\n\nData: PBTA_RNA study"))

# Imports + helpers (single file)
cells.append(load_code("imports.py"))

# Steps
for step in sorted(glob.glob(os.path.join(CELLS_DIR, "step_*.md"))):
    with open(step) as f:
        cells.append(md(f.read()))

for step in sorted(glob.glob(os.path.join(CELLS_DIR, "step_*.py"))):
    cells.append(load_code(os.path.basename(step)))

# Write notebook
notebook = {"nbformat":4,"nbformat_minor":5,
    "metadata":{"kernelspec":{"display_name":"Python (PBTA)","language":"python","name":"pbta_env"},
                "language_info":{"name":"python","version":"3.11.15"}},
    "cells":cells}

with open('/home/alon/menow_home_ass/clinical_analysis.ipynb','w') as f:
    json.dump(notebook, f, indent=1)

n_md = sum(1 for c in cells if c['cell_type']=='markdown')
n_code = sum(1 for c in cells if c['cell_type']=='code')
print(f"Notebook: {len(cells)} cells ({n_md} md + {n_code} code)")
