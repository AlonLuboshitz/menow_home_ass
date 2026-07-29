import json
import re

with open('clinical_analysis.ipynb') as f:
    nb = json.load(f)

cells = nb['cells']

# Mapping of step number -> (markdown cell index, list of code cell indices)
# Step identifiers from markdown titles
step_info = {
    1:  "Load & Profile Patient Data",
    2:  "Patient Demographics",
    3:  "Patient Survival Overview",
    4:  "Cancer Predispositions -- Prevalence & Demographics",
    5:  "Load & Profile Sample Data",
    6:  "Sample Cancer Type Distributions",
    7:  "Tumor Purity & Ploidy",
    8:  "Molecular Subtype Landscape",
    9:  "Sequencing Strategy & RNA Library",
    10: "Merge Patient + Sample Data",
    11: "Samples per Patient",
    12: "Survival by Cancer Group",
    13: "Survival by Molecular Subtype (Global)",
    14: "Age at Diagnosis by Cancer Group",
    15: "Sex Balance by Cancer Group",
    16: "Purity by Cancer Group & Tumor Type",
    17: "Predisposition vs Outcome",
    18: "CNS Region vs Cancer Group",
    19: "Generate Summary Report",
    20: "Summary Table of All Figures",
    '9a': "Multi-Cancer-Group Analysis",
}

# Purposes from markdown cells
step_purpose = {
    1:  "Comprehensive overview of the patient dataset -- size, column types, missingness levels.",
    2:  "Age distribution, sex balance, race/ethnicity makeup.",
    3:  "Outcome overview and Kaplan-Meier survival curves.",
    4:  "How common each predisposition is, and demographic patterns across predispositions.",
    5:  "Overview of the sample-level dataset.",
    6:  "Histological and anatomical breakdown of all samples.",
    7:  "Distribution of tumor purity and ploidy across samples.",
    8:  "The diversity of molecular subtypes and their relationship to cancer groups.",
    9:  "What sequencing methods and library prep were used.",
    10: "How well the two datasets connect.",
    11: "How many patients have single vs. multiple samples.",
    12: "OS and EFS stratified by major cancer groups.",
    13: "OS and EFS stratified by molecular subtype.",
    14: "Whether different cancer groups occur at different ages.",
    15: "Whether certain cancer groups show sex bias.",
    16: "Tumor purity differences across cancer groups and clinical states.",
    17: "Whether known cancer predisposition affects survival or age of onset.",
    18: "Anatomical distribution patterns of different cancer types.",
    19: "A markdown summary of all findings.",
    20: "Consolidated overview of every figure produced.",
    '9a': "How many patients have samples in different cancer groups.",
}

# Map code cell indices to steps based on content analysis
code_to_step = {
    # Step 1: Load & Profile Patient Data
    23: 1, 24: 1, 25: 1,
    # Step 2: Patient Demographics
    26: 2, 27: 2, 28: 2,
    # Step 3: Patient Survival Overview
    29: 3, 30: 3, 31: 3, 32: 3,
    # Step 4: Cancer Predispositions
    33: 4, 34: 4, 35: 4, 36: 4,
    # Step 5: Load & Profile Sample Data
    37: 5, 38: 5,
    # Step 6: Sample Cancer Type Distributions
    39: 6, 40: 6,
    # Step 7: Tumor Purity & Ploidy
    41: 7, 42: 7,
    # Step 8: Molecular Subtype Landscape
    43: 8, 44: 8, 45: 8,
    # Step 9: Sequencing Strategy & RNA Library
    46: 9, 47: 9,
    # Step 10: Merge Patient + Sample Data
    48: 10, 49: 10,
    # Step 11: Samples per Patient
    50: 11, 51: 11,
    # Step 12: Survival by Cancer Group
    52: 12, 53: 12, 54: 12,
    # Step 13: Survival by Molecular Subtype (Global)
    55: 13, 56: 13, 57: 13,
    # Step 14: Age at Diagnosis by Cancer Group
    58: 14, 59: 14,
    # Step 15: Sex Balance by Cancer Group
    60: 15, 61: 15,
    # Step 16: Purity by Cancer Group & Tumor Type
    62: 16,
    # Step 17: Predisposition vs Outcome
    63: 17, 64: 17,
    # Step 18: CNS Region vs Cancer Group
    65: 18, 66: 18,
    # Step 19: Generate Summary Report
    67: 19,
    # Step 20: Summary Table of All Figures
    68: 20, 69: 20,
    # Step 9a: Multi-Cancer-Group Analysis
    70: '9a', 71: '9a', 72: '9a',
}

# Markdown cell indices per step
# Cell 2 = Step 1, Cell 3 = Step 2, ..., Cell 21 = Step 20, Cell 22 = Step 9a
md_indices = {i: idx for i, idx in enumerate(range(2, 23), start=1)}
md_indices['9a'] = 22

# Build new cell order
# Cell 0: Title
# Cell 1: Preamble code
new_cells = [cells[0], cells[1]]

# Track steps in order: 1-20 then 9a
step_order = list(range(1, 21)) + ['9a']

for step in step_order:
    # Add markdown cell for this step
    md_idx = md_indices[step]
    new_cells.append(cells[md_idx])

    # Add code cells for this step
    for cell_idx, step_num in code_to_step.items():
        if step_num == step:
            code_cell = cells[cell_idx]
            # Prepend step comment
            title = step_info[step]
            purpose = step_purpose[step]
            comment = f"# Step {step}: {title} -- {purpose}\n"
            # Remove existing wrong comment line if present
            src = ''.join(code_cell['source'])
            src_lines = src.split('\n')
            # Remove lines that are the old incorrect comment
            cleaned_lines = [l for l in src_lines if not l.startswith('# Step 9a:')]
            new_src = comment + '\n'.join(cleaned_lines)
            if not new_src.endswith('\n'):
                new_src += '\n'
            code_cell['source'] = [new_src]
            new_cells.append(code_cell)

# Replace title cell
new_cells[0]['source'] = [
    "# PBTA_RNA Clinical Data Analysis\n",
    "\n",
    "**Exploratory analysis of pediatric brain tumor clinical data.**\n",
    "\n",
    "Data: PBTA_RNA study\n"
]

nb['cells'] = new_cells

with open('clinical_analysis.ipynb', 'w') as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)

print(f"Before: 73 cells")
print(f"After: {len(new_cells)} cells")
print("JSON valid: running validation...")

import json as j
with open('clinical_analysis.ipynb') as f:
    j.load(f)
print("JSON valid: YES")
