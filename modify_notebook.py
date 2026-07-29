import json
import re

NOTEBOOK = "/home/alon/menow_home_ass/clinical_analysis.ipynb"

with open(NOTEBOOK) as f:
    nb = json.load(f)

cells = nb["cells"]

# --- Task 1: Simplify title cell ---
title_cell = cells[0]
title_cell["source"] = [
    "# PBTA_RNA Clinical Data Analysis\n",
    "\n",
    "**Exploratory analysis of pediatric brain tumor clinical data.**\n",
    "\n",
    "Data: PBTA_RNA study\n",
]

# --- Task 2: Add step-identifying comments to code cells ---
step_pattern = re.compile(r"^## Step (\d+\w?): (.+)", re.MULTILINE)
purpose_pattern = re.compile(r"^\*\*Purpose:\*\* (.+)", re.MULTILINE)

modified_count = 0
active_step = None

for cell in cells:
    if cell["cell_type"] == "markdown":
        src_text = "".join(cell["source"])
        step_match = step_pattern.search(src_text)
        if step_match:
            purpose_match = purpose_pattern.search(src_text)
            if purpose_match:
                step_num = step_match.group(1)
                step_title = step_match.group(2).strip()
                purpose = purpose_match.group(1).strip()
                active_step = (step_num, step_title, purpose)
    elif cell["cell_type"] == "code" and active_step is not None:
        step_num, step_title, purpose = active_step
        comment = f"# Step {step_num}: {step_title} — {purpose}\n"
        source = cell["source"]
        if source and source[0] == comment:
            continue
        source.insert(0, comment)
        modified_count += 1

with open(NOTEBOOK, "w") as f:
    json.dump(nb, f, indent=1)
    f.write("\n")

# --- Verification ---
with open(NOTEBOOK) as f:
    json.load(f)

print(f"Modified {modified_count} code cells")
print("Valid JSON: yes")

# Print first 3 commented cells
verify_count = 0
with open(NOTEBOOK) as f:
    nb2 = json.load(f)
for cell in nb2["cells"]:
    if cell["cell_type"] == "code" and cell["source"] and cell["source"][0].startswith("# Step"):
        print(f"\n  [{verify_count+1}] {cell['source'][0].rstrip()}")
        print(f"      Next line: {cell['source'][1][:60].rstrip() if len(cell['source']) > 1 else '(empty)'}")
        verify_count += 1
        if verify_count >= 3:
            break
