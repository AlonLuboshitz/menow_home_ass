# Code Guidelines for Building `clinical_analysis.ipynb`

This document defines the development conventions for the PBTA_RNA clinical analysis notebook.  
The working agent MUST read and follow these guidelines before writing any code.

---

## 1. Python Environment (uv)

- **Tool:** `uv` (fast Python package manager, drop-in replacement for pip)
- **Python version:** >3.10
- **Setup steps (once):**
  ```bash
  cd /home/alon/menow_home_ass
  uv venv .venv --python 3.11
  source .venv/bin/activate
  ```
- **Installing packages:** Use `uv pip install` not `pip install`
  ```bash
  uv pip install pandas numpy plotly scipy jupyter
  ```
- **Freezing dependencies (after all packages installed):**
  ```bash
  uv pip freeze > requirements.txt
  ```
- **Notebook kernel:** After setting up the venv, register it as a Jupyter kernel:
  ```bash
  uv pip install ipykernel
  python -m ipykernel install --user --name=pbta_env --display-name="Python (PBTA)"
  ```
- **Never** use `!pip install` inside the notebook — all dependencies must be pre-installed before running.
- All generated output files (`.ipynb`, `.md`, `.txt`) belong in the project root `/home/alon/menow_home_ass/`.

## 2. Version Control (Git)

- **Commit convention:** Conventional Commits (`feat:`, `fix:`, `chore:`, `docs:`, `refactor:`)
- **One commit per logical change** — not one giant commit with everything.
- **Commit message format:**
  ```
  <type>: <short description>

  <optional body>
  ```
  Examples:
  - `chore: add code guidelines and environment setup`
  - `feat: implement clinical analysis notebook`
  - `docs: add summary report of findings`
- **Before each commit**, run `git status` to verify only intended files are staged.
- **Do NOT commit**:
  - `.venv/`, `__pycache__/`, `.ipynb_checkpoints/`
  - `PBTA_RNA/` (data files) and `PBTA_RNA.zip`
  - `.opencode/` (opencode configuration)
- Always commit from `/home/alon/menow_home_ass/`.
- Use `git add <file>` for specific files, never `git add .` blindly.

## 3. Coding Standards

- **Style:** Write clean, readable Python with comments for non-obvious logic.
- **Pandas:** Use method chaining where appropriate; avoid deep nesting.
- **Plotting:** Prefer Plotly (interactive) over matplotlib/seaborn. All plots must have: title, axis labels, legend (if applicable).
- **Validation:** Every step in the notebook must include a validation cell that prints missing-value counts or data integrity checks.
- **Independence:** Each notebook step should be independently runnable (repeat imports as needed).
- **File paths:** Always use absolute paths rooted at `/home/alon/menow_home_ass/` or relative paths from that directory.
- **Errors:** If a step fails, print a clear error message and continue — do not crash the notebook.

## 4. Workflow Sequence

The agent should follow this order:

1. Set up uv environment (if not already done)
2. Read the following reference files:
   - `/home/alon/menow_home_ass/instruction.md` (the 20-step notebook plan)
   - `/home/alon/menow_home_ass/basic_clinical_analysis_plan.md` (schema reference)
   - `/home/alon/menow_home_ass/referencess.md` (statistical methods)
3. Install all required libraries via `uv pip install`
4. Freeze dependencies: `uv pip freeze > requirements.txt`
5. Commit environment setup: `chore: set up Python environment and freeze dependencies`
6. Build `clinical_analysis.ipynb` following instruction.md
7. Commit the notebook: `feat: implement clinical analysis notebook with 20 steps`
8. Run the notebook end-to-end to verify no errors
9. If errors occur, fix and re-run until clean
10. Commit fixes: `fix: resolve notebook errors from verification run`
11. Build the summary: the notebook's Step 19 generates `basic_clinical_summary.md` automatically
