# Code Guidelines for PBTA_RNA Project

This document defines the universal development conventions for all notebooks and analyses in this project.
Any working agent MUST read and follow these guidelines before writing any code.

---

## 1. Python Environment (uv)

- **Tool:** `uv` (fast Python package manager, drop-in replacement for pip)
- **Python version:** >3.10
- **Setup steps (once, already done):**
  ```bash
  cd /home/alon/menow_home_ass
  uv venv .venv --python 3.11
  source .venv/bin/activate
  ```
- **Installing packages:** Use `uv pip install` not `pip install`
  ```bash
  uv pip install <package_name>
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
- All notebooks go under `notebooks/`. Generated outputs stay alongside their notebook.

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

## 3. Project Structure

```
/home/alon/menow_home_ass/
├── context/                ← All reference/instruction .md files for agents
│   ├── <analysis>_instruction.md
│   ├── <analysis>_plan.md
│   ├── code_guidelines.md
│   ├── assignment.md
│   └── referencess.md
├── notebooks/              ← All notebooks
│   ├── <analysis>.ipynb
│   └── <analysis>/
│       └── nb_cells/       ← Optional: per-step cell sources
├── analyses/               ← Empty dirs for future work
│   ├── statistical_analysis/
│   └── mrna_analysis/
├── .gitignore
├── requirements.txt
├── opencode.json
└── PBTA_RNA/               ← Raw data (gitignored .txt files)
```

## 4. Coding Standards

- **Style:** Write clean, readable Python with comments for non-obvious logic.
- **Pandas:** Use method chaining where appropriate; avoid deep nesting.
- **Plotting:** Prefer Plotly (interactive) over matplotlib/seaborn. All plots must have: title, axis labels, legend (if applicable).
- **Validation:** Every step in the notebook must include a validation cell that prints missing-value counts or data integrity checks.
- **Independence:** Each notebook step should be independently runnable (repeat imports as needed).
- **File paths:** Use absolute paths rooted at `/home/alon/menow_home_ass/` or relative paths from that directory.
- **Errors:** If a step fails, print a clear error message and continue — do not crash the notebook.

## 5. General Workflow Sequence

The agent should follow this order:

1. **Read context** — Read the relevant files from `context/`:
   - `<analysis>_instruction.md` (the step-by-step notebook plan)
   - `<analysis>_plan.md` (schema reference, if any)
   - `code_guidelines.md` (this file)
   - `referencess.md` (statistical methods, if relevant)
2. **Set up environment** — Activate existing venv; install any missing libraries via `uv pip install`
3. **Freeze dependencies** — `uv pip freeze > requirements.txt`
4. **Commit setup** — `chore: set up Python environment and freeze dependencies`
5. **Build the notebook** — Create `<analysis>.ipynb` in `notebooks/` following the instruction file
6. **Commit notebook** — `feat: implement <analysis> notebook`
7. **Verify** — Run the notebook end-to-end:
   ```bash
   jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=600 notebooks/<analysis>.ipynb --output notebooks/<analysis>_executed.ipynb
   ```
8. **Fix errors** if any, re-run until clean
9. **Commit fixes** — `fix: resolve notebook errors from verification run`
10. **Generate summary** — If the notebook produces a summary file, commit it too
