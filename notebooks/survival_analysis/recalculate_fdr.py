import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

# Read data
surv = pd.read_csv('/home/alon/menow_home_ass/notebooks/survival_analysis/survival_analysis_results.csv')
decile = pd.read_csv('/home/alon/menow_home_ass/notebooks/age_deciles/age_deciles_results.csv')

# Only per-group tests (exclude 'global')
surv_group = surv[surv['Group'] != 'global'].copy()

# Get all unique cancer groups that appear in either file
all_groups = sorted(set(surv_group['Group'].unique()) | set(decile['Group'].unique()))
print(f"Total cancer groups: {len(all_groups)}")

# For each group, collect per-group survival tests + decile tests
# Apply FDR within the group

rows = []
for group in all_groups:
    # Get survival tests for this group
    surv_g = surv_group[surv_group['Group'] == group].copy()
    
    # Get decile tests for this group
    decile_g = decile[decile['Group'] == group].copy()
    decile_g = decile_g.copy()
    # Add a Phase column to decile data for consistency
    decile_g['Phase'] = decile_g['Phase'].apply(lambda x: f'DECILE_{x}')
    
    # Combine into a single family for this group
    # We'll include all per-group KM log-rank tests + decile tests
    all_tests = pd.concat([surv_g, decile_g], ignore_index=True)
    
    if len(all_tests) == 0:
        continue
    
    # Get p-values, drop NaN
    p_vals = all_tests['p_value'].values
    valid_mask = pd.notna(p_vals)
    p_vals = p_vals[valid_mask]
    all_tests = all_tests[valid_mask].copy()
    
    if len(p_vals) == 0:
        continue
    
    # Apply FDR within this group
    _, q_vals, _, _ = multipletests(p_vals, method='fdr_bh')
    all_tests['Group_FDR'] = q_vals
    
    # Find FDR-significant tests
    sig = all_tests[all_tests['Group_FDR'] < 0.05]
    
    if len(sig) > 0:
        # Get N for this group - use the max N from any test
        n_val = int(all_tests['N'].max())
        
        # Build short test labels
        test_labels = []
        for _, r in sig.iterrows():
            comp = r.get('Comparison', str(r.get('Phase', '')))
            # Handle decile tests
            if 'DECILE' in str(r.get('Phase', '')):
                label = f"AGE_DECILE_{r['Phase'].replace('DECILE_', '')}"
            else:
                comp_str = str(comp)
                if 'AGE' in comp_str and 'EFS' in comp_str:
                    label = 'AGE × EFS'
                elif 'AGE' in comp_str and 'OS' in comp_str:
                    label = 'AGE × OS'
                elif 'TF' in comp_str and 'OS' in comp_str:
                    label = 'TF × OS'
                elif 'TF' in comp_str and 'EFS' in comp_str:
                    label = 'TF × EFS'
                elif 'TP' in comp_str and 'OS' in comp_str:
                    label = 'TP × OS'
                elif 'TP' in comp_str and 'EFS' in comp_str:
                    label = 'TP × EFS'
                elif 'SEX' in comp_str and 'OS' in comp_str:
                    label = 'SEX × OS'
                elif 'SEX' in comp_str and 'EFS' in comp_str:
                    label = 'SEX × EFS'
                elif 'PREDISPOSITION' in comp_str or 'Predisposition' in comp_str:
                    label = 'PRED'
                else:
                    label = comp_str[:30]
            test_labels.append(label)
        
        # Remove duplicates while preserving order
        seen = set()
        unique_labels = []
        for label in test_labels:
            if label not in seen:
                seen.add(label)
                unique_labels.append(label)
        
        rows.append({
            'Cancer_Group': group,
            'N_samples': n_val,
            'N_FDR_Significant': len(unique_labels),
            'Significant_Tests': ', '.join(unique_labels)
        })
        
        print(f"✅ {group:50s} N={n_val:<6d}  {len(unique_labels)} tests: {', '.join(unique_labels)}")

# Build summary
summary = pd.DataFrame(rows)
summary = summary.sort_values(['N_FDR_Significant', 'N_samples'], ascending=[False, False])

print(f"\n\n{'='*100}")
print(f"FDR-SIGNIFICANT RESULTS (FDR applied WITHIN each cancer group)")
print(f"{'='*100}")
print(f"Groups with at least one FDR-significant test: {len(summary)}")
print(f"\n{'Cancer Group':50s} {'N':>8s} {'Sig':>4s}  Tests")
print('-'*100)
for _, r in summary.iterrows():
    print(f"{r['Cancer_Group'][:49]:50s} {r['N_samples']:>8d} {r['N_FDR_Significant']:>4d}  {r['Significant_Tests'][:50]}")

# Save
out_dir = '/home/alon/menow_home_ass/notebooks/survival_analysis'
summary.to_csv(f'{out_dir}/fdr_significant_summary.csv', index=False)
print(f"\nSaved: {out_dir}/fdr_significant_summary.csv")

# Compare with old approach
print(f"\n\n=== COMPARISON: Old cross-group FDR vs New within-group FDR ===")
print(f"Old (cross-group FDR): 7 groups, 11 tests")
print(f"New (within-group FDR): {len(summary)} groups, {summary['N_FDR_Significant'].sum()} tests")
