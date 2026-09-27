#!/usr/bin/env python3
"""D4b discordance audit (ADDENDUM_4, descriptive): the cardiac-passing
candidates excluded from the both-tissue core - gate-(a) serum directions
x gate-(b) blood results. Reads committed artifacts only."""
import csv, json

CORE6 = {"miR-1-3p","miR-122-5p","miR-192-5p","miR-30c-5p","miR-145-5p","miR-194-5p"}

gate_a = {}
with open("results/h2_gate_a_passing.csv") as f:
    for r in csv.DictReader(f):
        gate_a[r["mirna"].replace("hsa-","")] = r

cardiac, blood = {}, {}
with open("results/h2_gate_b_enrichment.csv") as f:
    for r in csv.DictReader(f):
        if r["label"].endswith("|hipsc_GSE203525"):
            cardiac[r["mirna"]] = r
        elif r["label"].endswith("|blood_GSE244827"):
            blood[r["mirna"]] = r

rows = []
for m, c in sorted(cardiac.items(), key=lambda kv: float(kv[1]["fdr"])):
    if float(c["fdr"]) > 0.05 or m in CORE6:
        continue
    a = gate_a.get(m, {})
    b = blood.get(m, {})
    d = float(a.get("d_sev_mild", "nan"))
    rows.append({
        "mirna": m,
        "serum_direction_sev_vs_mild": "down" if d < 0 else "up",
        "cliff_d": round(d, 3),
        "cardiac_frac_DE_pred_dir": round(float(c["frac_DE"]), 4),
        "cardiac_fdr": float(c["fdr"]),
        "blood_frac_DE_pred_dir": round(float(b.get("frac_DE", "nan")), 4),
        "blood_fdr": float(b.get("fdr", "nan")),
        "blood_passes": str(float(b.get("fdr", "1")) <= 0.05),
    })

with open("results/d4b_discordance_audit.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print(json.dumps({"n_excluded_cardiac_passers": len(rows),
                  "members": [r["mirna"] for r in rows]}, indent=1))
