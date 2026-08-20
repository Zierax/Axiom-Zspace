import os, subprocess
HERE = "/mnt/d/Axioms/Axiom-Zspace-CODE/C99-Version/tests"
ROOT = "/mnt/d/Axioms/Axiom-Zspace-CODE/C99-Version"
BIN = os.path.join(ROOT, "bin", "zspace_card")

lines = {
    "period_days": 10.0, "transit_depth": 0.01, "transit_duration_hrs": 3.0,
    "t0_days": 1.0, "stellar_mass_solar": 1.0, "stellar_radius_solar": 1.0,
    "stellar_teff_k": 5778, "stellar_logg": 4.44, "planet_radius_earth": 10.0,
    "bls_snr": 15.0, "bls_fap": 0.0001, "even_odd_delta_sigma": 0.5,
    "shape_ratio": 1.0, "secondary_snr": 1.0, "secondary_depth_ratio": 0.1,
    "alias_secondary_ratio": 0.0, "coherent_evidence": 0,
    "centroid_sigma": 0.0, "limb_dark_u1": 0.45, "limb_dark_u2": 0.15,
    "s_periodicity": 0.85, "s_depth": 0.85, "s_limb": 0.85, "s_stellar": 0.85,
    "lc": {"t": [0.5, 1.0, 1.5], "f": [1.0, 0.99, 1.0]}
}

tmp = os.path.join(HERE, "tmp")
os.makedirs(tmp, exist_ok=True)
cand = os.path.join(tmp, "cand.txt")
with open(cand, "w") as fh:
    for k, v in lines.items():
        if k == "lc": continue
        fh.write(f"{k}={v}\n")

lc_path = os.path.join(tmp, "lc.csv")
with open(lc_path, "w") as fh:
    fh.write("time,flux,flux_err,model\n")
    for t, f in zip(lines["lc"]["t"], lines["lc"]["f"]):
        fh.write(f"{t},{f},0.001,1.0\n")

print("cand:", cand, os.path.exists(cand))
print("lc_path:", lc_path, os.path.exists(lc_path))

p = subprocess.run([BIN, cand, lc_path], capture_output=True, text=True, timeout=120)
print("returncode:", p.returncode)
if p.returncode != 0: print("stderr:", p.stderr[:500])
else: print("stdout:", p.stdout[:500])