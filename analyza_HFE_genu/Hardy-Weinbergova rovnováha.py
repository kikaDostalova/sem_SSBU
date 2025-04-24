import os
import pandas as pd
from scipy.stats import chi2

# nacitanie datasetu
base_dir = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(base_dir, "..", "SSBU25_dataset_cleaned.csv")
df = pd.read_csv(data_path, sep=";", encoding="utf-8-sig")

# vycistenie nazvov stlpcov
df.columns = df.columns.str.replace('\n', ' ', regex=False).str.strip()

# mapovanie mutacii
mutacie = {
    "H63D": "HFE C187G (H63D) [HFE]",
    "S65C": "HFE A193T (S65C) [HFE]",
    "C282Y": "HFE G845A (C282Y) [HFE]"
}

# mapovanie textovych hodnot na genotypy
genotyp_map = {
    "normal": "wt/wt",
    "heterozygot": "wt/mut",
    "homozygot": "mut/mut"
}

# funkcia na test HWE
def test_hwe(genotypes, label):
    # nahrad hodnoty podla mapy
    mapped = genotypes.str.strip().str.lower().map(genotyp_map)

    # pocet jednotlivych genotypov
    counts = mapped.value_counts()
    obs_wtwt = counts.get("wt/wt", 0)
    obs_wtmut = counts.get("wt/mut", 0)
    obs_mutmut = counts.get("mut/mut", 0)
    n = obs_wtwt + obs_wtmut + obs_mutmut

    if n == 0:
        print(f"\n⚠️ Nedostatok údajov pre {label}")
        return

    p = (2 * obs_wtwt + obs_wtmut) / (2 * n)
    q = 1 - p

    exp_wtwt = p ** 2 * n
    exp_wtmut = 2 * p * q * n
    exp_mutmut = q ** 2 * n

    chi2_stat = ((obs_wtwt - exp_wtwt) ** 2 / exp_wtwt if exp_wtwt > 0 else 0) + \
                ((obs_wtmut - exp_wtmut) ** 2 / exp_wtmut if exp_wtmut > 0 else 0) + \
                ((obs_mutmut - exp_mutmut) ** 2 / exp_mutmut if exp_mutmut > 0 else 0)

    p_value = 1 - chi2.cdf(chi2_stat, df=1)

    print(f"\n🧬 Mutácia: {label}")
    print(f"Pozorované genotypy (wt/wt, wt/mut, mut/mut): {obs_wtwt}, {obs_wtmut}, {obs_mutmut}")
    print(f"Očakávané genotypy: {exp_wtwt:.2f}, {exp_wtmut:.2f}, {exp_mutmut:.2f}")
    print(f"Chi² = {chi2_stat:.4f}, p-hodnota = {p_value:.4f}")

    if p_value < 0.05:
        print("❗ Odchýlka od Hardy-Weinbergovej rovnováhy (p < 0.05)")
    else:
        print("✅ Genotypy sú v Hardy-Weinbergovej rovnováhe")

# spusti test pre kazdu mutaciu
for mut_nazov, stlpec in mutacie.items():
    test_hwe(df[stlpec], mut_nazov)
