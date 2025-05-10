import os
import pandas as pd
from scipy.stats import chi2

# Pomocná mapa pre premenovanie genotypov
genotyp_label = {
    "wt/wt": "wt/wt (normal)",
    "wt/mut": "wt/mut (heterozygot)",
    "mut/mut": "mut/mut (mutant)"
}

def generate_hwe_tables(df):
    output_dir = "tabulky"
    os.makedirs(output_dir, exist_ok=True)

    df.columns = df.columns.str.replace('\n', ' ', regex=False).str.strip()

    mutacie = {
        "H63D": "HFE C187G (H63D) [HFE]",
        "S65C": "HFE A193T (S65C) [HFE]",
        "C282Y": "HFE G845A (C282Y) [HFE]"
    }

    genotyp_map = {
        "normal": "wt/wt",
        "heterozygot": "wt/mut",
        "homozygot": "mut/mut"
    }

    vysledky = []

    for mut_nazov, stlpec in mutacie.items():
        if stlpec not in df.columns:
            continue

        mapped = df[stlpec].str.strip().str.lower().map(genotyp_map)

        counts = mapped.value_counts()
        obs_wtwt = counts.get("wt/wt", 0)
        obs_wtmut = counts.get("wt/mut", 0)
        obs_mutmut = counts.get("mut/mut", 0)
        n = obs_wtwt + obs_wtmut + obs_mutmut

        if n == 0:
            continue

        p = (2 * obs_wtwt + obs_wtmut) / (2 * n)
        q = 1 - p

        exp_wtwt = p ** 2 * n
        exp_wtmut = 2 * p * q * n
        exp_mutmut = q ** 2 * n

        chi2_stat = ((obs_wtwt - exp_wtwt) ** 2 / exp_wtwt if exp_wtwt > 0 else 0) + \
                    ((obs_wtmut - exp_wtmut) ** 2 / exp_wtmut if exp_wtmut > 0 else 0) + \
                    ((obs_mutmut - exp_mutmut) ** 2 / exp_mutmut if exp_mutmut > 0 else 0)

        p_value = 1 - chi2.cdf(chi2_stat, df=1)

        vysledok_text = "Odchýlka" if p_value < 0.05 else "Súlad"

        vysledky.append({
            "Mutácia": mut_nazov,
            "Pozorované wt/wt (normal)": obs_wtwt,
            "Pozorované wt/mut (heterozygot)": obs_wtmut,
            "Pozorované mut/mut (mutant)": obs_mutmut,
            "Očakávané wt/wt (normal)": round(exp_wtwt, 2),
            "Očakávané wt/mut (heterozygot)": round(exp_wtmut, 2),
            "Očakávané mut/mut (mutant)": round(exp_mutmut, 2),
            "Chi²": round(chi2_stat, 4),
            "p-hodnota": round(p_value, 4),
            "Výsledok": vysledok_text
        })

    vysledky_df = pd.DataFrame(vysledky)
    vysledky_df.to_csv(os.path.join(output_dir, "hardy_weinberg_test.csv"), index=False, sep=";", encoding="utf-8-sig")
