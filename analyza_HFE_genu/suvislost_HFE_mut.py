import os
import pandas as pd
from scipy.stats import chi2_contingency
import matplotlib.pyplot as plt
import seaborn as sns

def generate_suvislosti_diag(df):
    output_dir_tab = "tabulky"
    output_dir_graf = "grafy"
    os.makedirs(output_dir_tab, exist_ok=True)
    os.makedirs(output_dir_graf, exist_ok=True)

    df.columns = df.columns.str.replace('\n', ' ', regex=False).str.strip()

    pecen_diag = ["K76.0", "K75.9"]
    df["pecen"] = df["diagnoza MKCH-10"].isin(pecen_diag)

    mutacie = {
        "H63D": "HFE C187G (H63D) [HFE]",
        "S65C": "HFE A193T (S65C) [HFE]",
        "C282Y": "HFE G845A (C282Y) [HFE]"
    }

    vysledky = []

    for mut_nazov, stlpec in mutacie.items():
        if stlpec not in df.columns:
            continue

        df[f"{mut_nazov}_mutovany"] = df[stlpec].apply(lambda x: str(x).strip().lower() != "normal")

        tab = pd.crosstab(df[f"{mut_nazov}_mutovany"], df["pecen"])

        # Ulož tabulku
        tab.to_csv(os.path.join(output_dir_tab, f"suvislost_{mut_nazov}.csv"), sep=";", encoding="utf-8-sig")

        # Ak 2x2, spočítaj chi2
        if tab.shape == (2, 2):
            chi2, p, _, _ = chi2_contingency(tab)
            vysledky.append({
                "Mutácia": mut_nazov,
                "Chi²": round(chi2, 4),
                "p-hodnota": round(p, 4),
                "Významnosť": "Významné" if p < 0.05 else "Nevýznamné"
            })

        # Vytvor graf
        sns.set(style="whitegrid")
        plt.figure(figsize=(6, 4))
        tab.plot(kind="bar", stacked=True)
        plt.title(f"Súvislosť {mut_nazov} vs. pečeňová diagnóza")
        plt.xlabel("Mutovaný genotyp")
        plt.ylabel("Počet pacientov")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir_graf, f"graf_suvislost_{mut_nazov}.png"))
        plt.close()

    # Ulož prehľadné výsledky všetkých mutácií
    if vysledky:
        pd.DataFrame(vysledky).to_csv(os.path.join(output_dir_tab, "vysledky_suvislosti.csv"), index=False, sep=";", encoding="utf-8-sig")
