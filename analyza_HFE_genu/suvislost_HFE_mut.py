import os
import pandas as pd
from scipy.stats import chi2_contingency

def generate_suvislosti_diag(df):
    output_dir = "tabulky"
    os.makedirs(output_dir, exist_ok=True)

    # Vyčistenie názvov stĺpcov
    df.columns = df.columns.str.replace('\n', ' ', regex=False).str.strip()

    # Príprava pečeňových diagnóz
    df["pecen_diag"] = df["diagnoza MKCH-10"].isin(["K76.0", "K75.9"])

    # Mapovanie mutácií
    mutacie = {
        "H63D": "HFE C187G (H63D) [HFE]",
        "S65C": "HFE A193T (S65C) [HFE]",
        "C282Y": "HFE G845A (C282Y) [HFE]"
    }

    # Mapovanie genotypov
    genotyp_map = {
        "normal": "wt/wt",
        "heterozygot": "wt/mut",
        "homozygot": "mut/mut"
    }

    genotyp_label = {
        "wt/wt": "wt/wt (normal)",
        "wt/mut": "wt/mut (heterozygot)",
        "mut/mut": "mut/mut (mutant)"
    }

    summary = []

    for mut_key, col_name in mutacie.items():
        if col_name not in df.columns:
            continue

        # Prevod na štandardné genotypy
        genotypy = df[col_name].str.strip().str.lower().map(genotyp_map)

        kont_tab = pd.crosstab(genotypy, df["pecen_diag"])
        kont_tab.index = kont_tab.index.map(lambda x: genotyp_label.get(x, x))
        kont_path = os.path.join(output_dir, f"suvislost_{mut_key}.csv")
        kont_tab.to_csv(kont_path, sep=";", encoding="utf-8-sig")

        pval = None
        chi2_stat = None

        if kont_tab.shape in [(3, 2), (2, 2)]:
            try:
                chi2_stat, pval, dof, expected = chi2_contingency(kont_tab)
                vyznamnost = "významné" if pval < 0.05 else "nevýznamné"
                summary.append({
                    "Mutácia": mut_key,
                    "Chi²": round(chi2_stat, 4),
                    "p-hodnota": round(pval, 4),
                    "Významnosť": vyznamnost
                })
            except ValueError as e:
                summary.append({
                    "Mutácia": mut_key,
                    "Chi²": "–",
                    "p-hodnota": "–",
                    "Významnosť": f"Chyba: {str(e)}"
                })
        else:
            summary.append({
                "Mutácia": mut_key,
                "Chi²": "–",
                "p-hodnota": "–",
                "Významnosť": "Neštandardný tvar tabuľky"
            })

    # Uloženie sumárnej tabuľky
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(os.path.join(output_dir, "vysledky_suvislosti.csv"), sep=";", index=False, encoding="utf-8-sig")
