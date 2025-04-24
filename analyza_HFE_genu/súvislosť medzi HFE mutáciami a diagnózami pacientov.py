import os
import pandas as pd
from scipy.stats import chi2_contingency

# zistenie absolutnej cesty k suboru (vzdy relativne k priecinku skriptu)
base_dir = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(base_dir, "..", "SSBU25_dataset_cleaned.csv")

# nacitanie datasetu
df = pd.read_csv(data_path, sep=";", encoding="utf-8-sig")

# vycistenie nazvov stlpcov (odstranenie \n a bielych znakov)
df.columns = df.columns.str.replace('\n', ' ', regex=False).str.strip()

# pecenove diagnozy, ktore skumame
pecen_diag = ["K76.0", "K75.9"]

# novy stlpec: True ak ma pacient pecenovu diagnozu
df["pecen"] = df["diagnoza MKCH-10"].isin(pecen_diag)

# mapovanie pre jednoduchsie pouzitie mutacii
mutacie = {
    "H63D": "HFE C187G (H63D) [HFE]",
    "S65C": "HFE A193T (S65C) [HFE]",
    "C282Y": "HFE G845A (C282Y) [HFE]"
}

# funkcia na analyzu jednej mutacie
def analyzuj_suvislost(mut_nazov, stlpec):
    # mutovany je vsetko okrem "normal" alebo "wt/wt"
    df[f"{mut_nazov}_mutovany"] = df[stlpec].apply(lambda x: x.strip().lower() != "normal")

    # kontingencna tabulka
    tab = pd.crosstab(df[f"{mut_nazov}_mutovany"], df["pecen"])
    print(f"\n📊 Suvislost medzi mutaciou {mut_nazov} a pecenovymi diagnozami:")
    print(tab)

    if tab.shape == (2, 2):
        chi2, p, _, _ = chi2_contingency(tab)
        print(f"Chi² = {chi2:.4f}, p-hodnota = {p:.4f}")
        if p < 0.05:
            print("❗ Statisticky vyznamna suvislost (p < 0.05)")
        else:
            print("✅ Ziadna statisticky vyznamna suvislost (p ≥ 0.05)")

# spustenie analyzy pre vsetky mutacie
for mut_nazov, stlpec in mutacie.items():
    analyzuj_suvislost(mut_nazov, stlpec)
