import pandas as pd
from scipy.stats import chi2_contingency

# nacitanie datasetu
df = pd.read_csv("SSBU25_dataset_cleaned.csv")

# peceňové diagnózy, ktoré skúmame
pecen_diag = ["K76.0", "K75.9"]

# pomocná funkcia – zistí, či má pacient pečeňovú diagnózu
df["pecen"] = df["diagnoza"].isin(pecen_diag)

# funkcia na analýzu pre jednu mutáciu
def analyzuj_suvislost(mutacia):
    # mutovaný = všetko okrem wt/wt
    df[f"{mutacia}_mutovany"] = df[mutacia].apply(lambda x: x != "wt/wt")

    # vytvor kontingenčnú tabuľku
    tab = pd.crosstab(df[f"{mutacia}_mutovany"], df["pecen"])
    print(f"\n📊 Súvislosť medzi mutáciou {mutacia} a pečeňovými diagnózami:")
    print(tab)

    if tab.shape == (2, 2):  # kontrola tvaru
        chi2, p, _, _ = chi2_contingency(tab)
        print(f"Chi² = {chi2:.4f}, p-hodnota = {p:.4f}")
        if p < 0.05:
            print("❗ Štatisticky významná súvislosť (p < 0.05)")
        else:
            print("✅ Žiadna štatisticky významná súvislosť (p ≥ 0.05)")

# analyzuj pre všetky 3 mutácie
for mut in ["H63D", "S65C", "C282Y"]:
    analyzuj_suvislost(mut)
