import pandas as pd
from scipy.stats import chi2

# nacitaj ocisteny dataset
df = pd.read_csv("SSBU25_dataset_cleaned.csv")


# definuj funkciu pre test HWE
def test_hwe(genotypes, label):
    # spocitaj genotypy
    counts = genotypes.value_counts()
    obs_wtwt = counts.get("wt/wt", 0)
    obs_wtmut = counts.get("wt/mut", 0)
    obs_mutmut = counts.get("mut/mut", 0)
    n = obs_wtwt + obs_wtmut + obs_mutmut

    # alelove frekvencie
    p = (2 * obs_wtwt + obs_wtmut) / (2 * n)
    q = 1 - p

    # ocakavane pocty
    exp_wtwt = p ** 2 * n
    exp_wtmut = 2 * p * q * n
    exp_mutmut = q ** 2 * n

    # chi^2 test
    chi2_stat = ((obs_wtwt - exp_wtwt) ** 2 / exp_wtwt) + \
                ((obs_wtmut - exp_wtmut) ** 2 / exp_wtmut) + \
                ((obs_mutmut - exp_mutmut) ** 2 / exp_mutmut)

    p_value = 1 - chi2.cdf(chi2_stat, df=1)

    # vysledok
    print(f"\n🧬 Mutácia: {label}")
    print(f"Pozorované genotypy (wt/wt, wt/mut, mut/mut): {obs_wtwt}, {obs_wtmut}, {obs_mutmut}")
    print(f"Očakávané genotypy: {exp_wtwt:.2f}, {exp_wtmut:.2f}, {exp_mutmut:.2f}")
    print(f"Chi² = {chi2_stat:.4f}, p-hodnota = {p_value:.4f}")

    if p_value < 0.05:
        print("❗ Odchýlka od Hardy-Weinbergovej rovnováhy (p < 0.05)")
    else:
        print("✅ Genotypy sú v Hardy-Weinbergovej rovnováhe")


# spusti test pre kazdu mutaciu
for mut in ["H63D", "S65C", "C282Y"]:
    test_hwe(df[mut], mut)
