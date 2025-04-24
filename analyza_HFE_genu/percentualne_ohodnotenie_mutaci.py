import os
import pandas as pd

# zistenie absolutnej cesty k suboru
base_dir = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(base_dir, "..", "SSBU25_dataset_cleaned.csv")

# nacitanie datasetu
df = pd.read_csv(data_path, sep=";", encoding="utf-8-sig")

# ocistenie nazvov stlpcov od \n a medzier
df.columns = df.columns.str.replace('\n', ' ', regex=False).str.strip()

# mapovanie kratkych nazvov na skutocne stlpce
mutacie = {
    "C282Y": "HFE G845A (C282Y) [HFE]",
    "H63D": "HFE C187G (H63D) [HFE]",
    "S65C": "HFE A193T (S65C) [HFE]"
}

# funkcia na vypocet percent pre jeden gen
def genotyp_percenta(mut_key):
    stlpec = mutacie[mut_key]
    print(f"\n📊 Genotypove zastupenie pre {mut_key}:")
    counts = df[stlpec].value_counts()
    total = counts.sum()
    for gtype, count in counts.items():
        percent = (count / total) * 100
        print(f"{gtype}: {count} pacientov ({percent:.2f}%)")

# zisti percenta pre kazdu mutaciu
for mut in mutacie:
    genotyp_percenta(mut)

# vypocet prenasacov (heterozygoti)
pren = df[
    (df[mutacie["C282Y"]] == "wt/mut") |
    (df[mutacie["H63D"]] == "wt/mut") |
    (df[mutacie["S65C"]] == "wt/mut")
]

# geneticka predispozicia na HH
predispozicia = df[
    (df[mutacie["C282Y"]] == "mut/mut") |  # homozygot
    (
        (df[mutacie["C282Y"]] == "wt/mut") &
        (df[mutacie["H63D"]] == "wt/mut")  # zlozeny heterozygot
    )
]

# vypis
print(f"\n🧬 Prenasaci: {len(pren)} pacientov ({(len(pren)/len(df))*100:.2f}%)")
print(f"🧬 Geneticka predispozicia na HH: {len(predispozicia)} pacientov ({(len(predispozicia)/len(df))*100:.2f}%)")
