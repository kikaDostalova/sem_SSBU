import pandas as pd

# nacitanie datasetu
df = pd.read_csv("SSBU25_dataset_cleaned.csv")

def genotyp_percenta(mutacia):
    print(f"\n📊 Genotypové zastúpenie pre {mutacia}:")
    counts = df[mutacia].value_counts()
    total = counts.sum()
    for gtype, count in counts.items():
        percent = (count / total) * 100
        print(f"{gtype}: {count} pacientov ({percent:.2f}%)")

# zisti percenta pre kazdu mutaciu
for mut in ["C282Y", "H63D", "S65C"]:
    genotyp_percenta(mut)

# vypocet prenasacov
pren = df[(df["C282Y"] == "wt/mut") |
          (df["H63D"] == "wt/mut") |
          (df["S65C"] == "wt/mut")]

# geneticka predispozicia
predispozicia = df[
    (df["C282Y"] == "mut/mut") |  # homozygot C282Y
    ((df["C282Y"] == "wt/mut") & (df["H63D"] == "wt/mut"))  # zlozeny heterozygot
]

print(f"\n🧬 Prenášači: {len(pren)} pacientov ({(len(pren)/len(df))*100:.2f}%)")
print(f"🧬 Genetická predispozícia na HH: {len(predispozicia)} pacientov ({(len(predispozicia)/len(df))*100:.2f}%)")
