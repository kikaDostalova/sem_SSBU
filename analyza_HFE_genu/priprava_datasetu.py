import pandas as pd
import numpy as np

# nacitanie CSV suboru so spravnym oddelovacom a ignorovanim zbytocnych prazdnych stlpcov
df = pd.read_csv("SSBU25_dataset.csv", sep=";", usecols=range(12))

# oprava hlaviciek stlpcov (jednoduchsie nazvy)
df.columns = [
    "id", "validovany_vysledok", "datum_odber", "cas_odber",
    "datum_prijem", "cas_prijem", "pohlavie", "vek",
    "diagnoza", "H63D", "S65C", "C282Y"
]

# spojenie datumov a casov do datetime
df["datetime_odber"] = pd.to_datetime(df["datum_odber"] + " " + df["cas_odber"], errors="coerce", dayfirst=True)
df["datetime_prijem"] = pd.to_datetime(df["datum_prijem"] + " " + df["cas_prijem"], errors="coerce", dayfirst=True)

# odstranenie povodnych stlpcov datumu/casu
df.drop(columns=["datum_odber", "cas_odber", "datum_prijem", "cas_prijem"], inplace=True)

# uprava veku – vezme prvu hodnotu pred ciarkou a prevedie na int
df["vek"] = df["vek"].astype(str).str.split(",").str[0].str.strip()
df["vek"] = pd.to_numeric(df["vek"], errors='coerce')

# standardizacia genotypov
genotyp_map = {
    "normal": "wt/wt",
    "heterozygot": "wt/mut",
    "mutant": "mut/mut",
    "homozygot": "mut/mut"  # pre istotu aj ak by sa vyskytol
}

for col in ["H63D", "S65C", "C282Y"]:
    df[col] = df[col].map(genotyp_map).fillna(df[col])

# vypis par zaznamov
print(df.head())

# ulozenie cisteho datasetu
df.to_csv("SSBU25_dataset_cleaned.csv", index=False)
