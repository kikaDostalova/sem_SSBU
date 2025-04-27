import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


# nacitanie dat
df = pd.read_csv("../SSBU25_dataset_cleaned.csv", sep=";", encoding="utf-8-sig")
df.columns = df.columns.str.replace('\n', ' ', regex=False).str.strip()

# oprava preklepu v nazve stlpca
if "pohavie" in df.columns and "pohlavie" not in df.columns:
    df = df.rename(columns={"pohavie": "pohlavie"})

# mapovanie mutacii
mutacie = {
    "H63D": "HFE C187G (H63D) [HFE]",
    "S65C": "HFE A193T (S65C) [HFE]",
    "C282Y": "HFE G845A (C282Y) [HFE]"
}

# mapovanie genotypov
genotyp_map = {
    "normal": "wt/wt",
    "heterozygot": "wt/mut",
    "homozygot": "mut/mut"
}

# aplikuj mapovanie na genotypy
for mut_key, col_name in mutacie.items():
    df[mut_key] = df[col_name].str.strip().str.lower().map(genotyp_map)

# extrakcia hlavneho veku
df["vek"] = df["vek"].str.replace(',', '.', regex=False).str.extract(r"(\d+\.\d+)").astype(float)

# pečenové diagnózy
df["pecen_diag"] = df["diagnoza MKCH-10"].isin(["K76.0", "K75.9"])

# priecinok na grafy
output_dir = "grafy"
os.makedirs(output_dir, exist_ok=True)
sns.set(style="whitegrid")

### 1. Rozdelenie genotypov
for mut in ["H63D", "S65C", "C282Y"]:
    plt.figure(figsize=(6, 4))
    sns.countplot(x=df[mut], order=["wt/wt", "wt/mut", "mut/mut"])
    plt.title(f"Rozdelenie genotypov – {mut}")
    plt.xlabel("Genotyp")
    plt.ylabel("Počet pacientov")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"rozdelenie_{mut}.png"))
    plt.show()

### 2. Genotyp vs. vek
for mut in ["H63D", "S65C", "C282Y"]:
    plt.figure(figsize=(6, 4))
    sns.boxplot(x=df[mut], y=df["vek"], order=["wt/wt", "wt/mut", "mut/mut"])
    plt.title(f"Vek pacientov podľa genotypu – {mut}")
    plt.xlabel("Genotyp")
    plt.ylabel("Vek")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"vek_vs_genotyp_{mut}.png"))
    plt.show()

### 3. Genotyp vs. pohlavie
for mut in ["H63D", "S65C", "C282Y"]:
    plt.figure(figsize=(6, 4))
    sns.countplot(data=df, x=mut, hue="pohlavie", order=["wt/wt", "wt/mut", "mut/mut"])
    plt.title(f"Genotyp vs. Pohlavie – {mut}")
    plt.xlabel("Genotyp")
    plt.ylabel("Počet pacientov")
    plt.legend(title="Pohlavie")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"pohlavie_vs_genotyp_{mut}.png"))
    plt.show()

### 4. Genotyp vs. pečeňová diagnóza
for mut in ["H63D", "S65C", "C282Y"]:
    plt.figure(figsize=(6, 4))
    sns.countplot(data=df, x=mut, hue="pecen_diag", order=["wt/wt", "wt/mut", "mut/mut"])
    plt.title(f"Genotyp vs. Ochorenia pečene – {mut}")
    plt.xlabel("Genotyp")
    plt.ylabel("Počet pacientov")
    plt.legend(title="Pečeňová diagnóza", labels=["Nie", "Áno"])
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"pecen_vs_genotyp_{mut}.png"))
    plt.show()
