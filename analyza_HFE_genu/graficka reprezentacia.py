import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# nacitaj ocisteny dataset
df = pd.read_csv("SSBU25_dataset_cleaned.csv")
sns.set(style="whitegrid")

# priecinok na ulozenie grafov
output_dir = "grafy"
os.makedirs(output_dir, exist_ok=True)

### 1. Rozdelenie genotypov – barploty
for mut in ["H63D", "S65C", "C282Y"]:
    plt.figure(figsize=(6, 4))
    sns.countplot(x=df[mut], order=["wt/wt", "wt/mut", "mut/mut"])
    plt.title(f"Rozdelenie genotypov – {mut}")
    plt.xlabel("Genotyp")
    plt.ylabel("Počet pacientov")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"rozdelenie_{mut}.png"))
    plt.show()

### 2. Genotyp vs. vek – boxploty
for mut in ["H63D", "S65C", "C282Y"]:
    plt.figure(figsize=(6, 4))
    sns.boxplot(x=df[mut], y=df["vek"], order=["wt/wt", "wt/mut", "mut/mut"])
    plt.title(f"Vek pacientov podľa genotypu – {mut}")
    plt.xlabel("Genotyp")
    plt.ylabel("Vek")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"vek_vs_genotyp_{mut}.png"))
    plt.show()

### 3. Genotyp vs. pohlavie – zhlukový stĺpcový graf
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
df["pecen_diag"] = df["diagnoza"].isin(["K76.0", "K75.9"])

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
