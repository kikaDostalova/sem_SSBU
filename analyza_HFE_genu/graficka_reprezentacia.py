import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

def generate_graphs(df):
    # Priečinok na uloženie grafov
    output_dir = "grafy"
    os.makedirs(output_dir, exist_ok=True)
    
    sns.set(style="whitegrid")

    # Oprava preklepu v názve stĺpca
    if "pohavie" in df.columns and "pohlavie" not in df.columns:
        df = df.rename(columns={"pohavie": "pohlavie"})

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

    # Aplikuj mapovanie na genotypy
    for mut_key, col_name in mutacie.items():
        if col_name in df.columns:
            df[mut_key] = df[col_name].str.strip().str.lower().map(genotyp_map)

    # Extrakcia veku
    df["vek"] = df["vek"].astype(str).str.replace(',', '.', regex=False).str.extract(r"(\d+\.\d+)").astype(float)

    # Pečeňové diagnózy
    df["pecen_diag"] = df["diagnoza MKCH-10"].isin(["K76.0", "K75.9"])

    ### 1. Rozdelenie genotypov
    for mut in ["H63D", "S65C", "C282Y"]:
        fig, ax = plt.subplots(figsize=(6, 4))
        sns.countplot(x=df[mut], order=["wt/wt", "wt/mut", "mut/mut"], ax=ax)
        ax.set_title(f"Rozdelenie genotypov – {mut}")
        ax.set_xlabel("Genotyp")
        ax.set_ylabel("Počet pacientov")
        plt.tight_layout()
        fig.savefig(os.path.join(output_dir, f"rozdelenie_{mut}.png"))
        plt.close(fig)

    ### 2. Genotyp vs. vek
    for mut in ["H63D", "S65C", "C282Y"]:
        fig, ax = plt.subplots(figsize=(6, 4))
        sns.boxplot(x=df[mut], y=df["vek"], order=["wt/wt", "wt/mut", "mut/mut"], ax=ax)
        ax.set_title(f"Vek pacientov podľa genotypu – {mut}")
        ax.set_xlabel("Genotyp")
        ax.set_ylabel("Vek")
        plt.tight_layout()
        fig.savefig(os.path.join(output_dir, f"vek_vs_genotyp_{mut}.png"))
        plt.close(fig)

    ### 3. Genotyp vs. pohlavie
    for mut in ["H63D", "S65C", "C282Y"]:
        fig, ax = plt.subplots(figsize=(6, 4))
        sns.countplot(data=df, x=mut, hue="pohlavie", order=["wt/wt", "wt/mut", "mut/mut"], ax=ax)
        ax.set_title(f"Genotyp vs. Pohlavie – {mut}")
        ax.set_xlabel("Genotyp")
        ax.set_ylabel("Počet pacientov")
        ax.legend(title="Pohlavie")
        plt.tight_layout()
        fig.savefig(os.path.join(output_dir, f"pohlavie_vs_genotyp_{mut}.png"))
        plt.close(fig)

    ### 4. Genotyp vs. pečeňová diagnóza
    for mut in ["H63D", "S65C", "C282Y"]:
        fig, ax = plt.subplots(figsize=(6, 4))
        sns.countplot(data=df, x=mut, hue="pecen_diag", order=["wt/wt", "wt/mut", "mut/mut"], ax=ax)
        ax.set_title(f"Genotyp vs. Ochorenia pečene – {mut}")
        ax.set_xlabel("Genotyp")
        ax.set_ylabel("Počet pacientov")
        ax.legend(title="Pečeňová diagnóza", labels=["Nie", "Áno"])
        plt.tight_layout()
        fig.savefig(os.path.join(output_dir, f"pecen_vs_genotyp_{mut}.png"))
        plt.close(fig)
