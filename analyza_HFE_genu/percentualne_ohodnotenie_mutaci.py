import os
import pandas as pd

def generate_percenta_genotypov(df):
    output_dir = "tabulky"
    os.makedirs(output_dir, exist_ok=True)

    df.columns = df.columns.str.replace('\n', ' ', regex=False).str.strip()

    mutacie = {
        "C282Y": "HFE G845A (C282Y) [HFE]",
        "H63D": "HFE C187G (H63D) [HFE]",
        "S65C": "HFE A193T (S65C) [HFE]"
    }

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

    vysledky = []

    for mut_key, stlpec in mutacie.items():
        if stlpec not in df.columns:
            continue

        genotypy = df[stlpec].str.strip().str.lower().map(genotyp_map)

        total = len(genotypy.dropna())
        counts = genotypy.value_counts()

        for genotyp, pocet in counts.items():
            vysledky.append({
                "Mutácia": mut_key,
                "Genotyp": genotyp_label.get(genotyp, genotyp),
                "Počet pacientov": pocet,
                "Percento": round((pocet / total) * 100, 2)
            })

    # Prenášači a predispozícia - celkové sumárne hodnoty
    pren = df[
        (df[mutacie["C282Y"]].str.strip().str.lower() == "heterozygot") |
        (df[mutacie["H63D"]].str.strip().str.lower() == "heterozygot") |
        (df[mutacie["S65C"]].str.strip().str.lower() == "heterozygot")
    ]
    predispozicia = df[
        (df[mutacie["C282Y"]].str.strip().str.lower() == "homozygot") |
        (
            (df[mutacie["C282Y"]].str.strip().str.lower() == "heterozygot") &
            (df[mutacie["H63D"]].str.strip().str.lower() == "heterozygot")
        )
    ]

    vysledky_df = pd.DataFrame(vysledky)

    # Exportuj tabuľku
    vysledky_df.to_csv(os.path.join(output_dir, "percenta_genotypov.csv"), index=False, sep=";", encoding="utf-8-sig")

    # Ulož sumár info o prenášačoch a predispozícii
    summary_path = os.path.join(output_dir, "prenasic_predispozicia.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(f"Prenášači: {len(pren)} pacientov ({(len(pren)/len(df))*100:.2f}%)\n")
        f.write(f"Genetická predispozícia na HH: {len(predispozicia)} pacientov ({(len(predispozicia)/len(df))*100:.2f}%)\n")
