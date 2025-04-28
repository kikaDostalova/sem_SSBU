import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

def generate_mkch10_analysis(df):
    output_dir_tab = "tabulky"
    output_dir_graf = "grafy"
    os.makedirs(output_dir_tab, exist_ok=True)
    os.makedirs(output_dir_graf, exist_ok=True)

    df.columns = df.columns.str.replace('\n', ' ', regex=False).str.strip()

    df["mkch_skupina"] = df["diagnoza MKCH-10"].str[0].str.upper()

    skupiny = {
        "K": "Gastro a pečeň",
        "E": "Metabolické",
        "D": "Hematológia",
        "B": "Infekcie",
        "C": "Novotvary",
        "Z": "Iné zdravotné faktory"
    }
    df["diag_skupina"] = df["mkch_skupina"].map(skupiny).fillna("Ostatné/Neznáme")

    if "validovany_vysledok" in df.columns:
        df["rok"] = pd.to_datetime(df["validovany_vysledok"], errors="coerce", dayfirst=True).dt.year
    else:
        df["rok"] = pd.NaT

    diag_vyvoj = df.groupby(["rok", "diag_skupina"]).size().reset_index(name="pocet")
    diag_vyvoj = diag_vyvoj.dropna(subset=["rok"])

    # Uloženie agregácie
    diag_vyvoj.to_csv(os.path.join(output_dir_tab, "vyvoj_skupin_diagnoz.csv"), index=False, sep=";", encoding="utf-8-sig")

    # Vizualizácia
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=diag_vyvoj, x="rok", y="pocet", hue="diag_skupina", marker="o")
    plt.title("Vývoj výskytu skupín diagnóz podľa rokov")
    plt.xlabel("Rok vyšetrenia")
    plt.ylabel("Počet pacientov")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir_graf, "vyvoj_skupin_diagnoz.png"))
    plt.close()

    # Detekcia zastaralých kódov
    zastarale_kody_info = {
        "K76.0": {
            "novy_kod": "E66.9 (Obezita, nešpecifikovaná)",
            "popis": "Tuková degenerácia pečene."
        },
        "K75.9": {
            "novy_kod": "K75.8 (Iné zápalové ochorenia pečene)",
            "popis": "Nešpecifikovaná zápalová choroba pečene."
        }
    }

    zastarale_records = []
    for kod, info in zastarale_kody_info.items():
        pocet = (df["diagnoza MKCH-10"] == kod).sum()
        if pocet > 0:
            zastarale_records.append({
                "Pôvodný kód": kod,
                "Počet pacientov": pocet,
                "Nový odporúčaný kód": info["novy_kod"],
                "Popis": info["popis"]
            })

    if zastarale_records:
        zastarale_df = pd.DataFrame(zastarale_records)
        zastarale_df.to_csv(os.path.join(output_dir_tab, "zastarale_kody_diagnoz.csv"), index=False, sep=";", encoding="utf-8-sig")
