import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

def generate_mkch_analysis(df):
    output_dir = "grafy"
    os.makedirs(output_dir, exist_ok=True)
    tabulky_dir = "tabulky"
    os.makedirs(tabulky_dir, exist_ok=True)

    # Vyčistenie názvov stĺpcov
    df.columns = df.columns.str.replace('\n', ' ', regex=False).str.strip()

    # Extrahovanie prvej skupiny z MKCH-10
    df["mkch_skupina"] = df["diagnoza MKCH-10"].str[0].str.upper()

    # Mapovanie skupín
    skupiny = {
        "K": "Gastro a pečeň",
        "E": "Metabolické",
        "D": "Hematológia",
        "B": "Infekcie",
        "C": "Novotvary",
        "Z": "Iné zdravotné faktory"
    }
    df["diag_skupina"] = df["mkch_skupina"].map(skupiny).fillna("Ostatné/Neznáme")

    # Extrakcia roka z prijmu vzorky
    if "prijem_vzorky" in df.columns:
        df["rok"] = pd.to_datetime(df["prijem_vzorky"], errors="coerce", dayfirst=True).dt.year
    else:
        df["rok"] = pd.NaT

    # Agregácia
    diag_vyvoj = df.groupby(["rok", "diag_skupina"]).size().reset_index(name="pocet")
    diag_vyvoj = diag_vyvoj.dropna(subset=["rok"])

    # Uloženie agregácie do CSV
    diag_vyvoj.to_csv(os.path.join(tabulky_dir, "vyvoj_skupin_diagnoz.csv"), index=False, sep=";", encoding="utf-8-sig")

    # Vizualizácia vývoja v čase
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=diag_vyvoj, x="rok", y="pocet", hue="diag_skupina", marker="o")
    plt.title("Vývoj výskytu skupín diagnóz podľa rokov")
    plt.xlabel("Rok vyšetrenia")
    plt.ylabel("Počet pacientov")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "vyvoj_skupin_diagnoz.png"))
    plt.close()

    # Detekcia zastaraných kódov
    zastarale_kody_info = {
        "K76.0": {
            "novy_kod": "E66.9 (Obezita, nešpecifikovaná)",  # Príklad
            "popis": "K76.0 označuje tukovú degeneráciu pečene. V modernej klasifikácii sa rozlišuje NAFLD."
        },
        "K75.9": {
            "novy_kod": "K75.8 (Iné zápalové ochorenia pečene)",
            "popis": "K75.9 je nešpecifikovaná zápalová choroba pečene. Dnes sa odporúča špecifikovať zápal."
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
        zastarale_df.to_csv(os.path.join(tabulky_dir, "zastarale_kody_diagnoz.csv"), index=False, sep=";", encoding="utf-8-sig")

