import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

def generate_mkch10_analysis():
    # Vytvorenie priečinku na uloženie grafov, ak neexistuje
    os.makedirs("grafy", exist_ok=True)

    # Načítanie dát
    df = pd.read_csv("sem_SSBU/SSBU25_dataset_cleaned.csv", sep=";", encoding="utf-8-sig")
    df.columns = df.columns.str.replace('\n', ' ', regex=False).str.strip()

    # Extrakcia písmena MKCH-10
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

    # Extrakcia roku
    if "validovany_vysledok" in df.columns:
        df["rok"] = pd.to_datetime(df["validovany_vysledok"], errors="coerce", dayfirst=True).dt.year
    else:
        df["rok"] = pd.NaT

    # Agregácia podľa rokov a skupín
    diag_vyvoj = df.groupby(["rok", "diag_skupina"]).size().reset_index(name="pocet")
    diag_vyvoj = diag_vyvoj.dropna(subset=["rok"])

    # Vizualizácia - ulož na disk
    sns.set(style="whitegrid")
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.lineplot(data=diag_vyvoj, x="rok", y="pocet", hue="diag_skupina", marker="o", ax=ax)
    plt.title("Vývoj výskytu skupín diagnóz podľa rokov")
    plt.xlabel("Rok vyšetrenia")
    plt.ylabel("Počet pacientov")
    plt.tight_layout()
    fig.savefig("grafy/graf_diag_vyvoj_v_case.png")
    plt.close(fig)

    # Kontrola neštandardných/zastaraných kódov
    podozrive_mask = ~df["diagnoza MKCH-10"].fillna("").str.match(r"^[A-Z][0-9]{2}(\.[0-9])?$")
    podozrive = df[podozrive_mask]

    # Návrat dataframe podozrivých kódov pre prípadné ďalšie spracovanie
    return podozrive[["id", "diagnoza MKCH-10"]]
