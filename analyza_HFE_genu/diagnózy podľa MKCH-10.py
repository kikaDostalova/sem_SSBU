import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


# nacitanie dat
df = pd.read_csv("../SSBU25_dataset_cleaned.csv", sep=";", encoding="utf-8-sig")
df.columns = df.columns.str.replace('\n', ' ', regex=False).str.strip()

# extrahuj zaciatocne pismeno z MKCH-10 kodu
df["mkch_skupina"] = df["diagnoza MKCH-10"].str[0].str.upper()

# mapovanie skupin
skupiny = {
    "K": "Gastro a pečeň",
    "E": "Metabolické",
    "D": "Hematológia",
    "B": "Infekcie",
    "C": "Novotvary",
    "Z": "Iné zdravotné faktory"
}
df["diag_skupina"] = df["mkch_skupina"].map(skupiny).fillna("Ostatné/Neznáme")

# extrakcia roku z validovaneho vysledku
if "validovany_vysledok" in df.columns:
    df["rok"] = pd.to_datetime(df["validovany_vysledok"], errors="coerce", dayfirst=True).dt.year
else:
    df["rok"] = pd.NaT

# agregacia podla rokov a skupin
diag_vyvoj = df.groupby(["rok", "diag_skupina"]).size().reset_index(name="pocet")
diag_vyvoj = diag_vyvoj.dropna(subset=["rok"])

# vizualizacia
sns.set(style="whitegrid")
plt.figure(figsize=(10, 6))
sns.lineplot(data=diag_vyvoj, x="rok", y="pocet", hue="diag_skupina", marker="o")
plt.title("Vývoj výskytu skupín diagnóz podľa rokov")
plt.xlabel("Rok vyšetrenia")
plt.ylabel("Počet pacientov")
plt.tight_layout()
plt.savefig("graf_diag_vyvoj_v_case.png")
plt.show()

# kontrola neplatných kódov MKCH-10
podozrive_mask = ~df["diagnoza MKCH-10"].fillna("").str.match(r"^[A-Z][0-9]{2}(\.[0-9])?$")
podozrive = df[podozrive_mask]

print("\n❗ Možné neštandardné alebo zastarané kódy diagnóz:")
print(podozrive["diagnoza MKCH-10"].value_counts())
