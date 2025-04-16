import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# nacitanie dat
df = pd.read_csv("SSBU25_dataset_cleaned.csv")

# extrahuj zaciatocne pismeno z MKCH-10 kodu
df["mkch_skupina"] = df["diagnoza"].str[0].str.upper()

# mapovanie skupín
skupiny = {
    "K": "Gastro a pečeň",
    "E": "Metabolické",
    "D": "Hematológia",
    "B": "Infekcie",
    "C": "Novotvary",
    "Z": "Iné zdravotné faktory"
}

df["diag_skupina"] = df["mkch_skupina"].map(skupiny).fillna("Ostatné/Neznáme")

# extrahuj rok
df["rok"] = pd.to_datetime(df["datetime_odber"], errors="coerce").dt.year

# počet podľa skupín a rokov
diag_vyvoj = df.groupby(["rok", "diag_skupina"]).size().reset_index(name="pocet")

# vykreslenie vývoja
plt.figure(figsize=(10, 6))
sns.lineplot(data=diag_vyvoj, x="rok", y="pocet", hue="diag_skupina", marker="o")
plt.title("Vývoj výskytu skupín diagnóz podľa rokov")
plt.xlabel("Rok vyšetrenia")
plt.ylabel("Počet pacientov")
plt.tight_layout()
plt.savefig("graf_diag_vyvoj_v_case.png")
plt.show()

# kontrola podozrivých/zastaraných kódov
podozrive = df[~df["diagnoza"].str.match(r"^[A-Z][0-9]{2}(\.[0-9])?$")]
print("\n❗ Možné neštandardné alebo zastarané kódy diagnóz:")
print(podozrive["diagnoza"].value_counts())
