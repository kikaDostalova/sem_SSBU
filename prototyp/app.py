import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy.stats import chi2_contingency, chi2
import export_report
st.set_page_config(page_title="Analýza HFE génu", layout="wide")
sns.set(style="whitegrid")

st.title("Analýza HFE génu a genetická predispozícia na hemochromatózu")

# Mutácie a mapovanie stĺpcov
mutacie_vzory = {
    "C282Y": "G845A",
    "H63D": "C187G",
    "S65C": "A193T"
}

# Labely pre genotypy
genotyp_label = {
    "wt/wt": "wt/wt (normal)",
    "wt/mut": "wt/mut (heterozygot)",
    "mut/mut": "mut/mut (mutant)"
}

# Pomocná funkcia na premapovanie výsledkov na genotypy
def premapuj_na_genotyp(series):
    return series.map({
        "normal": "wt/wt",
        "heterozygot": "wt/mut",
        "mutant": "mut/mut",
        "mutácia": "mut/mut",
        "patogénna": "mut/mut"
    }).fillna("neznamy")

# 📥 Upload a čistenie datasetu v Sidebar
st.sidebar.header("📥 Nahrajte CSV súbor")
uploaded_file = st.sidebar.file_uploader("Vyberte súbor:", type=["csv"])

if uploaded_file:
    with st.spinner("🔄 Čistenie datasetu..."):
        # Načítanie súboru
        df = pd.read_csv(uploaded_file, sep=";", encoding="utf-8-sig", dtype=str)

        # Vytvorenie nových dátumových stĺpcov
        df["validovany_vysledok"] = df.iloc[:, 1].fillna('') + " " + df.iloc[:, 2].fillna('')
        df["prijem_vzorky"] = df.iloc[:, 3].fillna('') + " " + df.iloc[:, 4].fillna('')

        # Odstránenie nepotrebných stĺpcov
        drop_cols = ['Unnamed: 3', 'Unnamed: 5'] + [col for col in df.columns if "Unnamed" in col]
        df.drop(columns=drop_cols + [df.columns[1], df.columns[2], df.columns[3], df.columns[4]], inplace=True, errors="ignore")

        # Odstránenie prázdnych ID
        df = df[df["id"].notna() & (df["id"].str.strip() != "")]

        # Validácia veku
        def is_valid_age(age_str):
            try:
                age_main = float(age_str.split(",")[0].replace(",", "."))
                return age_main <= 120
            except:
                return False

        df = df[df["vek"].apply(is_valid_age)]
        df["vek"] = df["vek"].astype(str).str.replace(",", ".").astype(float)

        # Normalize názvy stĺpcov
        df.columns = df.columns.str.replace('\\n', ' ', regex=False).str.strip()

        # Dynamické mapovanie názvov stĺpcov pre mutácie
        mutacie = {}
        for mut_short, expected_part in mutacie_vzory.items():
            matching_cols = [col for col in df.columns if expected_part in col and "HFE" in col]
            if matching_cols:
                mutacie[mut_short] = matching_cols[0]

        st.success("✅ Dataset bol úspešne načítaný a očistený.")
else:
    st.warning("⬆️ Nahrajte súbor vľavo v sidebare, aby sa zobrazili ďalšie sekcie.")
    st.stop()

#  Vyber sekciu v Sidebar
vyber_sekciu = st.sidebar.radio(
    "Vyberte sekciu",
    ("Úvodná analýza", "Grafy", "Analýza diagnóz podľa MKCH-10")
)

# =================== Sekcia: Úvodná analýza ===================
if vyber_sekciu == "Úvodná analýza":

    st.header("🧬 Hardy-Weinbergova rovnováha")

    def hardy_weinberg_test(genotypy):
        obs = genotypy.value_counts().reindex(["wt/wt", "wt/mut", "mut/mut"], fill_value=0)
        n = obs.sum()
        if n == 0:
            return None
        p = (2*obs["wt/wt"] + obs["wt/mut"]) / (2*n)
        q = 1 - p
        exp = [p**2 * n, 2*p*q * n, q**2 * n]
        chi2_stat = sum((obs - exp)**2 / exp)
        pval = 1 - chi2.cdf(chi2_stat, df=1)
        df_result = pd.DataFrame({
            "Pozorované": obs,
            "Očakávané": exp,
            "p-hodnota": [pval]*3
        })
        df_result.index = df_result.index.map(lambda x: genotyp_label.get(x, x))
        return df_result

    for mut, col in mutacie.items():
        st.subheader(f"Mutácia {mut}")
        if col in df.columns:
            genotypy = premapuj_na_genotyp(df[col])
            result = hardy_weinberg_test(genotypy)
            if result is not None:
                st.dataframe(result)
            else:
                st.warning(f"Mutácia {mut}: Nedostatočné dáta pre test.")
        else:
            st.warning(f"Stĺpec pre {mut} ({col}) sa nenašiel v datasete.")

    st.header("📊 Percentá genotypov a prenášači")

    selected = st.selectbox("Vyber mutáciu", list(mutacie.keys()))
    col = mutacie[selected]

    if col in df.columns:
        genotypy = premapuj_na_genotyp(df[col])
        counts = genotypy.value_counts(normalize=True).round(3) * 100
        counts.index = counts.index.map(lambda x: genotyp_label.get(x, x))
        total = len(genotypy)
        pren = genotypy.isin(["wt/mut", "mut/mut"]).sum()
        pred = genotypy.eq("mut/mut").sum()
        st.markdown(f"- **Prenášači**: {pren} pacientov ({(pren/total)*100:.1f} %)")
        st.markdown(f"- **Predispozícia**: {pred} pacientov ({(pred/total)*100:.1f} %)")
        st.dataframe(counts.rename("Percentá").to_frame())
    else:
        st.warning(f"Mutácia {selected} nie je dostupná v datasete.")

    st.header("📚 Súvislosť HFE mutácií s pečeňovými diagnózami")

    df["pecen_diag"] = df["diagnoza MKCH-10"].isin(["K76.0", "K75.9"])

    for mut, col in mutacie.items():
        if col in df.columns:
            st.subheader(f"Mutácia {mut}")
            genotypy = premapuj_na_genotyp(df[col])
            kont_tab = pd.crosstab(genotypy, df["pecen_diag"])
            kont_tab.index = kont_tab.index.map(lambda x: genotyp_label.get(x, x))
            st.dataframe(kont_tab)

            if kont_tab.shape[1] == 2:
                chi2_stat, pval, dof, expected = chi2_contingency(kont_tab)
                st.markdown(f"**p-hodnota**: `{pval:.4f}` {'(významné)' if pval < 0.05 else '(nevýznamné)'}")
            else:
                st.warning("Tabuľka nemá správny tvar na výpočet chi2 testu.")
        else:
            st.warning(f"Stĺpec pre {mut} ({col}) sa nenašiel v datasete.")

# =================== Sekcia: Grafy ===================
elif vyber_sekciu == "Grafy":
    st.header("📊 Grafy - rozdelenie genotypov, veku, pohlavia a diagnóz")

    st.subheader("Vyber mutáciu")
    selected_mut = st.selectbox("Vyber mutáciu", list(mutacie.keys()), key="grafy")

    col = mutacie[selected_mut]

    if col in df.columns:
        genotypy = premapuj_na_genotyp(df[col])

        # 1. Rozdelenie genotypov
        st.subheader("📊 Rozdelenie genotypov")
        fig1, ax1 = plt.subplots()
        sns.countplot(x=genotypy.map(lambda x: genotyp_label.get(x, x)), ax=ax1, order=[genotyp_label["wt/wt"], genotyp_label["wt/mut"], genotyp_label["mut/mut"]])
        ax1.set_xlabel("Genotyp")
        ax1.set_ylabel("Počet pacientov")
        st.pyplot(fig1)

        # 2. Vzťah genotypu a veku
        st.subheader("📈 Vek podľa genotypu")
        fig2, ax2 = plt.subplots()
        sns.boxplot(x=genotypy.map(lambda x: genotyp_label.get(x, x)), y=df["vek"], ax=ax2, order=[genotyp_label["wt/wt"], genotyp_label["wt/mut"], genotyp_label["mut/mut"]])
        ax2.set_xlabel("Genotyp")
        ax2.set_ylabel("Vek")
        st.pyplot(fig2)

        # 3. Vzťah genotypu a pohlavia
        st.subheader("🚻 Pohlavie podľa genotypu")
        fig3, ax3 = plt.subplots()
        sns.countplot(x=genotypy.map(lambda x: genotyp_label.get(x, x)), hue=df["pohavie"], ax=ax3, order=[genotyp_label["wt/wt"], genotyp_label["wt/mut"], genotyp_label["mut/mut"]])
        ax3.set_xlabel("Genotyp")
        ax3.set_ylabel("Počet pacientov")
        st.pyplot(fig3)

        # 4. Vzťah genotypu a diagnózy
        st.subheader("🩺 Diagnózy podľa genotypu (pečeňové ochorenia)")
        df["pecen_diag"] = df["diagnoza MKCH-10"].isin(["K76.0", "K75.9"])
        fig4, ax4 = plt.subplots()
        sns.countplot(x=genotypy.map(lambda x: genotyp_label.get(x, x)), hue=df["pecen_diag"], ax=ax4, order=[genotyp_label["wt/wt"], genotyp_label["wt/mut"], genotyp_label["mut/mut"]])
        ax4.set_xlabel("Genotyp")
        ax4.set_ylabel("Počet pacientov")
        ax4.legend(title="Pečeňové ochorenie", labels=["Nie", "Áno"])
        st.pyplot(fig4)
    else:
        st.warning(f"Mutácia {selected_mut} nie je dostupná v datasete.")

# =================== Sekcia: Analýza MKCH-10 ===================
elif vyber_sekciu == "Analýza diagnóz podľa MKCH-10":
    st.header("📋 Analýza diagnóz podľa MKCH-10 a ich vývoj v čase")

    # Skontroluj či existuje stĺpec s diagnózami
    if "diagnoza MKCH-10" not in df.columns:
        st.error("Dataset neobsahuje stĺpec 'diagnoza MKCH-10'.")
        st.stop()

    # 1. Roztriedenie podľa prvého písmena MKCH-10
    df["mkch10_prve_pismeno"] = df["diagnoza MKCH-10"].str[0]

    skupiny = {
        "K": "Gastro a pečeň",
        "E": "Metabolické",
        "D": "Hematológia",
        "B": "Infekcie",
        "C": "Novotvary",
        "Z": "Iné zdravotné faktory"
    }

    df["skupina_diag"] = df["mkch10_prve_pismeno"].map(skupiny).fillna("Ostatné")

    # 2. Vývoj v čase (extrahovanie roka z prijem_vzorky)
    try:
        df["rok_prijmu"] = pd.to_datetime(df["prijem_vzorky"], errors="coerce").dt.year
    except:
        st.warning("Nebol nájdený stĺpec 'prijem_vzorky' alebo je chybný formát.")
        df["rok_prijmu"] = None

    # Zobrazenie počtu diagnóz podľa rokov
    st.subheader("📅 Vývoj počtu diagnostikovaných pacientov v čase")
    if df["rok_prijmu"].notna().any():
        diag_time = df.groupby(["rok_prijmu", "skupina_diag"]).size().unstack(fill_value=0)
        st.line_chart(diag_time)
    else:
        st.warning("Chýbajú informácie o dátumoch prijatia vzoriek.")

    # 3. Zobrazenie tabuľky diagnóz a skupín
    st.subheader("📋 Diagnózy a ich zaradenie do skupín")
    st.dataframe(df[["id", "diagnoza MKCH-10", "skupina_diag", "rok_prijmu"]])

    # 4. Skontrolovanie zastaraných kódov (voliteľné)
    # Slovník pre zastarané kódy
    zastarale_kody_info = {
        "K76.0": {
            "novy_kod": "E66.9 (Obezita, nešpecifikovaná)",  # Príklad
            "popis": "K76.0 označuje tukovú degeneráciu pečene. V modernej klasifikácii sa rozlišuje nealkoholová tuková choroba pečene (NAFLD)."
        },
        "K75.9": {
            "novy_kod": "K75.8 (Iné zápalové ochorenia pečene)",  # Príklad
            "popis": "K75.9 je nešpecifikovaná zápalová choroba pečene. Dnes sa odporúča podrobnejšie klasifikovať zápal podľa príčiny."
        }
    }

    st.subheader("Zastarané diagnózy podľa MKCH-10")

    df["zastara_diag"] = df["diagnoza MKCH-10"].isin(zastarale_kody_info.keys())

    if df["zastara_diag"].any():
        for kod, info in zastarale_kody_info.items():
            pocet = (df["diagnoza MKCH-10"] == kod).sum()
            if pocet > 0:
                st.markdown(f"**Kód {kod}**")
                st.markdown(f"- Počet pacientov: **{pocet}**")
                st.markdown(f"- Odporúčaný nový kód: **{info['novy_kod']}**")
                st.markdown(f"- Popis: {info['popis']}")
                st.markdown("---")
    else:
        st.success("✅ V datasete sa nenachádzajú žiadne zastarané MKCH-10 kódy.")

# # ===== Export sekcia v Sidebare =====
#     st.sidebar.header("📤 Export reportu")

#     if st.sidebar.button("Exportovať report (.docx)"):
#         with st.spinner("Generuje sareport..."):
#             report_path = export_report.generate_report(df)
#             with open(report_path, "rb") as f:
#                 st.sidebar.download_button(
#                     label="Stiahni report",
#                     data=f,
#                     file_name="HFE_gene_analysis_report.docx",
#                     mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document"
#                 )
# 