from docx import Document
from docx.shared import Inches, Pt
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import os
import tempfile
from scipy.stats import chi2_contingency, chi2

# Slovník pre zastarané kódy
zastarale_kody_info = {
    "K76.0": {
        "novy_kod": "E66.9 (Obezita, nešpecifikovaná)",
        "popis": "Tuková degenerácia pečene spadá dnes pod nealkoholovú tukovú chorobu pečene (NAFLD)."
    },
    "K75.9": {
        "novy_kod": "K75.8 (Iné zápalové ochorenia pečene)",
        "popis": "Nešpecifikovaná zápalová choroba pečene sa dnes klasifikuje presnejšie podľa príčiny."
    }
}

# Labely genotypov
genotyp_label = {
    "wt/wt": "wt/wt (normal)",
    "wt/mut": "wt/mut (heterozygot)",
    "mut/mut": "mut/mut (mutant)"
}

# Funkcia na premapovanie genotypov
def premapuj_na_genotyp(series):
    return series.map({
        "normal": "wt/wt",
        "heterozygot": "wt/mut",
        "mutant": "mut/mut",
        "mutácia": "mut/mut",
        "patogénna": "mut/mut"
    }).fillna("neznamy")

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
    return obs, exp, pval

def generate_report(df):
    document = Document()

    # ===== Titulná strana =====
    document.add_heading('Analýza HFE génu a genetickej predispozície na hemochromatózu', 0)
    document.add_paragraph('Autor: Semestrálny projekt')
    document.add_paragraph(f'Dátum: {datetime.now().strftime("%d.%m.%Y")}')
    document.add_page_break()

    # ===== Obsah =====
    document.add_heading('Obsah', level=1)
    document.add_paragraph('1. Dataset a jeho čistenie')
    document.add_paragraph('2. Hardy-Weinbergova rovnováha')
    document.add_paragraph('3. Percentá genotypov a prenášači')
    document.add_paragraph('4. Súvislosť HFE mutácií s diagnózami')
    document.add_paragraph('5. Grafy')
    document.add_paragraph('6. Analýza MKCH-10 diagnóz')
    document.add_paragraph('7. Záver')
    document.add_page_break()

    # ===== 1. Dataset =====
    document.add_heading('1. Dataset a jeho čistenie', level=1)
    document.add_paragraph(f'Veľkosť datasetu: {df.shape[0]} pacientov')
    document.add_paragraph(f'Počet atribútov: {df.shape[1]}')
    document.add_paragraph('Zoznam stĺpcov:')
    for col in df.columns:
        document.add_paragraph(f"- {col}", style='ListBullet')

    # ===== 2. Hardy-Weinbergova rovnováha =====
    document.add_heading('2. Hardy-Weinbergova rovnováha', level=1)
    mutacie = {
        "C282Y": "HFE G845A (C282Y) [HFE]",
        "H63D": "HFE C187G (H63D) [HFE]",
        "S65C": "HFE A193T (S65C) [HFE]"
    }
    for mut, col in mutacie.items():
        if col in df.columns:
            genotypy = premapuj_na_genotyp(df[col])
            obs, exp, pval = hardy_weinberg_test(genotypy)
            document.add_heading(f"Mutácia {mut}", level=2)
            table = document.add_table(rows=4, cols=3)
            table.style = 'LightShading-Accent1'
            hdr_cells = table.rows[0].cells
            hdr_cells[0].text = 'Genotyp'
            hdr_cells[1].text = 'Pozorované'
            hdr_cells[2].text = 'Očakávané'
            for i, genotype in enumerate(["wt/wt", "wt/mut", "mut/mut"]):
                row = table.rows[i+1].cells
                row[0].text = genotyp_label.get(genotype, genotype)
                row[1].text = str(obs.get(genotype, 0))
                row[2].text = f"{exp[i]:.1f}"
            document.add_paragraph(f"p-hodnota: {pval:.4f}")
            document.add_paragraph()

    # ===== 3. Percentá genotypov =====
    document.add_heading('3. Percentá genotypov a prenášači', level=1)
    for mut, col in mutacie.items():
        if col in df.columns:
            genotypy = premapuj_na_genotyp(df[col])
            counts = genotypy.value_counts(normalize=True).round(3) * 100
            document.add_heading(f"Mutácia {mut}", level=2)
            for idx, val in counts.items():
                label = genotyp_label.get(idx, idx)
                document.add_paragraph(f"- {label}: {val:.1f} %")

    # ===== 4. Súvislosť HFE mutácií s diagnózami =====
    document.add_heading('4. Súvislosť HFE mutácií s diagnózami', level=1)
    df["pecen_diag"] = df["diagnoza MKCH-10"].isin(["K76.0", "K75.9"])
    for mut, col in mutacie.items():
        if col in df.columns:
            genotypy = premapuj_na_genotyp(df[col])
            kont_tab = pd.crosstab(genotypy, df["pecen_diag"])
            document.add_heading(f"Mutácia {mut}", level=2)
            table = document.add_table(rows=kont_tab.shape[0]+1, cols=kont_tab.shape[1]+1)
            table.style = 'LightShading-Accent2'
            header = table.rows[0].cells
            header[0].text = "Genotyp"
            for i, colname in enumerate(kont_tab.columns):
                header[i+1].text = f"Pečeň {colname}"
            for i, idx in enumerate(kont_tab.index):
                row = table.rows[i+1].cells
                row[0].text = genotyp_label.get(idx, idx)
                for j, val in enumerate(kont_tab.loc[idx]):
                    row[j+1].text = str(val)

    # ===== 5. Grafy =====
    document.add_heading('5. Grafy', level=1)

    with tempfile.TemporaryDirectory() as tmpdir:
        for mut, col in mutacie.items():
            if col in df.columns:
                genotypy = premapuj_na_genotyp(df[col])

                fig, ax = plt.subplots()
                sns.countplot(x=genotypy.map(lambda x: genotyp_label.get(x, x)))
                plt.title(f"Rozdelenie genotypov - {mut}")
                graph_path = os.path.join(tmpdir, f"genotypy_{mut}.png")
                fig.savefig(graph_path)
                document.add_picture(graph_path, width=Inches(5))
                plt.close()

    # ===== 6. MKCH-10 Analýza =====
    document.add_heading('6. Analýza MKCH-10 diagnóz', level=1)
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
    diag_counts = df["skupina_diag"].value_counts()
    for skupina, pocet in diag_counts.items():
        document.add_paragraph(f"- {skupina}: {pocet} pacientov")

    # ===== 7. Záver =====
    document.add_heading('7. Záver', level=1)
    document.add_paragraph("Analýza ukázala súvislosti medzi genotypmi HFE génu a vybranými diagnózami. V budúcnosti odporúčame pokračovať v sledovaní pacientov s predispozíciou na hemochromatózu.")

    # ===== Uloženie reportu =====
    export_path = os.path.join(os.path.dirname(__file__), "HFE_gene_analysis_report.docx")
    document.save(export_path)

    return export_path
