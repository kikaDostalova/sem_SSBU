import os
import pandas as pd
from docx import Document
from docx.shared import Inches

def generate_report():
    os.makedirs("report", exist_ok=True)

    doc = Document()
    doc.add_heading('Analýza HFE génu a genetická predispozícia na hemochromatózu', 0)

    doc.add_paragraph("Semestrálna práca – analytická časť.\n")
    doc.add_page_break()

    # Obsah
    doc.add_heading('Obsah', level=1)
    doc.add_paragraph("1. Základné informácie o datasete\n2. Grafy\n3. Hardy-Weinbergove testy\n4. Percentá genotypov\n5. Súvislosť mutácií s pečeňovými diagnózami\n6. Diagnózy MKCH-10")
    doc.add_page_break()

    # 1. Základné informácie o datasete
    doc.add_heading('1. Základné informácie o datasete', level=1)
    df = pd.read_csv("sem_SSBU/SSBU25_dataset_cleaned.csv", sep=";", encoding="utf-8-sig")
    doc.add_paragraph(f"Počet riadkov: {df.shape[0]}")
    doc.add_paragraph(f"Počet stĺpcov: {df.shape[1]}")
    doc.add_paragraph("Stĺpce: " + ", ".join(df.columns))
    doc.add_page_break()

    # 2. Grafy
    doc.add_heading('2. Grafy', level=1)
    graf_folder = "grafy"
    if os.path.exists(graf_folder):
        for img in sorted(os.listdir(graf_folder)):
            if img.endswith(".png"):
                doc.add_paragraph(img.replace(".png", "").replace("_", " ").capitalize())
                doc.add_picture(os.path.join(graf_folder, img), width=Inches(5))
    doc.add_page_break()

    # 3. Hardy-Weinberg testy
    doc.add_heading('3. Hardy-Weinbergove testy', level=1)
    hwe_path = "tabulky/hardy_weinberg_test.csv"
    if os.path.exists(hwe_path):
        hwe_df = pd.read_csv(hwe_path, sep=";")
        for idx, row in hwe_df.iterrows():
            doc.add_paragraph(f"Mutácia {row['Mutácia']} - p-hodnota: {row['p-hodnota']} - Výsledok: {row['Výsledok']}")
    doc.add_page_break()

    # 4. Percentá genotypov
    doc.add_heading('4. Percentuálne rozdelenie genotypov', level=1)
    percenta_path = "tabulky/percenta_genotypov.csv"
    if os.path.exists(percenta_path):
        percenta_df = pd.read_csv(percenta_path, sep=";")
        for idx, row in percenta_df.iterrows():
            doc.add_paragraph(f"{row['Mutácia']} - {row['Genotyp']}: {row['Počet pacientov']} pacientov ({row['Percentá']}%)")
    doc.add_page_break()

    # 5. Súvislosť s pečeňovými diagnózami
    doc.add_heading('5. Súvislosť HFE mutácií a pečeňových diagnóz', level=1)
    suvislost_path = "tabulky/suvislost_mutacie_pecen.csv"
    if os.path.exists(suvislost_path):
        suvislost_df = pd.read_csv(suvislost_path, sep=";")
        for idx, row in suvislost_df.iterrows():
            doc.add_paragraph(f"Mutácia {row['Mutácia']} - p-hodnota: {row['p-hodnota']} - Výsledok: {row['Výsledok']}")
    doc.add_page_break()

    # 6. Vývoj diagnóz MKCH-10
    doc.add_heading('6. Diagnózy MKCH-10 a vývoj v čase', level=1)
    mkch_graf = "grafy/graf_diag_vyvoj_v_case.png"
    if os.path.exists(mkch_graf):
        doc.add_picture(mkch_graf, width=Inches(5))

    # Uloženie
    report_path = "report/Analyza_HFE_gen.docx"
    doc.save(report_path)
    return report_path
