import pandas as pd

# Nacitanie CSV so spravnym oddelovacom
df = pd.read_csv("SSBU25_dataset.csv", sep=";", engine="python", encoding="utf-8-sig", dtype=str)

# Spojenie datum+cas pre 'validovany vysledok' a 'prijem vzorky'
df['validovany_vysledok'] = df.iloc[:, 1].fillna('') + ' ' + df.iloc[:, 2].fillna('')
df['prijem_vzorky'] = df.iloc[:, 3].fillna('') + ' ' + df.iloc[:, 4].fillna('')

# Vymazanie povodnych datum+cas stlpcov + unnamed stlpce
stlpce_na_mazanie = ['Unnamed: 3', 'Unnamed: 5'] + [col for col in df.columns if "Unnamed" in col]
df = df.drop(columns=stlpce_na_mazanie + [df.columns[1], df.columns[2], df.columns[3], df.columns[4]])

# Odstranenie riadkov s prazdnym ID
df = df[df['id'].notna() & (df['id'].str.strip() != "")]

# Validacia veku: kontrolujeme len prvu zlozku (napr. 59 z "59,70")
def is_valid_age(age_str):
    try:
        age_main = float(age_str.split(',')[0].replace(',', '.'))
        return age_main <= 120
    except:
        return False

df = df[df['vek'].apply(is_valid_age)]

# Ulozenie upraveneho datasetu
df.to_csv("SSBU25_dataset_cleaned.csv", index=False, sep=";", encoding="utf-8-sig")
