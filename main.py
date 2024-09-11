import json
import requests
import tkinter as tk
from tkinter import messagebox, scrolledtext
from tkinter import ttk
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageTk
import io
from admet_ai import ADMETModel

model = ADMETModel()

def filter_results(results, non_important_features):
    filtered_results = []
    for result in results:
        filtered_result = {key: value for key, value in result.items() if key not in non_important_features}
        filtered_results.append(filtered_result)
    return filtered_results

def smiles_to_iupac(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "Неверный SMILES"
        iupac_name = Chem.MolToInchi(mol)
        return iupac_name
    except Exception as e:
        return "Ошибка при обработке SMILES"

def iupac_to_smiles(iupac_name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{iupac_name}/property/CanonicalSMILES/JSON"
        response = requests.get(url)
        data = response.json()
        smiles = data['PropertyTable']['Properties'][0]['CanonicalSMILES']
        return smiles
    except Exception as e:
        return None

def get_common_name(smiles):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/synonyms/JSON"
        response = requests.get(url)
        data = response.json()

        synonyms = data['InformationList']['Information'][0].get('Synonym', [])
        for synonym in synonyms:
            if len(synonym.split()) == 1:
                return synonym
        return None
    except Exception as e:
        return None

def get_pubchem_properties(smiles):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/property/MolecularWeight,LogP/JSON"
        response = requests.get(url)
        data = response.json()


        if 'PropertyTable' in data:
            properties = data['PropertyTable']['Properties'][0]
            return properties
        else:
            return None
    except Exception as e:
        return None

def generate_structure_image(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            img = Draw.MolToImage(mol)
            return img
        else:
            return None
    except Exception as e:
        return None


def search_compound(input_value, search_type):
    if search_type == "SMILES":
        return input_value
    elif search_type == "IUPAC":

        smiles = iupac_to_smiles(input_value)
        return smiles
    else:
        return None


def run_predictions():
    input_value = search_text.get("1.0", "end-1c")
    search_type = search_type_var.get()
    
    if not input_value:
        messagebox.showwarning("Ошибка", "Введите хотя бы один SMILES, IUPAC или название!")
        return


    processed_value = search_compound(input_value, search_type)
    if processed_value is None:
        messagebox.showwarning("Ошибка", "Не удалось найти соединение.")
        return


    pubchem_data = get_pubchem_properties(processed_value)
    common_name = get_common_name(processed_value)  # Ищем тривиальное название
    
    if pubchem_data:

        if common_name:
            searched_name.set(f"Найдено соединение: {common_name}")
        else:
            searched_name.set(f"Найдено соединение (IUPAC): {smiles_to_iupac(processed_value)}")
        
        show_structure_image(processed_value)
        display_pubchem_results(processed_value, pubchem_data)
    else:

        non_important_features = [
            "molecular_weight_drugbank_approved_percentile", "logP_drugbank_approved_percentile",
            "hydrogen_bond_acceptors_drugbank_approved_percentile", "hydrogen_bond_donors_drugbank_approved_percentile",
            "Lipinski_drugbank_approved_percentile", "QED_drugbank_approved_percentile",
            "Solubility_AqSolDB_drugbank_approved_percentile", "VDss_Lombardo_drugbank_approved_percentile"
        ]

        try:

            result = model.predict(smiles=processed_value)


            filtered_result = filter_results([result], non_important_features)
            
            searched_name.set(f"Соединение не найдено, предсказание для: {smiles_to_iupac(processed_value)}")
            show_structure_image(processed_value)
            display_predictions(processed_value, filtered_result[0])

        except Exception as e:
            messagebox.showerror("Ошибка", f"Ошибка во время предсказания: {str(e)}")


def display_pubchem_results(smiles, properties):
    result_text.config(state=tk.NORMAL)
    result_text.delete(1.0, tk.END)
    
    result_text.insert(tk.END, f"Соединение найдено в PubChem по SMILES: {smiles}\n", "header")
    result_text.insert(tk.END, "-"*50 + "\n")
    
    for key, value in properties.items():
        result_text.insert(tk.END, f"{key}: ", "property")
        result_text.insert(tk.END, f"{value}\n", "value")
    result_text.insert(tk.END, "-"*50 + "\n\n")
    
    result_text.config(state=tk.DISABLED)


def display_predictions(search_value, result):
    result_text.config(state=tk.NORMAL)
    result_text.delete(1.0, tk.END)
    
    result_text.insert(tk.END, f"Соединение не найдено в PubChem. Выполнено предсказание для: {search_value}\n", "header")
    result_text.insert(tk.END, "-"*50 + "\n")

    for key, value in result.items():
        result_text.insert(tk.END, f"{key}: ", "property")
        result_text.insert(tk.END, f"{value}\n", "value")
    result_text.insert(tk.END, "-"*50 + "\n\n")
    
    result_text.config(state=tk.DISABLED)


def show_structure_image(smiles):
    img = generate_structure_image(smiles)
    if img:

        img = img.resize((200, 200))
        img = ImageTk.PhotoImage(img)


        if structure_label.winfo_exists():
            structure_label.config(image='')

        structure_label.img = img
        structure_label.config(image=img)


root = tk.Tk()
root.title("ADMET Prediction")


search_label = tk.Label(root, text="Введите значение для поиска (SMILES, IUPAC, или название):", font=("Arial", 12))
search_label.pack(pady=5)

search_text = tk.Text(root, height=2, width=50, font=("Arial", 12))
search_text.pack(pady=5)


search_type_var = tk.StringVar(value="SMILES")
search_type_label = tk.Label(root, text="Выберите метод поиска:", font=("Arial", 12))
search_type_label.pack(pady=5)

search_type_combo = ttk.Combobox(root, textvariable=search_type_var, values=["SMILES", "IUPAC"], state="readonly", font=("Arial", 12))
search_type_combo.pack(pady=5)


predict_button = tk.Button(root, text="Запуск", command=run_predictions, font=("Arial", 12, "bold"), bg="#4CAF50", fg="white")
predict_button.pack(pady=10)


searched_name = tk.StringVar()
searched_name_label = tk.Label(root, textvariable=searched_name, font=("Arial", 12, "italic"))
searched_name_label.pack(pady=5)

structure_label = tk.Label(root)
structure_label.pack(pady=5)


result_text = scrolledtext.ScrolledText(root, height=20, width=80, font=("Arial", 10))
result_text.pack(pady=5)
result_text.config(state=tk.DISABLED)


exit_button = tk.Button(root, text="Выйти", command=root.quit, font=("Arial", 12), bg="#FF5733", fg="white")
exit_button.pack(pady=10)


result_text.tag_configure("header", font=("Arial", 12, "bold"), foreground="#4CAF50")
result_text.tag_configure("property", font=("Arial", 10, "bold"), foreground="#333333")
result_text.tag_configure("value", font=("Arial", 10, "underline"), foreground="#333333")


root.mainloop()
