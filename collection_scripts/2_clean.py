import json

def extract_drug_status_and_indication(text):
    drug_status = ""
    representative_indication = ""

    if "Drug Status:" in text and "Representative Indication:" in text:
        ds_index = text.find("Drug Status:")
        ri_index = text.find("Representative Indication:")
        if ds_index != -1 and ri_index != -1:
            drug_status_part = text[ds_index:ri_index].split("Drug Status:")[1].strip()
            drug_status = drug_status_part.split("\n")[0].strip()

            ri_part = text[ri_index:].split("Representative Indication:")[1].strip()
            representative_indication = ri_part.split("\n")[0].strip()

    return drug_status, representative_indication

def load_data(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data

def save_data(file_path, data):
    with open(file_path, 'w') as f:
        json.dump(data, f, indent=4)

def clean_data(data):
    for item in data:
        drug_status, representative_indication = extract_drug_status_and_indication(item["drug_status"])
        item["drug_status"] = drug_status
        item["representative_indication"] = representative_indication
    return data

if __name__ == '__main__':
    input_file_path = 'data/adc_data.json'
    output_file_path = 'data/adc_data_cleaned.json'

    data = load_data(input_file_path)
    cleaned_data = clean_data(data)
    save_data(output_file_path, cleaned_data)

    print(f"Cleaned data has been saved to {output_file_path}")
