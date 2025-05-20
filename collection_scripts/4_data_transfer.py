import pandas as pd
import json
from datetime import datetime
import os

USE_SAMPLE = False

with open('data/adc_detail.json', 'r') as file:
    data = json.load(file)

rows = []
columns = []

if USE_SAMPLE:
    sample_data = list(data.items())[:20]
else:
    sample_data = list(data.items())

for adc_id, details in sample_data:
    row = {}
    
    for key, value in details.get('search', {}).items():
        col_name = f'search.{key}'
        if col_name not in columns:
            columns.append(col_name)
        if isinstance(value, dict):
            row[col_name] = json.dumps(value)
        else:
            row[col_name] = value
    
    for key, value in details.get('adc', {}).items():
        col_name = f'adc.{key}'
        if col_name not in columns:
            columns.append(col_name)
        if isinstance(value, dict):
            row[col_name] = json.dumps(value)
        else:
            row[col_name] = value
    
    for key, value in details.get('antibody', {}).items():
        col_name = f'antibody.{key}'
        if col_name not in columns:
            columns.append(col_name)
        if isinstance(value, dict):
            row[col_name] = json.dumps(value)
        else:
            row[col_name] = value
    
    for key, value in details.get('payload', {}).items():
        col_name = f'payload.{key}'
        if col_name not in columns:
            columns.append(col_name)
        if isinstance(value, dict):
            row[col_name] = json.dumps(value)
        else:
            row[col_name] = value
    
    for key, value in details.get('linker', {}).items():
        col_name = f'linker.{key}'
        if col_name not in columns:
            columns.append(col_name)
        if isinstance(value, dict):
            row[col_name] = json.dumps(value)
        else:
            row[col_name] = value

    row['adc_id'] = adc_id
    rows.append(row)

df = pd.DataFrame(rows, columns=columns + ['adc_id'])

search_df = df[['adc_id'] + [col for col in df.columns if col.startswith('search.')]]
adc_df = df[['adc_id'] + [col for col in df.columns if col.startswith('adc.')]]
antibody_df = df[['adc_id'] + [col for col in df.columns if col.startswith('antibody.')]]
payload_df = df[['adc_id'] + [col for col in df.columns if col.startswith('payload.')]]
linker_df = df[['adc_id'] + [col for col in df.columns if col.startswith('linker.')]]

if USE_SAMPLE:
    output_dir = 'data/final/sample/'
else:
    output_dir = 'data/final/'

os.makedirs(output_dir, exist_ok=True)

timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
output_filename = os.path.join(output_dir, f"adc_details_{timestamp}.xlsx")

with pd.ExcelWriter(output_filename, engine='openpyxl') as writer:
    df.to_excel(writer, index=False, sheet_name='Merged')
    search_df.to_excel(writer, index=False, sheet_name='Search')
    adc_df.to_excel(writer, index=False, sheet_name='ADC')
    antibody_df.to_excel(writer, index=False, sheet_name='Antibody')
    payload_df.to_excel(writer, index=False, sheet_name='Payload')
    linker_df.to_excel(writer, index=False, sheet_name='Linker')

print(f"save as: {output_filename}")
