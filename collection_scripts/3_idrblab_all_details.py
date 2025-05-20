import requests
from bs4 import BeautifulSoup
import json
import os
import time

def fetch_data(url=""):
    headers = {
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7",
        "Accept-Language": "zh-CN,zh;q=0.9",
        "Cache-Control": "max-age=0",
        "Proxy-Connection": "keep-alive",
        "Upgrade-Insecure-Requests": "1",
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126.0.0.0 Safari/537.36"
    }
    try:
        response = requests.get(url, headers=headers, verify=False)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Failed to fetch {url}: {e}", flush=True)
        return None

    soup = BeautifulSoup(response.text, 'html.parser')
    first_unit = soup.select_one('div.div-unit')
    if first_unit:
        return extract_table_data(first_unit)
    else:
        return None

def extract_table_data(soup):
    data = {}
    rows = soup.select('div.div-content table tbody tr')

    for row in rows:
        th_element = row.find('th')
        if th_element:
            th = th_element.get_text(strip=True)
            if th == "Pharmaceutical Properties":
                properties = {}
                td_elements = row.find_all('td')
                if len(td_elements) >= 4:
                    prop_name1 = td_elements[0].get_text(strip=True)
                    prop_value1 = td_elements[1].get_text(strip=True)
                    prop_name2 = td_elements[2].get_text(strip=True)
                    prop_value2 = td_elements[3].get_text(strip=True)
                    properties[prop_name1] = prop_value1
                    properties[prop_name2] = prop_value2

                prop_rows = row.find_next_siblings('tr')
                for prop_row in prop_rows:
                    td_elements = prop_row.find_all('td')
                    if len(td_elements) >= 4:
                        prop_name1 = td_elements[0].get_text(strip=True)
                        prop_value1 = td_elements[1].get_text(strip=True)
                        prop_name2 = td_elements[2].get_text(strip=True)
                        prop_value2 = td_elements[3].get_text(strip=True)
                        properties[prop_name1] = prop_value1
                        properties[prop_name2] = prop_value2
                    else:
                        break
                data[th] = properties

    for row in rows:
        th_element = row.find('th')
        if not th_element:
            continue

        th = th_element.get_text(strip=True)
        if 'ID' in th:
            td = row.find('td').get_text(strip=True) if row.find('td') else ""
            data[th] = td
        elif th == "Structure":
            continue
        elif th == "Pharmaceutical Properties":
            continue
        elif th == "Indication":
            indications = []
            td_element = row.find('td')
            if td_element:
                child_indications = td_element.find_all('div', class_='child-indication')
                for indication_div in child_indications:
                    spans = indication_div.find_all('span')
                    if len(spans) >= 2:
                        indications.append({
                            "indication": spans[0].get_text(strip=True),
                            "status": spans[1].get_text(strip=True)
                        })
            data[th] = indications
        else:
            td = row.find('td').get_text(strip=True) if row.find('td') else ""
            link = row.find('a')['href'] if row.find('a') else ""
            if link:
                data[th] = {"text": td, "link": link}
            else:
                data[th] = td

    return data

def load_data(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data

def save_data(file_path, data):
    with open(file_path, 'w') as f:
        json.dump(data, f, indent=4)

if __name__ == '__main__':
    domain = "http://adcdb.idrblab.net"
    input_file_path = "data/adc_data_cleaned.json"
    output_file_path = "data/adc_detail.json"
    batch_size = 10

    data_list = load_data(input_file_path)

    if os.path.exists(output_file_path):
        all_details = load_data(output_file_path)
    else:
        all_details = {}

    antibody_cache = {v["Antibody Name"]["text"]: v["antibody"] for k, v in all_details.items() if "antibody" in v and isinstance(v["antibody"].get("Antibody Name"), dict) and "text" in v["antibody"]["Antibody Name"]}
    payload_cache = {v["Name"]["text"]: v["payload"] for k, v in all_details.items() if "payload" in v and isinstance(v["payload"].get("Name"), dict) and "text" in v["payload"]["Name"]}
    linker_cache = {v["Linker Name"]["text"]: v["linker"] for k, v in all_details.items() if "linker" in v and isinstance(v["linker"].get("Linker Name"), dict) and "text" in v["linker"]["Linker Name"]}

    remaining_items = [item for item in data_list if item.get("adc_id") not in all_details]
    total_items = len(remaining_items)
    batch_data = {}

    for index, item in enumerate(remaining_items, start=1):
        adc_id = item.get("adc_id")
        if not adc_id:
            continue

        details = {"adc": {}, "antibody": {}, "payload": {}, "linker": {}, "search": item}
        url_mapping = {
            "adc": item.get("adc_info_link"),
            "antibody": item.get("antibody_info_link"),
            "payload": item.get("payload_info_link"),
            "linker": item.get("linker_info_link")
        }

        for key, url in url_mapping.items():
            if url:
                if key == "antibody" and item.get("antibody_name") and item.get("antibody_name") in antibody_cache:
                    details[key] = antibody_cache[item.get("antibody_name")]
                    print(f"Using cached data for Antibody Name {item.get('antibody_name')} in ADC ID {adc_id}", flush=True)
                elif key == "payload" and item.get("payload_name") and item.get("payload_name") in payload_cache:
                    details[key] = payload_cache[item.get("payload_name")]
                    print(f"Using cached data for Payload Name {item.get('payload_name')} in ADC ID {adc_id}", flush=True)
                elif key == "linker" and item.get("linker_name") and item.get("linker_name") in linker_cache:
                    details[key] = linker_cache[item.get("linker_name")]
                    print(f"Using cached data for Linker Name {item.get('linker_name')} in ADC ID {adc_id}", flush=True)
                else:
                    full_url = domain + url
                    fetched_data = fetch_data(full_url)
                    if fetched_data:
                        details[key] = fetched_data
                        if key == "antibody" and item.get("antibody_name"):
                            antibody_cache[item.get("antibody_name")] = fetched_data
                        elif key == "payload" and item.get("payload_name"):
                            payload_cache[item.get("payload_name")] = fetched_data
                        elif key == "linker" and item.get("linker_name"):
                            linker_cache[item.get("linker_name")] = fetched_data

        if details:
            batch_data[adc_id] = details
            print(f"Saved details for batch ending with ADC ID {adc_id} ({index}/{total_items})", flush=True)

        if index % batch_size == 0 or index == total_items:
            all_details.update(batch_data)
            save_data(output_file_path, all_details)
            batch_data.clear()

    print(f"Updated data has been saved to {output_file_path}", flush=True)
