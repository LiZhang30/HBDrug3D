import requests
import json
import time
import random
from bs4 import BeautifulSoup


def search_html(params="", page=0):
    headers = {
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7",
        "Accept-Language": "zh-CN,zh;q=0.9",
        "Cache-Control": "max-age=0",
        "Proxy-Connection": "keep-alive",
        "Upgrade-Insecure-Requests": "1",
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126.0.0.0 Safari/537.36"
    }
    url = "http://adcdb.idrblab.net/search/result/adc"
    params = {
        "search_api_fulltext": params,
        "page": page
    }
    response = requests.get(url, headers=headers, params=params, verify=False)
    return response.text


def parse_html_to_json(html_content=""):
    if (html_content==""):
        return []
    soup = BeautifulSoup(html_content, 'html.parser')
    parent_div = soup.find('div', class_='col-md-12')
    units = parent_div.find_all('div', class_='div-unit-search')
    results = []

    for unit in units:
        adc_id = unit.select_one('span:contains("ADC ID:")').text.split(': ')[1].strip() if unit.select_one('span:contains("ADC ID:")') else ""
        adc_name = unit.select_one('span:contains("ADC Name:")').contents[-1].strip() if unit.select_one('span:contains("ADC Name:")') else ""
        adc_info_link = unit.select_one('a:contains("ADC Info")')['href'] if unit.select_one('a:contains("ADC Info")') else ""
        drug_status = unit.select_one('div:contains("Drug Status:")').text.split(': ')[1].strip() if unit.select_one('div:contains("Drug Status:")') else ""
        representative_indication = unit.select_one('div:contains("Representative Indication:")').text.split(': ')[1].strip() if unit.select_one('div:contains("Representative Indication:")') else ""
        antibody_name = unit.select_one('span:contains("Antibody Name:")').contents[-1].strip() if unit.select_one('span:contains("Antibody Name:")') else ""
        antibody_info_link = unit.select_one('a:contains("Antibody Info")')['href'] if unit.select_one('a:contains("Antibody Info")') else ""
        payload_name = unit.select_one('span:contains("Payload Name:")').contents[-1].strip() if unit.select_one('span:contains("Payload Name:")') else ""
        payload_info_link = unit.select_one('a:contains("Payload Info")')['href'] if unit.select_one('a:contains("Payload Info")') else ""
        linker_name = unit.select_one('span:contains("Linker Name:")').contents[-1].strip() if unit.select_one('span:contains("Linker Name:")') else ""
        linker_info_link = unit.select_one('a:contains("Linker Info")')['href'] if unit.select_one('a:contains("Linker Info")') else ""

        data = {
            "adc_id": adc_id,
            "adc_name": adc_name,
            "adc_info_link": adc_info_link,
            "drug_status": drug_status,
            "representative_indication": representative_indication,
            "antibody_name": antibody_name,
            "antibody_info_link": antibody_info_link,
            "payload_name": payload_name,
            "payload_info_link": payload_info_link,
            "linker_name": linker_name,
            "linker_info_link": linker_info_link
        }

        results.append(data)
    return results


def load_existing_data(file_path):
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        return []

def save_data(file_path, data):
    with open(file_path, 'w') as f:
        json.dump(data, f, indent=4)

if __name__ == '__main__':
    file_path = 'data/adc_data.json'
    existing_data = load_existing_data(file_path)
    existing_adc_ids = {entry['adc_id'] for entry in existing_data}
    loop = 'abcdefghijklmnopqrstuvwxyz0123456789'

    for char in range(0, loop.__len__()):
        search_char = loop[char]
        page = 0

        while True:
            html_content = search_html(params=search_char, page=page)
            results = parse_html_to_json(html_content)

            if not results:
                break

            new_results = [result for result in results if result['adc_id'] not in existing_adc_ids]

            if new_results:
                existing_data.extend(new_results)
                existing_adc_ids.update(result['adc_id'] for result in new_results)
                save_data(file_path, existing_data)
                print(f"Saved {len(new_results)} new results for {search_char} page {page}")

            page += 1
            time.sleep(random.uniform(0, 2))