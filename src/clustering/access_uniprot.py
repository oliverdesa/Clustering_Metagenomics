import requests
from xml.etree import ElementTree
import xml.etree.ElementTree as ET

uniprot_ids = ["A0A0S2HUS6", "D6YTN7", "A0A1X9MM83"]

def extract_pdb_ids(xml_content):
    # Parse the XML content
    root = ET.fromstring(xml_content)

    # Define the namespace used in UniProt XML
    ns = {'uniprot': 'http://uniprot.org/uniprot'}

    # Find all 'dbReference' elements with type 'PDB'
    pdb_elements = root.findall(".//{http://uniprot.org/uniprot}dbReference[@type='PDB']")

    # Extract the 'id' attribute from each 'dbReference' element
    pdb_ids = [element.attrib['id'] for element in pdb_elements]

    return pdb_ids

with open('C:\\Users\\odesa\\Desktop\\CRCFinal\\CRC-Final\\D6YTN7.xml', 'r') as file:
    xml_content = file.read()

    pdb_ids = extract_pdb_ids(xml_content)

    print(pdb_ids)


def test_uniprot_request(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
    response = requests.get(url)

    # Parse the XML response
    tree = ET.fromstring(response.content)

    # Extract PDB IDs using the extract_pdb_ids function
    pdb_ids = extract_pdb_ids(tree)

    print(f"PDB IDs for UniProt ID {uniprot_id}: {pdb_ids}")

def download_uniprot_xml(uniprot_id, file_name):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
    response = requests.get(url)

    if response.status_code == 200:
        with open(file_name, 'w') as file:
            file.write(response.text)
        print(f"XML file for UniProt ID {uniprot_id} has been saved as {file_name}")
    else:
        print(f"Failed to download XML for UniProt ID {uniprot_id}. HTTP Status Code: {response.status_code}")

# download_uniprot_xml("D6YTN7", "D6YTN7.xml")


# for uniprot_id in uniprot_ids:
#     response = requests.get(f"https://www.uniprot.org/uniprot/{uniprot_id}.xml")
#     tree = ElementTree.fromstring(response.content)
#     # Extract PDB IDs from the XML (depends on XML structure)
    
#     print(tree)

    # pdb_ids = extract_pdb_ids(tree)

    # for pdb_id in pdb_ids:
    #     pdb_response = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", stream=True)
    #     with open(f"{pdb_id}.pdb", 'wb') as file:
    #         for chunk in pdb_response.iter_content(chunk_size=128):
    #             file.write(chunk)
