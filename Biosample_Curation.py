#write a function to access biosample record from NCBI

from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd

# Provide your email to NCBI Entrez
Entrez.email = "abrilz2296@gmail.com"

def get_biosample_accession(genbank_accessions):
    handle = Entrez.elink(dbfrom="nucleotide", db="biosample", id=genbank_accession)
    record = Entrez.read(handle)
    #print(record)
    handle.close()
    biosample_accessions = []
    if record[0]['LinkSetDb']:
        links = record[0]['LinkSetDb'][0]['Link']
        biosample_accessions = [link['Id'] for link in links]
    return biosample_accessions


def get_biosample_record(biosample_accession):
    handle = Entrez.efetch(db="biosample", id=biosample_accession, retmode="xml")
    xml_data = handle.read()
    handle.close()
    return xml_data

def extract_information(xml_data, genbank_accession):
    # Parse the XML data
    root = ET.fromstring(xml_data)

    # Find the relevant sample attributes
    sample_attributes = root.findall(".//Attribute[@attribute_name]")

    # Extract information from the sample attributes
    extracted_info = {"GenBank Accession": genbank_accession}  # Include GenBank accession number
    for attr in sample_attributes:
        attribute_name = attr.attrib.get('attribute_name')
        attribute_value = attr.text
        extracted_info[attribute_name] = attribute_value

    return extracted_info

genbank_accessions = ['K02025','KP849471','OY757791']
#excel_file = "BVBRC_genome-4.xlsx"
#df = pd.read_excel(excel_file)
#genbank_accessions = df['GenBank Accessions'].tolist()

all_extracted_info = []

print("Fetching BioSample information...")

for i, genbank_accession in enumerate(genbank_accessions, 1):
    print(f"Processing GenBank accession {genbank_accession} ({i}/{len(genbank_accessions)})...")
    try:
        biosample_accessions = get_biosample_accession([genbank_accession])
        for accession in biosample_accessions:
            biosample_xml = get_biosample_record(accession)
            extracted_info = extract_information(biosample_xml, genbank_accession)
            all_extracted_info.append(extracted_info)
        print(f"BioSample information retrieved for GenBank accession {genbank_accession}.")
    except Exception as e:
        print(f"Error processing GenBank accession {genbank_accession}: {e}")

print("Creating DataFrame...")
df = pd.DataFrame(all_extracted_info)

# Move 'GenBank Accession' column to the first position (column 'A')
#genbank_column = df.pop('GenBank Accession')  # Remove 'GenBank Accession' column
#df.insert(0, 'GenBank Accession', genbank_column)  # Insert 'GenBank Accession' column at the beginning

# Write DataFrame to Excel
excel_filename = "biosample_information.xlsx"
df.to_excel(excel_filename, index=False)

print(f"Results saved to {excel_filename}")