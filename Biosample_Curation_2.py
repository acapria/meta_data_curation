from Bio import Entrez, SeqIO
import csv
import pandas as pd
from urllib.error import HTTPError
import xml.etree.ElementTree as ET


# extract the accessions from the genome group csv file
def extract_accessions(file_path):
    df = pd.read_csv(file_path)
    accession_numbers = df["GenBank Accessions"].tolist()
    return accession_numbers


# fetch the genbank record
def fetch_genbank_records(accession_numbers):
    Entrez.email = "bvbrcacapria@gmail.com"  # put your email
    records = []
    for accession_number in accession_numbers:
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            records.append(record)
            handle.close()
        except Exception as e:
            print(f"Error fetching GenBank record for accession {accession_number}: {e}")
    return records


# extract the host and labhost metadata from a genbank record
def extract_host_info(record):
    host_info = ""
    lab_host_info = ""
    for feature in record.features:
        if feature.type == "source":
            qualifiers = feature.qualifiers
            if "host" in qualifiers:
                host_info = ', '.join(qualifiers["host"])
            if "lab_host" in qualifiers:
                lab_host_info = ', '.join(qualifiers["lab_host"])
    return host_info, lab_host_info


def handle_missing_info(accession_number):
    all_biosample_accessions = []
    try:
        handle = Entrez.elink(dbfrom="nucleotide", db="biosample", id=accession_number)
        record = Entrez.read(handle)
        handle.close()
        biosample_accessions = []
        if record and record[0].get('LinkSetDb', []):
            links = record[0]['LinkSetDb'][0].get('Link', [])
            biosample_accessions.extend(link['Id'] for link in links)
        if biosample_accessions:  # Only append if BioSample accessions exist
            all_biosample_accessions.extend(biosample_accessions)
            # print(all_biosample_accessions)
    except HTTPError as e:
        print(f"HTTP Error occurred for GenBank accession {accession_number}: {e}")
    except Exception as e:
        print(f"An error occurred for GenBank accession {accession_number}: {e}")
    return all_biosample_accessions


def get_biosample_record(biosample_accession):
    try:
        handle = Entrez.efetch(db="biosample", id=biosample_accession, retmode="xml")
        xml_data = handle.read()
        handle.close()
        return xml_data
    except HTTPError as e:
        print(f"HTTP Error occurred for BioSample accession {biosample_accession}: {e}")
        return None


def extract_information(xml_data, genbank_accession):
    if xml_data is None:
        return None
    root = ET.fromstring(xml_data)
    sample_attributes = root.findall(".//Attribute[@harmonized_name]")
    extracted_info = {"GenBank Accession": genbank_accession}
    print(xml_data)
    for attr in sample_attributes:
        harmonized_name = attr.attrib.get('harmonized_name')
        # Only interested in the following attribute names:
        # Attribute names for host: host common name or common name
        # Attribute names for lab host: lab_host or lab host
        # Additional fields to add to excel: host scientific name or scientific_name
        if 'host' in harmonized_name or 'isolation_source' in harmonized_name or 'lab_host' in harmonized_name or 'host_common_name' in harmonized_name:
            # print(attribute_name)
            attribute_value = attr.text
            # print(attribute_value)
            extracted_info[harmonized_name] = attribute_value
    # print(extracted_info)
    return extracted_info


import csv


def main(extract_biosample_info=False):
    file_path = "/Users/rbhattac/Desktop/BVBRC_genome-6.csv"  # path to your genome group csv file
    # extract the accession numbers
    accession_numbers = extract_accessions(file_path)
    # accession_numbers=['KT844544']

    # Fetch GenBank records
    genbank_records = fetch_genbank_records(accession_numbers)

    # open the output file to write
    with open("/Users/rbhattac/Desktop/results_code.csv", "w", newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        # Write header
        csv_writer.writerow(
            ["Accession", "Genbank Host", "Genbank Lab Host", "Biosample Host Common Name", "Biosample Lab Host",
             "Biosample Host Scientific Name", "Biosample Isolation Source"])
        for record in genbank_records:
            accession = record.id
            # Extract host and lab host information
            host_info, lab_host_info = extract_host_info(record)  # Extracted directly from GenBank record
            # Check if both GenBank host and lab host info exist and if user wants to extract BioSample info
            if (not host_info or not lab_host_info or extract_biosample_info):
                biosample_accessions = handle_missing_info(accession)
                biosample_host_info = ''  # Initialize host_info
                biosample_lab_host_info = ''  # Initialize lab_host_info
                biosample_isolation_source = ''
                biosample_scientific_names = []  # Initialize biosample_scientific_names
                # Loop through each biosample accession
                for biosample_accession in biosample_accessions:
                    biosample_record = get_biosample_record(biosample_accession)
                    if biosample_record:
                        extracted_info = extract_information(biosample_record, accession)
                        if extracted_info:
                            # Check if attribute is host common name, common name
                            if 'host_common_name' in extracted_info:
                                biosample_host_info = extracted_info['host_common_name']
                            # elif 'common name' in extracted_info:
                            #   biosample_host_info = extracted_info['common name']
                            if 'isolation_source' in extracted_info:
                                biosample_isolation_source = (extracted_info['isolation_source'])
                            # elif 'isolation_source' in extracted_info:
                            #   biosample_isolation_source=(extracted_info['isolation_source'])
                            # Check if attribute is lab_host
                            if 'lab_host' in extracted_info:
                                biosample_lab_host_info = extracted_info['lab_host']
                            # elif 'laboratory_host' in extracted_info:
                            # biosample_lab_host_info = extracted_info['laboratory_host']
                            # Check if attribute is host scientific name or scientific_name
                            if 'host' in extracted_info:
                                biosample_scientific_names.append(extracted_info['host'])
                            # elif 'scientific_name' in extracted_info:
                            #   biosample_scientific_names.append(extracted_info['scientific_name'])
                            # elif 'host' in extracted_info:
                            # biosample_scientific_names.append(extracted_info['host'])
                            # elif 'specific_host' in extracted_info:
                            # biosample_scientific_names.append(extracted_info['specific_host'])

                        else:
                            print(f"No information extracted for BioSample accession {biosample_accession}")
                    else:
                        print(f"No record found for BioSample accession {biosample_accession}")
                # Write data to CSV
                csv_writer.writerow([accession, host_info, lab_host_info, biosample_host_info, biosample_lab_host_info,
                                     ", ".join(biosample_scientific_names), biosample_isolation_source])
            else:
                # Write data to CSV if host and lab host info are already available
                csv_writer.writerow([accession, host_info, lab_host_info, "", "", ""])


if __name__ == "__main__":
    extract_biosample_info = input(
        "Do you want to extract information from BioSample even if both GenBank host and lab host exist? (yes/no): ").lower()
    if extract_biosample_info == "yes":
        main(True)
    else:
        main(False)
