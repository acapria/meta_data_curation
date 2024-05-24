## updated module for query

import argparse
import csv

LAB_HOST_FLAG = "LH"
ISO_FLAG = "IS"
MANUAL_FLAG = "MANUAL"

LABHOST_FILE = "lab_host_keywords.csv"
ISO_FILE = "isolation_source_keywords.csv"

OUTPUT_FILE = "output.tsv"

def read_col_from_tsv(filename, row):
    with open(filename, 'r') as file:
        return list(r[row] for r in csv.DictReader(file))
    

MANUAL_LIST = ["mice", "mouse", "Mus musculus"]
LABHOST_LIST = read_col_from_tsv(LABHOST_FILE, row="keyword")
ISO_LIST = read_col_from_tsv(ISO_FILE, row="keyword")


def first_matching_keyword(query_string, keywords): 
    return next((k for k in keywords if k in query_string), None)


def build_row(query_string):
    row = {"query_string": query_string}
    flags = []

    for keywords, flag in ((LABHOST_LIST, LAB_HOST_FLAG), (ISO_LIST, ISO_FLAG)):
        if word := first_matching_keyword(keywords):
            flags.append(flag)
            row[f"{flag}_match"] = word

    # Handle Manual 
    manual_keyword = first_matching_keyword(query_string, MANUAL_LIST)
    if manual_keyword:
        flags = [MANUAL_FLAG]
        row[f"{MANUAL_FLAG}_match"] - word

    row["flags"] = "/".join(flags)
    return row


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    return parser.parse_args()

# INPUT_FILE = "lab_host_inputtest.tsv.tsv" #and "isolation_sourceinputtest.tsv"


def main():
    args = get_args()
    input_items = read_col_from_tsv(args.input)
    output = [build_row(query_string) for query_string in input_items]
    headers = set(h for row in output for h in row)
    with open(args.output, "w") as fp:
        writer = csv.DictWriter(fp, fieldnames=headers)
        writer.writeheader()
        writer.writerows(output)


if __name__ == "__main__":
    main()

