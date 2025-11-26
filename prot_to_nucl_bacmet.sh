#!/bin/bash
# Extract BacMet nucleotide sequences from public databases

python3 <<'EOF'
import subprocess
import pandas as pd

dbfetch = "/storage/bioinfo/webservice-clients/python/dbfetch.py"
interviews_df = pd.read_csv('BACMET_ID1.txt', sep='\t')
rows = []
for row in interviews_df.values:
    ID = row[0]
    jj = subprocess.check_output([dbfetch, 'fetchBatch', 'UniParc', ID]).decode("utf-8")
    for item in jj.split("\n"):
        if "EMBLWGS" in item:
            kkk = list(item.strip().split(" "))
            sstrmatch = [s for s in kkk if "id=" in s]
            emblwgs = sstrmatch[0].strip('"')
            rows.append([ID, emblwgs])
df = pd.DataFrame(rows, columns=["ID", "EMBL"])
df.to_csv('ID_EMBL.csv', index=False)

interviews_df = pd.read_csv('ID_EMBL.csv', sep=',')
for row in interviews_df.values:
    ID = row[1]
    jj = subprocess.check_output([dbfetch, 'fetchBatch', 'ENA_Coding', ID, "fasta", "raw"]).decode("utf-8")
    with open("BACMET_PRED_FULL.txt", "a+") as file_object:
        file_object.write(jj)
EOF
