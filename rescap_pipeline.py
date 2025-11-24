#!/usr/bin/env python3
"""
Genomic Data Pipeline Script (2024 Edition)
Executes: Fastp QC → Bowtie2 Alignment → SeqKit Stats → KMA Mapping → SeqKit Normalization
Language: English
Author: Miguel
"""

import subprocess
import sys
import logging
from pathlib import Path
from datetime import datetime

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f'pipeline_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
        logging.StreamHandler(sys.stdout)
    ]
)

class GenomicPipeline:
    def __init__(self, threads=24):
        self.threads = threads

    def run_cmd(self, cmd, description):
        logging.info(f"Running: {description}")
        try:
            subprocess.run(cmd, shell=True, check=True)
            logging.info(f"Completed: {description}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error in: {description}\n{e}")
            sys.exit(1)

    def fastp_qc(self, fq_dir="."):
        r1_files = sorted(Path(fq_dir).glob("*R1.fastq.gz"))
        r2_files = sorted(Path(fq_dir).glob("*R2.fastq.gz"))
        for r1, r2 in zip(r1_files, r2_files):
            cmd = f"fastp -i {r1} -I {r2} -o {r1}.out1.gz -O {r2}.out2.gz --thread {self.threads}"
            self.run_cmd(cmd, f"Fastp QC: {r1.name}")

    def bowtie2_build(self, fasta, db_name):
        cmd = f"bowtie2-build {fasta} {db_name}"
        self.run_cmd(cmd, "Bowtie2 DB Build")

    def bowtie2_align(self, db_path, fq_dir="."):
        r1_files = sorted(Path(fq_dir).glob("*R1.fastq.gz"))
        r2_files = sorted(Path(fq_dir).glob("*R2.fastq.gz"))
        for r1, r2 in zip(r1_files, r2_files):
            cmd = f"bowtie2 -x {db_path} -1 {r1} -2 {r2} -p {self.threads} -q --al-conc-gz {r1}.bow.gz > tmp.txt"
            self.run_cmd(cmd, f"Bowtie2 Align: {r1.name}")

    def seqkit_stats(self, pattern="*bow*", output="StatsTable.txt"):
        cmd = f"seqkit stats {pattern} -T -j {self.threads} -o {output}"
        self.run_cmd(cmd, "SeqKit Stats")

    def kma_index(self, fasta, db_name):
        cmd = f"/storage/bioinfo/kma_2022/kma/kma_index -i {fasta} -o {db_name}"
        self.run_cmd(cmd, "KMA Index build")

    def kma_mapping(self, bow_dir, db_path, suffix="kma"):
        bow1_files = sorted(Path(bow_dir).glob("*.bow.1.gz"))
        bow2_files = sorted(Path(bow_dir).glob("*.bow.2.gz"))
        for f1, f2 in zip(bow1_files, bow2_files):
            out_prefix = f"{f1.stem}.{suffix}"
            cmd = f"/storage/bioinfo/kma_2022/kma/kma -ipe {f1} {f2} -o {out_prefix} -t_db {db_path} -t {self.threads}"
            self.run_cmd(cmd, f"KMA mapping: {f1.name}")

    def seqkit_normalization(self, bow_dir, sample_size=13020097):
        bow_files = sorted(Path(bow_dir).glob("*bow*"))
        for bf in bow_files:
            cmd = f"seqkit sample {bf} -n {sample_size} -j {self.threads} -o {bf}.meannorm"
            self.run_cmd(cmd, f"SeqKit normalization: {bf.name}")

if __name__ == "__main__":
    pipeline = GenomicPipeline(threads=24)
    # Adjust these paths and flags as needed for your workflow:
    pipeline.fastp_qc("./fastq")
    pipeline.bowtie2_build("all_db.fasta", "BOWTIE2_DB_2024")
    pipeline.bowtie2_align("/storage/bioinfo/ResCap_analysis/dbs/bowtie_db/full_platform", "./fastq")
    pipeline.seqkit_stats("*bow*", "StatsTable.txt")
    pipeline.kma_index("all_db_2022.fsa", "kma_db_2022")
    pipeline.kma_mapping("./bowtie_out", "/storage/bioinfo/DB/KMA_2024/originales/library_kma_2024", "kma")
    pipeline.seqkit_normalization("./bowtie_out", 13020097)
