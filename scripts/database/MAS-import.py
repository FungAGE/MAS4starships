import mysql.connector
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import os
import argparse

import gzip
import re
from tqdm import tqdm
# import csv

parser = argparse.ArgumentParser(description="MAS SQL import")
parser.add_argument("-f", "--features", dest="feature_file", help="Path to features file", required=True)
parser.add_argument("-a", "--prot", dest="aa_file", help="Path to proteome", required=True)

fasta_group = parser.add_mutually_exclusive_group(required=True)
fasta_group.add_argument("-n", "--fasta", dest="fa_file", help="Path to fasta file")
fasta_group.add_argument("-nd", "--fasta_dir", dest="fa_dir", help="Path to directory containing fasta files")

gff_group = parser.add_mutually_exclusive_group(required=True)
gff_group.add_argument("-g", "--gff", dest="gff_file", help="Path to the input GFF file")
gff_group.add_argument("-gd", "--gff_dir", dest="gff_dir", help="Path to directory containing GFF files")

# Parse the arguments
args = parser.parse_args()

if args.fa_file:
    fa_files = [args.fa_file]
elif args.fa_dir:
    fa_files = [os.path.join(args.fa_dir, filename) for filename in os.listdir(args.fa_dir)]

if args.gff_file:
    gff_files = [args.gff_file]
elif args.gff_dir:
    gff_files = [os.path.join(args.gff_dir, filename) for filename in os.listdir(args.gff_dir)]

# python bin/MAS-import.py -n ../Starships/ships/fna/blastdb/concatenated.fa -d /home/adrian/Downloads/splitGff -g /home/adrian/Downloads/splitGff/altals1_s00058.gff -a /home/adrian/Downloads/test.faa.gz -f ../Starships/metadata/ships/starfish/output/mycodb.final.starships.feat 
# python bin/MAS-import.py -n ../Starships/ships/fna/blastdb/concatenated.fa -g /home/adrian/Downloads/splitGff/altals1_s00058.gff -a /home/adrian/Downloads/test.faa.gz -f ../Starships/metadata/ships/starfish/output/mycodb.final.starships.feat 

star_taxa = pd.read_csv("/home/adrian/Systematics/Starship_Database/Starships/MTDB/mycodb2899.ome2species.txt", sep="\t",names=['Code','Species'])
man_taxa = pd.read_csv("/home/adrian/Systematics/Starship_Database/Starships/metadata/ships/manual/Starships.fulltaxa.csv", sep=",")[['Code', 'Species']]
taxa = pd.concat([star_taxa,man_taxa], ignore_index=True)
taxa.set_index('Code', inplace=True)
taxa_dict = taxa.to_dict(orient='index')

# Establish a connection to the MySQL database
connection = mysql.connector.connect(
    host="0.0.0.0",
    port="3307",
    user="root",
    password="changeme",
    database="mas"
)

# Create a cursor object to execute SQL queries
cursor = connection.cursor()

def is_compressed(filename):
    compressed_extensions = ['.gz', '.zip']
    _, file_extension = os.path.splitext(filename)
    if file_extension in compressed_extensions:
        mode = 'rb'
    else:
        mode = 'r'
    return mode

def parse_features(feature_file):
    print(f"parsing {feature_file}...")
    df = pd.read_csv(feature_file, sep="\t")
    df.set_index('starshipID', inplace=True)
    feature_dict = df.to_dict(orient='dict')
    return feature_dict

# Function to find key name for a specific value
def find_key_for_value(dictionary, value):
    if dictionary is None or value is None:
        return None
    for key, values in dictionary.items():
        if values is not None and value in values:
            return key
    return None

def parse_gff(gff_files,feature_dict):
    pattern = r"Alias=([^;]+)"
    gff_dict = defaultdict(list)
    for gff_file in gff_files:
        print(f"parsing {gff_file}...")
        # mode = is_compressed(gff_file)
        with gzip.open(gff_file, "rt") as gff:
            for line in gff:
                parts = line.strip().split('\t')
                chr = parts[0]
                ship_ids = [find_key_for_value(feature_dict['#contigID'], chr)]
                if len(ship_ids) != 0:
                    for ship_id in ship_ids:
                        try:
                            if int(feature_dict['elementBegin'][ship_id]) <= int(parts[3]) <= int(feature_dict['elementEnd'][ship_id]):
                                # gene_id = parts[8].split("=")[1]
                                gene_id = re.search(pattern, parts[8]).group(1)
                                
                                ship_start = int(feature_dict['elementBegin'][ship_id])
                                ship_end = int(feature_dict['elementEnd'][ship_id])
                                coord1 = int(parts[3])
                                coord2 = int(parts[4])
                                start = abs(ship_start - coord1) + 1
                                end = abs(coord2 - coord1) + start
                                strand = parts[6]
                                out = {'chr': chr, 'start': start, 'stop': end, 'type': "gene", 'strand': strand, 'ship_id': ship_id, 'gene_id': gene_id}
                                gff_dict[parts[8]].append(out) 
                                # break
                        except KeyError:
                            pass
    return gff_dict

def parse_fa(fa_files):
    geno_id = 0
    fa_dict = defaultdict(list)
    for fa_file in fa_files:
        print(f"parsing {fa_file}...")
        with open(fa_file, "r") as fa:
            fa_records = SeqIO.parse(fa, "fasta")
            for fa_record in fa_records:
                fa_id = str(fa_record.id).replace("|-", "").replace("|+", "").replace("|", "")

                if "_s" in fa_id or any(char.isupper() for char in fa_id): # FIXME: hacky solution
                    key = fa_id.split("_")[0]
                else:
                    key = fa_id
                try:
                    org = str(taxa_dict[key]['Species'])
                except KeyError:
                    org = "Unknown"
                fa_seq = str(fa_record.seq)
                geno_id += 1
                out = {'geno_id': geno_id, 'fa_id': fa_id, 'seq': fa_seq, 'org': org}            
                fa_dict[fa_id].append(out)
    return fa_dict

def retrieve_match(table,col,key,query):
    sql_check = f"SELECT {col} FROM {table} WHERE {key} = %s"
    cursor.execute(sql_check, (query))

def sql_check_existing(table,key,query):
    check = f"SELECT COUNT(*) FROM {table} WHERE {key} = %s"
    cursor.execute(check, (query))
    count = cursor.fetchone()[0]
    return count

def insert_genome(table, geno_id, genome_name, genome_sequence, organism):
    # sql_check = sql_check_existing(table, key, genome_name)
    # if sql_check == 0:
    sql = f"INSERT INTO {table} (id, genome_name, genome_sequence, organism) VALUES (%s, %s, %s, %s)"
    cursor.execute(sql, (geno_id, genome_name, genome_sequence, organism))
    connection.commit()

def insert_feature(table, start, stop, type, strand, ann_id, geno_id):
    # sql_check = sql_check_existing(table, key, ann_id)
    # if sql_check == 0:
    sql = f"INSERT INTO {table} (start,stop,type,strand,annotation_id,genome_id) VALUES (%s, %s, %s, %s, %s, %s)"
    cursor.execute(sql, (start, stop, type, strand, ann_id, geno_id))
    connection.commit()

def insert_annotation(table,id,gene_id, sequence):
    # sql_check = sql_check_existing(table, key, gene_id)
    # if sql_check == 0:
    sql = f"INSERT INTO {table} (id, annotation, sequence, flag, assigned_to_id) VALUES (%s, %s, %s, %s, %s)"
    cursor.execute(sql, (id,gene_id, sequence, "0", "1"))
    connection.commit()

def save_tsv(output_file, *strs):
    with open(output_file, 'w', newline='') as tsvfile:
        out_str = '\t'.join(map(str, strs))
        tsvfile.write(out_str + '\n')

def mas_import(aa_file,gff_files,fa_files,feature_file):
    feature_dict = parse_features(feature_file)
    fa_records = parse_fa(fa_files)
    gff_records = parse_gff(gff_files,feature_dict)
    ann_id = 0
    # mode = is_compressed(aa_file)
    with gzip.open(aa_file, "rt") as aa:
        # for aa_record in tqdm(SeqIO.parse(aa, 'fasta'), desc="Iterating through ", unit=" fasta records"):
        for aa_record in SeqIO.parse(aa, 'fasta'):
            aa_id = str(aa_record.id)
            if aa_id in gff_records:
                gff_sub = gff_records[find_key_for_value(gff_records['gene_id'], aa_id)][0]
                ship_id = gff_sub["ship_id"]
                seq = str(aa_record.seq)
                if ship_id is not None:
                    if ship_id in fa_records:
                        ann_id += 1
                        fa_sub = fa_records[ship_id][0]
                        # save_tsv("genomes.tsv",fa_sub["geno_id"],fa_sub["fa_id"],fa_sub["seq"], fa_sub["org"])
                        insert_genome("genome_genome",fa_sub["geno_id"],fa_sub["fa_id"], fa_sub["seq"], fa_sub["org"])
                        print(f"{ship_id} sequence added")
                        # save_tsv("annotations.tsv",aa_id,seq,"0","0")
                        insert_annotation("genome_annotation",ann_id,aa_id,seq)
                        print(f"\t{ann_id} annotation added")
                        # save_tsv("features.tsv",gff_sub["start"],gff_sub["stop"],gff_sub["type"],gff_sub["strand"],ann_id,fa_sub["geno_id"])
                        insert_feature("genome_feature",gff_sub["start"],gff_sub["stop"],gff_sub["type"],gff_sub["strand"],ann_id,fa_sub["geno_id"])
                        print(f"\t\t{ann_id} feature added")
                    # else:
                        # print(f"{ship_id} fasta not found")
                # else:
                    # print(f"{os.basename(ship_id)} not found")
            # else:
                # print(f"{aa_id} not found in GFF")

mas_import(args.aa_file,gff_files,fa_files,args.feature_file)

# Close cursor and connection
cursor.close()
connection.close()
