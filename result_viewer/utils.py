import csv
from typing import List, Dict, Any

def interproscan_tsv_to_dict(tsv_file) -> List[Dict[str, Any]]:
    """
    Parse InterProScan TSV output into a dictionary format
    
    TSV Format columns:
    0. Protein accession
    1. Sequence MD5 digest
    2. Sequence length
    3. Analysis (database)
    4. Signature accession
    5. Signature description
    6. Start location
    7. Stop location
    8. Score
    9. Status
    10. Date
    11. InterPro accession
    12. InterPro description
    13. GO annotations
    14. Pathways annotations
    """
    results = []
    
    with open(tsv_file) as f:
        tsv_reader = csv.reader(f, delimiter='\t')
        
        # Group results by signature accession to combine locations
        signature_results = {}
        
        for row in tsv_reader:
            if len(row) < 11:  # Basic validation
                continue
                
            signature_acc = row[4]
            
            if signature_acc not in signature_results:
                signature_results[signature_acc] = {
                    'database': row[3],
                    'accession': signature_acc,
                    'name': row[5],  # Signature description
                    'description': row[12] if len(row) > 12 and row[12] else row[5],  # Use InterPro description if available
                    'interpro_acc': row[11] if len(row) > 11 else '',
                    'go_terms': row[13].split('|') if len(row) > 13 and row[13] else [],
                    'pathways': row[14].split('|') if len(row) > 14 and row[14] else [],
                    'locations': []
                }
            
            # Add location information
            signature_results[signature_acc]['locations'].append({
                'start': int(row[6]),
                'end': int(row[7]),
                'score': row[8] if row[8] != '-' else '',
                'evalue': '',  # TSV doesn't include e-value
                'status': row[9]
            })
    
    # Convert dictionary to list
    results = list(signature_results.values())
    
    # Sort results by start position of first location
    results.sort(key=lambda x: x['locations'][0]['start'] if x['locations'] else 0)
    
    return results 