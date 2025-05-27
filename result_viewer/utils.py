import csv
from typing import List, Dict, Any
import xml.etree.ElementTree as ET

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

def interproscan_xml_to_dict(xml_file) -> List[Dict[str, Any]]:
    """
    Parse InterProScan XML output into a dictionary format
    
    Args:
        xml_file: File object containing InterProScan XML output
        
    Returns:
        List of dictionaries containing parsed matches
    """
    tree = ET.parse(xml_file)
    root = tree.getroot()
    
    # Define namespace
    ns = {'ipro': 'http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5'}
    
    results = []
    signature_results = {}
    
    # Find all matches
    for match in root.findall('.//ipro:match', ns):
        signature = match.find('ipro:signature', ns)
        if signature is None:
            continue
            
        signature_acc = signature.get('ac')
        
        if signature_acc not in signature_results:
            signature_results[signature_acc] = {
                'database': match.get('database'),
                'accession': signature_acc,
                'name': signature.get('name', ''),
                'description': signature.get('desc', ''),
                'interpro_acc': '',  # Will be updated if available
                'go_terms': [],
                'pathways': [],
                'locations': []
            }
            
            # Get InterPro reference if available
            entry_ref = signature.find('ipro:entry', ns)
            if entry_ref is not None:
                signature_results[signature_acc]['interpro_acc'] = entry_ref.get('ac', '')
                signature_results[signature_acc]['description'] = entry_ref.get('desc', signature_results[signature_acc]['description'])
                
                # Get GO terms
                for go_term in entry_ref.findall('.//ipro:go-xref', ns):
                    signature_results[signature_acc]['go_terms'].append(go_term.get('id'))
                    
                # Get pathway annotations
                for pathway in entry_ref.findall('.//ipro:pathway-xref', ns):
                    signature_results[signature_acc]['pathways'].append(pathway.get('id'))
        
        # Add location information
        for location in match.findall('ipro:locations/ipro:location', ns):
            loc_data = {
                'start': int(location.get('start')),
                'end': int(location.get('end')),
                'score': location.get('score', ''),
                'evalue': location.get('evalue', ''),
                'status': 'T'  # Default status for matches in XML
            }
            signature_results[signature_acc]['locations'].append(loc_data)
    
    # Convert dictionary to list
    results = list(signature_results.values())
    
    # Sort results by start position of first location
    results.sort(key=lambda x: x['locations'][0]['start'] if x['locations'] else 0)
    
    return results 