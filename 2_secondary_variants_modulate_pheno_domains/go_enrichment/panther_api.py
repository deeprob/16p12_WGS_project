# http://pantherdb.org/services/tryItOut.jsp?url=%2Fservices%2Fapi%2Fpanther

import requests
import json
import pandas as pd

def panther(input_genes, reference_genes=None, organism='9606', annotation_dataset='GO:0008150', test_type='FISHER', correction_type='FDR'):
    '''
        Function to query Panther API. More info at Panther:
        http://pantherdb.org/services/tryItOut.jsp?url=%2Fservices%2Fapi%2Fpanther
        
        input_genes:
            list of genes to query
        
        reference_genes:
            (optional) list of reference genes
        
        organism: 
            all supported organisms listed at http://pantherdb.org/services/oai/pantherdb/supportedgenomes
            default is 9606 (human)
            
        annotation_dataset: 
            all supported annotations listed at http://pantherdb.org/services/oai/pantherdb/supportedannotdatasets
            default is GO:0008150 (Biological Process)
            ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP (Biological Process Slim)
            ANNOT_TYPE_ID_PANTHER_PATHWAY (Panther Pathway)
            ANNOT_TYPE_ID_REACTOME_PATHWAY (Reactome Pathway)
            
        test_type: FISHER, BINOMIAL
        
        correction_type: FDR, BONFERRONI, NONE
    '''
    
    # stringify input genes
    input_genes     = ','.join(input_genes)
    
    # stringify reference genes if provided
    if reference_genes != None:
        reference_genes = ','.join(reference_genes)
    
    #API details
    url = 'http://pantherdb.org/services/oai/pantherdb/enrich/overrep'

    query = {'geneInputList':input_genes,
                       'organism': organism,
                       'refInputList': reference_genes,
                       'refOrganism': organism,
                       'annotDataSet':annotation_dataset,
                       'enrichmentTestType': test_type,
                       'correction': correction_type
                      }
    
    # if no reference genes provided remove from query
    if reference_genes == None:
        del query['refInputList']
        del query['refOrganism']

    #Making http post request
    response = requests.post(url, data=query)

    res = response.json()
    res = res['results']
    res = pd.DataFrame(res['result'])
    res['term_id'] = res['term'].apply(lambda s: s['id'] if 'id' in s else 'N/A')
    res['term_label'] = res['term'].apply(lambda s: s['label'] if 'label' in s else 'N/A')
    res = res.drop('term', axis=1)
    # res = res[res.fdr < 0.05].copy()
    return res
