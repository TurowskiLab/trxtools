import pandas as pd
import requests

go_mappings = {
    'biological_process': "GO:0008150",
    'molecular_function': "GO:0003674"
}
available_test_types = ["FISHER", "BINOMIAL"]
available_corrections = ["FDR", "BONFERRONI", "NONE"]

def get_enrichment(query_genes, organism, use_reference_set=False, ref_genes=None, ref_organism=None, go_dataset='biological_process', test_type="FISHER", correction="FDR"):
    '''
    Run a GO term enrichment test using PANTHER API
    :param query_genes: list() of sequence identifiers of queried genes (e.g. transcript ids, gene ids)
    :param organism: str() taxid of query species (e.g. "9606" for H. sapiens)
    :param use_reference_set: Use a custom set of rerence (background) genes? Default False. If True, ref_genes and ref_species need to be specifed.
    :param ref_genes: *optional* list() of reference genes. Specifying None (default) will use the whole genome of species specified in organism. When passing a list, ref_organism taxid must also be provided.
    :param ref_species: str() taxid of reference species, required when ref_genes is not None
    :param go_dataset: str() which annotation dataset to query, "biological_process" or "molecular_function"
    :param test_type: str() statistical test to use. Available: "FISHER" (default), "BINOMIAL"
    :param correction: str() multiple testing correction method. Available: "FDR" (default), "BONFERRONI", "NONE"
    :return: unfiltered DataFrame of results
    '''
    ## check if passed parameters are correct
    if go_dataset not in go_mappings.keys():
        raise Exception('Incorrect go_dataset value specified. Available values: "biological_process", "molecular_function"')
    if test_type not in available_test_types:
        raise Exception('Incorrect test_type value specified. Available values: "FISHER", "BINOMIAL"')
    if correction not in available_corrections:
        raise Exception('Incorrect correction value specified. Available values: "FDR", "BONFERRONI", "NONE"')
    ## format provided lists into strings and construct query
    seq_ids = ",".join(query_genes)
    if use_reference_set:
        if ref_genes is None:
            raise Exception("Reference gene list needs to be provided when use_reference_set is True")
        if ref_organism is None:
            raise Exception("Reference organism taxid needs to be provided when use_reference_set is True")
        ref_ids = ",".join(ref_genes)
        query = {
        "geneInputList": seq_ids,
        "organism": organism,
        "refInputList": ref_ids,
        "refOrganism": ref_organism,
        "annotDataSet": go_mappings[go_dataset],
        "enrichmentTestType": test_type,
        "correction": correction
    }
    else:
        query = {
            "geneInputList": seq_ids,
            "organism": organism,
            "annotDataSet": go_mappings[go_dataset],
            "enrichmentTestType": test_type,
            "correction": correction
        }
    response = requests.post('http://pantherdb.org/services/oai/pantherdb/enrich/overrep', data=query)
    ## extract data to a dataframe
    return pd.json_normalize(response.json()['results']['result'])
