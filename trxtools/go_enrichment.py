author="Jan Mikołajczyk"

import pandas as pd
import requests
import warnings

go_mappings = {
    'biological_process': "GO:0008150",
    'molecular_function': "GO:0003674"
}
available_test_types = ["FISHER", "BINOMIAL"]
available_corrections = ["FDR", "BONFERRONI", "NONE"]

def get_enrichment(query_genes, organism, go_dataset='biological_process', use_reference_set=False, ref_genes=None, ref_organism=None, test_type="FISHER", correction="FDR"):
    '''Run a GO term enrichment test using PANTHER API

    :param query_genes: List of sequence identifiers of queried genes (e.g. transcript ids, gene ids)
    :type query_genes: list
    :param organism: Taxid of query species (e.g. "9606" for H. sapiens)
    :type organism: str
    :param go_dataset: Which annotation dataset to query, "biological_process" or "molecular_function"
    :type go_dataset: str
    :param use_reference_set: Use a custom set of rerence (background) genes? Default False. If True, ref_genes and ref_species need to be specifed.
    :type use_reference_set: bool
    :param ref_genes: *optional* list of reference genes. Specifying None (default) will use the whole genome of species specified in organism. When passing a list, ref_organism taxid must also be provided.
    :type ref_genes: list
    :param ref_species: Taxid of reference species, required when ref_genes is not None
    :type ref_species: str
    :param test_type: Which tatistical test to use. Available: "FISHER" (default), "BINOMIAL"
    :type test_type: str
    :param correction: Which multiple testing correction method to use. Available: "FDR" (default), "BONFERRONI", "NONE"
    :type correction: str
    :returns: Unfiltered DataFrame of results.
    :rtype: pandas.DataFrame:
    '''
    ## return empty df if empty list provided as input
    if len(query_genes) == 0:
        dummy_df = pd.DataFrame({
            "number_in_list": [], "fold_enrichment": [], "fdr": [],
            "expected": [], "number_in_reference": [], "pValue": [],
            "plus_minus": [], "term.id": [], "term.label": []
        })
        warnings.warn("Empty list provided, returning empty DataFrame.")
        return dummy_df
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
    response = requests.post('https://pantherdb.org/services/oai/pantherdb/enrich/overrep', data=query)
    ## extract data to a dataframe
    try:
        out_df = pd.json_normalize(response.json()['results']['result'])
    except KeyError:
        print("Response from server:")
        print(response.json())
        raise ValueError("Improper request data provided, check output for details")
    return out_df
