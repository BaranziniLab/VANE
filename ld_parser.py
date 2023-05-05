import requests, sys


def get_ld_snps(rsid, ld_treshold, population=None):
    """Populations can be extracted from here:  https://rest.ensembl.org/documentation/info/variation_populations"""
    if population == None:

        population = "1000GENOMES:phase_3:CEU"
        
    server = "https://rest.ensembl.org"
    ext = "/ld/human/{rsid}/{population}?window_size=500;r2={ld_treshold}".format(rsid=rsid, population = population, ld_treshold=ld_treshold)

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if r.status_code == 200:

        decoded = r.json()
        decoded = process_api_call(decoded)
    
        return decoded

    else:

        return []
        
    #if not r.ok:
       # r.raise_for_status()
       # sys.exit()
        
    #decoded = r.json()

   # decoded = process_api_call(decoded)
    
   # return decoded

def process_api_call(api_results):

    snp_list = []
    for x in api_results:

        snp_in_ld = x['variation2']
        snp_list.append(snp_in_ld)
    
    return snp_list

