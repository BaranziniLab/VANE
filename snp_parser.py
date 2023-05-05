import sqlite3
import pandas as pd 
import myvariant

def get_position_snp(rsid, verbose):
    """Query MyVariant to get the HG38 position at each rsID. """
    mv = myvariant.MyVariantInfo()
    query = mv.query('dbsnp.rsid:{}'.format(rsid),  assembly='hg38', fields='dbsnp')

    if query['total'] > 0:

        chromosome = query['hits'][0]['dbsnp']['chrom']
        alt = query['hits'][0]['dbsnp']['alt']

 
        ref = query['hits'][0]['dbsnp']['ref']
        

        try:
            hg38 = str(query['hits'][0]['dbsnp']['hg38']['start'])
        except:
            if verbose:
                print('{rsid} not found position'.format(rsid=rsid))
            hg38 = -1

        df_snp = pd.DataFrame(zip([chromosome], [hg38], [ref], [alt], [rsid]), columns = ['chromosome', 'position', 'ref', 'alt', 'rsid'])
    else:

        if verbose:
            print("Variant {variant} not found in MyVariant.info".format(variant=rsid))
        df_snp = pd.DataFrame(columns=['chromosome', 'position', 'ref', 'alt', 'rsid'])

    return df_snp

def create_df_for_locus(list_rs_ids, tagging_snp, verbose):
    """Accept List of rsids, and outputs a clean DF with positions at hg38, ref and alt allele.
    Tagging snp: is a rsid string"""
    list_of_df = []

    for variant in list_rs_ids:
        
        df_rsid = get_position_snp(variant, verbose)
        list_of_df.append(df_rsid)
    
    if len(list_of_df) > 0:

        df_result = pd.concat(list_of_df)
        df_result.dropna(inplace=True)
        df_result = df_result.loc[df_result['position'] != -1]
        df_result['position'] = df_result['position'].astype(int)
        df_result['tagging_snp'] = tagging_snp
        variant_number = len(df_result)

    else:
        df_result = pd.DataFrame()
    if verbose:
        print('There are total of {numb} mapped variants in LD with the tagging SNP'.format(numb = variant_number))
    return df_result


def query_snp_opentargets(df):

    """df: DF from create_df_for_locus()
    returns a dataframe with the SNPs in opentargets and the probability of affecting genes"""

    sql_dir = "base_files/clean_snp2gene.db"
    cnx = sqlite3.connect(sql_dir)

    position_list = tuple(df['position'].to_list())
    chromosome_list = tuple(df['chromosome'].to_list())
    SQL_QUERY = """SELECT * FROM snp2gene WHERE position IN {pos} AND chr_id in {chrom}""".format(pos = position_list, chrom = chromosome_list)
    df_results = pd.read_sql_query(SQL_QUERY, cnx)
    df_results.rename(columns={'chr_id':'chromosome'}, inplace=True)
    return df_results 

def merge_and_clean(snp_2_gene, df_with_position):
    """snp2gene: DF returned from query_snp_opentargets()
    df_with_position: DF returned from create_df_for_locus()
    returns a dataframe with each row containing the SNP and the highest probability per each gene"""
    snp_2_gene = snp_2_gene.merge(df_with_position, on=['chromosome', 'position'])
    snp_2_gene = snp_2_gene[['rsid', 'gene_id', 'position', 'chromosome', 'ref', 'alt', 'overall_score', 'tagging_snp']]
    idx = snp_2_gene.groupby('rsid')['overall_score'].transform(max) == snp_2_gene['overall_score']
    snp_2_gene = snp_2_gene[idx]

    gene_group = snp_2_gene.groupby('gene_id').count().reset_index()
    gene_group = gene_group[['gene_id', 'rsid']]
    return snp_2_gene, gene_group