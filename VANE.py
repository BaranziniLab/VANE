from re import S
import pandas as pd
from ld_parser import get_ld_snps, process_api_call
from snp_parser import get_position_snp, create_df_for_locus, merge_and_clean
import sqlite3
from regulome_db_parsing import regulome_parse_snp, score_assign, organ_df, cell_df
import seaborn as sns
import networkx as nx
from pyvis.network import Network
import os
import random as rd
import matplotlib.pyplot as plt
import numpy as np

class VANE:

    def __init__(self, list_of_variants = None, mapping_table_dir = None, sql_dir = None):
        """list of variants: List composed of tagging SNPs in rsid format."""
        #if sql_dir == None:
            #self.sql_dir = "/pool0/home/secret/snp2gene_clean/snp2_gene_sql/clean_snp2gene.db"
        
        if mapping_table_dir == None:
            mapping_table = pd.read_table('/pool0/home/secret/regulomedb/streamlit_app/data_tables/gene_protein_ensembl_map_table.tsv')

        self.mapping_table = mapping_table
        self.tagging_snps = list_of_variants
        self.snp_2_gene = None
        self.regulome_db = None


    def process_tag_variants(self, ld, population=None, verbose=False):
        """Populations can be extracted from here:  https://rest.ensembl.org/documentation/info/variation_populations, if None = "1000GENOMES:phase_3:CEU" 
        verbose: True/False"""

        tagging_snp_list = self.tagging_snps
        df_all_positons = []

        for rsid in tagging_snp_list:

            list_of_snps = get_ld_snps(rsid, ld, population)
            list_of_snps.append(rsid)
            df_with_position = create_df_for_locus(list_of_snps, rsid, verbose)
            df_all_positons.append(df_with_position)
        
        df_all_positons = pd.concat(df_all_positons)

        return df_all_positons

    def process_list_finemapped(self, list_of_locus_name, list_of_finemapped, verbose=False):

        """In case that you want to use fine-mapped SNPs for any given LOCUS, you can use this function.
        list_of_locus_name = List for how you want to call the locus. 
        list_of_finemapped = This is a list of lists. For each locus you want a list of the finemapped SNPs in rsid format.
        Both list_of_locus_name, and list_of_fiemapped SNPs need to be in order.
        
        Example: 
        You have three locus, locus 1, locus 2, and locus 3. You have three SNPs rs2, rs3, and rs4, finemapped to locus 1, 
        rs42, rs32, and rs11 finemapped to Locus 2, and rs111, rs321, and rs33, finemapped to locus 3.

        This how it should be formatted:

        list_of_locus_name = ['locus_1', 'locus_2', 'locus_3']
        list_of_finemapped = [['rs2', 'rs3', 'rs4'], ['rs42', 'rs32', 'rs11'], ['rs111', 'rs321', 'rs33']]"""

        df_all_positons = []

        for locus, fine_mapped_snps in zip(list_of_locus_name, list_of_finemapped):
          
            df_with_position = create_df_for_locus(fine_mapped_snps, locus, verbose=verbose)
            df_all_positons.append(df_with_position)
        
        df_all_positons = pd.concat(df_all_positons)

        return df_all_positons

        
    def query_snp_opentargets(self, df, sqlt3_dir):

        """df: DF from create_df_for_locus()
        returns a dataframe with the SNPs in opentargets and the probability of affecting genes"""
        sql_dir = sqlt3_dir
        #sql_dir = "/pool0/home/secret/snp2gene_clean/snp2_gene_sql/clean_snp2gene.db"
        
        cnx = sqlite3.connect(sql_dir)

        position_list = tuple(df['position'].to_list())
        chromosome_list = tuple(df['chromosome'].to_list())
        SQL_QUERY = """SELECT * FROM snp2gene WHERE position IN {pos} AND chr_id in {chrom}""".format(pos = position_list, chrom = chromosome_list)
        df_results = pd.read_sql_query(SQL_QUERY, cnx)
        df_results.rename(columns={'chr_id':'chromosome'}, inplace=True)

        return df_results 

    def merge_and_clean(self, snp_2_gene, df_with_position, threshold = 0):
        """snp2gene: DF returned from query_snp_opentargets()
        df_with_position: DF returned from create_df_for_locus()
        threshold: Amount of rsIDs that need to map to a gene to be used in the analysis.

        returns a dataframe with each row containing the SNP and the highest probability per each gene"""

        snp_2_gene_clean = snp_2_gene.copy()
        snp_2_gene_clean = snp_2_gene_clean.merge(df_with_position, on=['chromosome', 'position'])
        snp_2_gene_clean = snp_2_gene_clean[['rsid', 'gene_id', 'position', 'chromosome', 'ref', 'alt', 'overall_score', 'tagging_snp']]
        idx = snp_2_gene_clean.groupby('rsid')['overall_score'].transform(max) == snp_2_gene_clean['overall_score']
        snp_2_gene_clean = snp_2_gene_clean[idx]

        gene_group = snp_2_gene_clean.groupby('gene_id').count().reset_index()
        gene_group = gene_group[['gene_id', 'rsid']]

        genes_ = gene_group.loc[gene_group['rsid'] > threshold]['gene_id'].to_list()
        snp_2_gene_clean = snp_2_gene_clean.loc[snp_2_gene_clean['gene_id'].isin(genes_)]

        self.snp_2_gene_clean = snp_2_gene_clean
        
        return gene_group

        
    def regulomedb_scoring(self, snp_2_gene_all=None, save=True, verbose=False):
        """Take the SNP 2 GENE cleaned and merge dataframe, query regulomeDB, and score the chromatin state."""
        if snp_2_gene_all == None:
            
            snp_2_gene_all = self.snp_2_gene_clean

        total_number = len(snp_2_gene_all)
        list_of_results = []
        a = 0
        for index, row in snp_2_gene_all.iterrows():

            a += 1
            if verbose:
                print("{number} SNPs processed out of {total}".format(number = a, total = total_number))
            chromomosome = row['chromosome']
            position = row['position']
            gene = row['gene_id'] 
            df_results = regulome_parse_snp(position, chromomosome)
            df_results['gene'] = gene
            list_of_results.append(df_results)
        
        result_df_all = pd.concat(list_of_results)
        result_df_all = score_assign(result_df_all)

        if save:

            self.regulome_db = result_df_all
            return 

        else:

            return result_df_all


    def process_cells(self, cells, result_df_all=None, save=True):
        """Takes a cell list, and returns a processed files with the cells selected.
        It also saves in to the attribute df_for_network, a dataframe that can be used to generate networks"""


        if result_df_all == None:
            result_df_all = self.regulome_db

        gene_list = result_df_all['gene'].unique()
        results_cells_per_gene = []
        
        for gene in gene_list:

            gene_df = result_df_all.loc[result_df_all['gene'] == gene]
            cells_scored, df_cells = cell_df(cells, gene_df)
            cells_scored['gene'] = gene
            results_cells_per_gene.append(cells_scored)
        
        cell_results = pd.concat(results_cells_per_gene).reset_index()

        if save:
            self.cell_results = cell_results

        result_to_plot = self.cell_results.merge(self.mapping_table, on='gene')
        plot_df = result_to_plot[['normalized score', 'protein_gene', 'cell_specific']]

        plot_df.drop_duplicates(inplace=True)
        plot_df_pivoted = plot_df.pivot('protein_gene', 'cell_specific')
        self.df_for_network = plot_df_pivoted


        return cell_results

    def process_organs(self, organs, result_df_all=None, save=True):

        if result_df_all == None:
            result_df_all = self.regulome_db

        organ_lists = []
        gene_list = result_df_all['gene'].unique()

        for gene in gene_list:

            gene_df = result_df_all.loc[result_df_all['gene'] == gene]
            organ_scored, df_organ = organ_df(organs, gene_df)
            organ_scored['gene'] = gene
            organ_lists.append(organ_scored)

        organs_results = pd.concat(organ_lists).reset_index()

        if save:
            self.organ_results = organs_results

        return organs_results

    def plot_cells(self, scale_heatmap=False, standarize_gene = False):
        
        """Scale_heatmap: If True scale heatmap from -1 to 1."""
        result_to_plot = self.cell_results.merge(self.mapping_table, on='gene')

        plot_df = result_to_plot[['normalized score', 'protein_gene', 'cell_specific']]

        plot_df.drop_duplicates(inplace=True)
        plot_df_pivoted = plot_df.pivot('protein_gene', 'cell_specific')
        gene_standarized_cell = ((plot_df_pivoted.T - plot_df_pivoted.T.mean()) / plot_df_pivoted.T.std()).T
        
        self.cell_results_full = plot_df_pivoted
        
        if standarize_gene:

            standard_image_cell = sns.heatmap(gene_standarized_cell, cmap='vlag')
            self.standardized_gene = standard_image_cell
        else:

            if scale_heatmap:
                standard_image_cell = sns.heatmap(plot_df_pivoted, cmap='vlag', vmin= -1, vmax=1)
            else:
                standard_image_cell = sns.heatmap(plot_df_pivoted, cmap='vlag')
            
        return standard_image_cell

    def plot_organs(self, scale_heatmap=False, standarize_gene = False):
        
        """Scale_heatmap: If True scale heatmap from -1 to 1. If False, auto scale."""

        organs_results = self.organ_results.merge(self.mapping_table, on='gene')
        organs_results.drop_duplicates(inplace=True)
        organ_plot = organs_results[['organ_specific', 'normalized score', 'protein_gene']]
        organ_plot.drop_duplicates(subset=['protein_gene', 'organ_specific'], inplace=True)
        organ_plot = organ_plot.pivot('protein_gene', 'organ_specific')
        self.organ_df = organ_plot

        if standarize_gene:
            gene_standarized_organs = ((organ_plot.T - organ_plot.T.mean()) / organ_plot.T.std()).T
            standard_image_organ = sns.heatmap(gene_standarized_organs, cmap='vlag', vmin= -1, vmax=1)
        else:

            if scale_heatmap:
                standard_image_organ = sns.heatmap(organ_plot, cmap='vlag', vmin= -1, vmax=1)

            else:
                standard_image_organ = sns.heatmap(organ_plot, cmap='vlag')

        return standard_image_organ


class VANE_network():

    def __init__(self, network_folder, df_for_network, mapping_table = None, cell_dicts = None):
        """list of variants: List composed of tagging SNPs in rsid format.
        df_for_network is df from process_self, can be called with self.df_for_network() """
        
        self.df_for_network = df_for_network.reset_index()

        self.network_folder = network_folder

        self.cell_dicts =cell_dicts
        if cell_dicts == None:
            self.cell_dicts = {'B':'B-Cell', 'T':'T-Cell', 'mono':'Monocytes', 'astrocyte':'Astrocytes', 'oligodendrocyte':'Oligodendrocytes'}
        

        self.mapping_table = mapping_table
        if mapping_table == None:
            self.mapping_table = pd.read_table('base_files/gene_protein_ensembl_map_table.tsv')


    def create_cellular_network(self, cell, epigenome_threshold = 0, custom=False):
        cell_name = self.cell_dicts[cell]

        #load graph
        cell_graph = self.load_cell_network(cell_name)
        self.original_graph = cell_graph
        #Get subset 
        df_protein_expressed = self.df_for_network
        if custom:
            list_genes = df_protein_expressed.loc[df_protein_expressed[cell] >= epigenome_threshold]['cell_specific'].to_list()            
        else:
            list_genes = df_protein_expressed[('protein_gene','')].loc[df_protein_expressed['normalized score', cell] > epigenome_threshold].to_list()
        print(list_genes)

        if len(list_genes) == 0:

            print('List of Gene empty, try with a lower Epigenomic threshold')

        else:

            proteins_to_extract = self.mapping_table.loc[self.mapping_table['protein_gene'].isin(list_genes)]['protein'].to_list()
            cell_subgraph = cell_graph.subgraph(proteins_to_extract)
            self.graph = cell_subgraph

            return cell_subgraph

        
    def load_cell_network(self, cell_name):
        """cell name comes from the form
        cell dict is a cell from regulome db to cells from protein atlas
        """
        network_folder = self.network_folder
        cell_list = os.listdir(self.network_folder)
        cell_file_dir = [cell for cell in cell_list if cell_name.lower() in cell.lower()]
        cell_graph = nx.read_gml(self.network_folder + cell_file_dir[0])

        return cell_graph


    def plot_network(self, html_name):
        """once you create a cellular network with create_cellular_network(), you can plot the network using pyvis.
        This function returns a HTML file that can be opened any browser
        html_name = any name for the html"""
        nt = Network('900px', width='100%',
                directed = False)

        nt.from_nx(self.graph)
        nt.show_buttons()

        html_name = html_name + '.html'

        nt.show(html_name)

        return print('Network generate at {html_name}'.format(html_name = html_name))

    def calculate_p_value(self, n = 1000):
        """You need to run creat_cellular_network(), before running this function.
        It returns a set with the 5 and 95 CI, and one-tailed pval."""
        graph_node_sample = self.graph.number_of_nodes()
        graph_edge_original = self.graph.number_of_edges()

        list_numb_nodes = []
        list_numb_edges = []

        for x in range(n):

            sample = rd.sample(self.original_graph.nodes, graph_node_sample)
            graph_sample = self.original_graph.subgraph(sample)
            number_of_nodes = graph_sample.number_of_nodes()
            number_of_edges = graph_sample.number_of_edges()

            list_numb_nodes.append(number_of_nodes)
            list_numb_edges.append(number_of_edges)

        ci_5 = np.percentile(list_numb_edges, 5)
        ci_95 = np.percentile(list_numb_edges, 95)
        
        number_bigger = len([x for x in list_numb_edges if x > graph_edge_original])
        pval = (1+number_bigger)/(len(list_numb_edges)+1)

        return (ci_5, ci_95), pval

    def plot_bootstrap(self):
        
        graph_node_sample = self.graph.number_of_nodes()
        graph_edge_original = self.graph.number_of_edges()

        list_numb_nodes = []
        list_numb_edges = []

        for x in range(1000):

            sample = rd.sample(self.original_graph.nodes, graph_node_sample)
            graph_sample = self.original_graph.subgraph(sample)
            number_of_nodes = graph_sample.number_of_nodes()
            number_of_edges = graph_sample.number_of_edges()

            list_numb_nodes.append(number_of_nodes)
            list_numb_edges.append(number_of_edges)

        
        sns.set_theme()
        plt.axvline(x=graph_edge_original, ymin=0, ymax=0.95,color='r')
        figure = sns.histplot(list_numb_edges, binwidth=0.3)
        figure.set(xlabel='Edge Count', ylabel='Frequency in bootstrapping')

        return figure

        