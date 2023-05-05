import pandas as pd
import requests
import re
from bs4 import BeautifulSoup



def regulome_parse_snp(position, chromosome):

    """Integer Position at the Build CH38
    Chromosome is just the chromosome number or position"""

    position_b = position
    position_a = position - 1 
    r = requests.get("https://beta.regulomedb.org/regulome-search?regions=chr{chr}%3A{position_a}-{position_b}&genome=GRCh38/thumbnail=chromatin".format(position_a = position_a,
    position_b = position_b, chr = chromosome))
    soup = BeautifulSoup(r.text, 'html.parser')
    table = soup.find('table', class_='table-sortable')
    df = pd.read_html(str(table))[0]

    return df

def score_assign(df_results):

    """Takes the regulome df and return a scored df"""

    df_results.dropna(inplace=True)
    df_results.loc[df_results['Chromatin state'] == 'Quiescent/Low', 'score'] = 0
    df_results.loc[df_results['Chromatin state'] == 'Weak transcription', 'score'] = 0
    df_results.loc[df_results['Chromatin state'] == 'Strong transcription', 'score'] = +1
    df_results.loc[df_results['Chromatin state'].str.contains('Repressed'), 'score'] = -1
    df_results.loc[df_results['Chromatin state'].str.contains('repressed'), 'score'] = -1
    df_results.loc[df_results['Chromatin state'].str.contains('enhancer'), 'score'] = +1
    df_results.loc[df_results['Chromatin state'].str.contains('Enhancer'), 'score'] = +1
    df_results.loc[df_results['Chromatin state'].str.contains('TSS'), 'score'] = +1
    df_results.loc[df_results['Chromatin state'].str.contains('Heterochromatin'), 'score'] = -1
    df_results.loc[df_results['Chromatin state'].str.contains('ZNF'), 'score'] = +0

    return df_results


def organ_df(organs, regulome_results):
    """Regolume results: Take a score regulome DF results (from score_assign())
    organs: Take a list of organs"""
    df_results = regulome_results
    list_of_organ = []
    for organ in organs:

        organ_results = df_results.loc[df_results['Organ'].str.contains(organ)]
        organ_results['organ_specific'] = organ
        list_of_organ.append(organ_results)
    
    df_organ = pd.concat(list_of_organ)

    df_scores = pd.DataFrame()

    df_scores['count'] = df_organ.groupby('organ_specific').count()['Chromatin state']
    df_scores['score'] = df_organ.groupby('organ_specific').sum()['score']
    df_scores['normalized score'] = df_scores['score'] / df_scores['count']


    return df_scores, df_organ

def cell_df(cells, regulome_results):
    """Regolume results: Take a score regulome DF results (from score_assign())
    cells: Take a list of cells"""
    df_results = regulome_results
    list_of_cell = []

    for cell in cells:

        cell_results = df_results.loc[df_results['Biosample'].str.contains(cell)]
        cell_results['cell_specific'] = cell
        list_of_cell.append(cell_results)
    
    df_cell = pd.concat(list_of_cell)

    df_scores = pd.DataFrame()

    df_scores['count'] = df_cell.groupby('cell_specific').count()['Chromatin state']
    df_scores['score'] = df_cell.groupby('cell_specific').sum()['score']
    df_scores['normalized score'] = df_scores['score'] / df_scores['count']

    return df_scores, df_cell