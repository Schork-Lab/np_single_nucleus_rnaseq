#Python Libraries
import os

#Python Packages
import pandas as pd
from numpy import sum

def load_mapped_data(idx_dir, sample_map):
    '''
    Loads in the reads mapped data from a directory containing
    sample specific idxstats file generated using samtools idxstats.

    Returns a pandas Dataframe with indices as sample names based on filename,
    and three columns: ERCC, GENOME, UNMAPPED detailing the number of reads in each category

    :param idx_dir: path to the directory containing .idxstats files
    :type idx_dir: str
    :param sample_map: Sample map for changing names for figure creation
                        and readability
    :type sample_map: dict with key,value as sample_id, sample_name
    '''
    mapped_stats = {}
    for filename in os.listdir(idx_dir):
        if filename.endswith('.idxstats'):
            sample = filename.split('.')[0]
            filename = os.path.join(idx_dir, filename)
            df = pd.read_table(filename,
                                sep='\t',
                                names=['contig', 'length',
                                      'reads_mapped', 'reads_unmapped'])
            is_ercc = df.contig.str.contains('ERCC')
            ercc_mapped = df[is_ercc]['reads_mapped'].sum()
            other_mapped = df[~is_ercc]['reads_mapped'].sum()
            unmapped = df['reads_unmapped'].sum()
            mapped_stats[sample] = (ercc_mapped, other_mapped, unmapped)
    mapped_df = pd.DataFrame(mapped_stats).T
    mapped_df.columns = ['ERCC', 'GENOME', 'UNMAPPED']
    mapped_df = mapped_df.ix[sample_map]
    mapped_df.index = [sample_map[idx] for idx in mapped_df.index]
    return mapped_df

def load_tpms(rsem_dir, sample_map):
    '''

    :param rsem_dir: Path to the directory with the rsem genes.results files
    :type rsem_dir: str
    :param sample_map: Sample map for changing names for figure creation
                        and readability
    :type sample_map: dict with key,value as sample_id, sample_name
    '''

    def read_sample(filename, sample_name=None):
        '''
        Reads in RSEM generated .genes.results files, extracts
        TPM column and returns TPPM column only dataframe renamed to sample_name

        :param filename: path to the .genes.results file
        :type filename: str
        :param sample_name: name of sample for
        :type sample_name:
        '''

        if not sample_name:
            sample_name = os.path.basename(filename).split('.')[0]
        df = pd.read_table(filename, index_col=0)
        df[sample_name] = df['TPM']
        return df[[sample_name]]

    sample_dfs = []
    for fn in os.listdir(rsem_dir):
        sample_id = fn.split('.')[0]
        if fn.endswith('.genes.results') and sample_id in sample_map:
            sample_name = sample_map[sample_id]
            sample_df = read_sample(os.path.join(rsem_dir, fn),
                                    sample_name)
            sample_dfs.append(sample_df)

    tpm_df = pd.concat(sample_dfs, axis=1)

    return tpm_df




def filter_df(tpm_df,
              genes_only=True,
              expressed_in_multiple=False,
              neuronal=True, non_neuronal=True, total=True):
    '''

    :param tpm_df:
    :type tpm_df:
    :param genes_only:
    :type genes_only:
    :param neuronal:
    :type neuronal:
    :param non_neuronal:
    :type non_neuronal:
    :param total:
    :type total:
    '''

    if genes_only:
        tpm_df = tpm_df.ix[[index for index in tpm_df.index
                            if index.find('ENSG') != -1]]

    if expressed_in_multiple:
        tpm_df = tpm_df[sum(tpm_df > 0, axis=1) > 1]

    columns = []
    if neuronal:
        columns += [column for column in tpm_df.columns
                    if column.find('Neuronal nucleus') != -1]
    if non_neuronal:
        columns += [column for column in tpm_df.columns
                    if column.find('Non-neuronal nucleus') != -1]
    if total:
        columns += [column for column in tpm_df.columns
                    if column.find('Total RNA') != -1]

    return tpm_df[columns].sort(axis=1)

def get_low_mid_high_genes(control_tpm, cutoffs=(0.3333333, 0.6666667)):
    '''

    :param control_tpm:
    :type control_tpm:
    :param cutoffs:
    :type cutoffs:
    '''
    control_expressed = control_tpm[control_tpm > 0]
    first_quantile = control_expressed.quantile(cutoffs[0])
    second_quantile = control_expressed.quantile(cutoffs[1])
    low_cutoff = control_expressed <= first_quantile
    mid_cutoff = ((control_expressed > first_quantile) &
                  (control_expressed < second_quantile))
    high_cutoff = (control_expressed >= second_quantile)
    low_expressed = control_expressed[low_cutoff]
    mid_expressed = control_expressed[ mid_cutoff]
    high_expressed = control_expressed[high_cutoff]
    return low_expressed, mid_expressed, high_expressed

def label_expression(tpm_df, low, mid, high):
    def label_type(gene_id):
        if gene_id in low:
            return "Low"
        elif gene_id in mid:
            return "Mid"
        elif gene_id in high:
            return "High"
        else:
            return "Novel"

    tpm_df['Type'] = tpm_df.index.map(label_type)
    return tpm_df

def calculate_relative_expression(tpm_df, low, mid, high):

    tpm_df = label_expression(tpm_df, low, mid, high)
    value_counts = {}
    for sample_name in tpm_df.columns:
        sample_df = tpm_df[tpm_df[sample_name] > 0]
        value_counts[sample_name] = sample_df['Type'].value_counts()
    return pd.DataFrame(value_counts).T.fillna(0)

def load_gene_body_coverage(coverage_file, sample_map):
    coverage_df = pd.read_table(coverage_file, index_col=0).T
    sample_name = lambda x: x.split('.')[0][1:]
    coverage_df.columns = coverage_df.columns.map(sample_name)
    coverage_df = coverage_df[[column for column in coverage_df.columns
                                        if column in sample_map]]
    coverage_df.columns = coverage_df.columns.map(lambda x: sample_map[x])
    coverage_range = coverage_df.max()-coverage_df.min()
    normalized_df = (coverage_df - coverage_df.min())/coverage_range
    return normalized_df, coverage_df

def load_bedtools_coverage(filename, min_reads=1, min_length=100):
    columns = ['Chromosome', 'Start', 'End', 'GTF-INFO', 'Reads', 'Basesgt0',
                'Length', 'Fractiongt0']
    coverage_df = pd.read_table(filename, names=columns)
    coverage_df = coverage_df[(coverage_df.Reads >= min_reads) &
                              (coverage_df.Length >= min_length)]
    return coverage_df
