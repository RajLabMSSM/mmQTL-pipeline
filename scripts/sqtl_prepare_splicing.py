import pickle
import pandas as pd
import numpy as np
import argparse
import subprocess
import sys
import os
import gzip
import contextlib
from datetime import datetime
import tempfile
import shutil
import glob
from sklearn.decomposition import PCA

## Leafcutter SQTLs - adapted for mmQTL

# jack humphrey - changes to original script by Francois Auget:

# renamed script from cluster_prepare_fastqtl.py to sqtl_prepare_splicing.py
# junc files are assumed not to be gzipped
# junc files are produced by regtools, not bam2junc
# clustering is done with leafcutter_cluster_regtools.py
# added --checkChrom to clustering call to ignore any contig chromosome names
# explicitly call python2 with the leafcutter scripts
# all leafcutter scripts and map_genes.R script are assumed to be in the same folder (--leafcutter_dir)
# line where glob looks for qq files was looking within --output_dir now looks in os.path.dirname(--prefix) - this is same as outFolder
# swapped out write_bed function for the one in the eQTL expression prepare_bed.py script which works with python3

# added sample_participant_lookup so that output leafcutter.bed.gz will have participant IDs, not sample IDs
# junc_files folder is removed at the end of the script to save space
# sorted.gz files from clustering, intermediate phenotype and bed files are also deleted.

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)


def gtf_to_bed(annotation_gtf, feature='gene'):
    """
    Parse genes from GTF, create placeholder DataFrame for BED output
    """
    chrom = []
    start = []
    end = []
    gene_id = []
    with open(annotation_gtf, 'r') as gtf:
        for row in gtf:
            row = row.strip().split('\t')
            if row[0][0]=='#' or row[2]!=feature: continue # skip header
            chrom.append(row[0])

            # TSS: gene start (0-based coordinates for BED)
            if row[6]=='+':
                start.append(np.int64(row[3])-1)
                end.append(np.int64(row[3]))
            elif row[6]=='-':
                start.append(np.int64(row[4])-1)  # last base of gene
                end.append(np.int64(row[4]))
            else:
                raise ValueError('Strand not specified.')

            gene_id.append(row[8].split(';',1)[0].split(' ')[1].replace('"',''))

    bed_df = pd.DataFrame(
        data={'chr':chrom, 'start':start, 'end':end, 'gene_id':gene_id},
        columns=['chr', 'start', 'end', 'gene_id'],
        index=gene_id)
    return bed_df


#def write_bed(bed_df, output_name):
#    """Write DataFrame to BED"""
#    bgzip = subprocess.Popen('bgzip -c > '+output_name,
#        stdin=subprocess.PIPE, shell=True)
#    bed_df.to_csv(bgzip.stdin, sep='\t', index=False)
#    stdout, stderr = bgzip.communicate()
#    subprocess.check_call('tabix -f '+output_name, shell=True)

def write_bed(bed_df, output_name):
    """
    Write DataFrame to BED format
    """
    #assert bed_df.columns[0]=='chr' and bed_df.columns[1]=='start' and bed_df.columns[2]=='end'
    # header must be commented in BED format
    header = bed_df.columns.values.copy()
    #header[0] = '#chr'
    bed_df.to_csv(output_name, sep='\t', index=False, header=True)
    subprocess.check_call('bgzip -f '+output_name, shell=True, executable='/bin/bash')
    subprocess.check_call('tabix -f '+output_name+'.gz', shell=True, executable='/bin/bash')



if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Run leafcutter clustering, prepare for FastQTL')
    parser.add_argument('junc_files_list', help='File with paths to ${sample_id}.junc files')
    parser.add_argument('exons', help='Exon definitions file, with columns: chr, start, end, strand, gene_id, gene_name')
    parser.add_argument('genes_gtf', help='Collapsed gene annotation in GTF format')
    parser.add_argument('prefix', help='Prefix for output files (sample set ID)')
    parser.add_argument('sample_participant_lookup', help='Lookup table linking samples to participants')
    parser.add_argument('--min_clu_reads', default='50', type=str, help='Minimum number of reads supporting each cluster')
    parser.add_argument('--min_clu_ratio', default='0.001', type=str, help='Minimum fraction of reads in a cluster that support a junction')
    parser.add_argument('--max_intron_len', default='500000', type=str, help='Maximum intron length')
    parser.add_argument('--num_pcs', default=5, type=int, help='Number of principal components to calculate')
    parser.add_argument('--coord_mode', default="TSS", type = str, help = "Whether to set the coordinates as the gene TSS (TSS) or the middle of the cluster (cluster_middle) ")
    parser.add_argument('--leafcutter_dir', default='/opt/leafcutter',
                        help="leafcutter directory, containing 'clustering' directory")
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
    args = parser.parse_args()

    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] leafcutter clustering')
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    with cd(args.output_dir):
        print('  * decompressing and renaming junc files')
        with open(args.junc_files_list) as f:
            all_junc_files = f.read().strip().split('\n')

        junc_dir = os.path.join(args.output_dir, 'junc_files')
        if not os.path.exists(junc_dir):
            os.mkdir(junc_dir)
        # read in sample participant lookup
        sample_participant_lookup_s = pd.read_csv(args.sample_participant_lookup, sep='\t', index_col=0, dtype=str, squeeze=True)

        # match sample_id to participant_id, use this as naming the junctions
        # do on final leafcutter bed
        sample_ids = []
        junc_files = []

        for f in all_junc_files:
            sample_id = os.path.split(f)[1].split('.')[0]
            # only if sample is in the sample key
            if sample_id in sample_participant_lookup_s.keys():
                sample_ids.append(sample_id)
                junc_files.append(f)
            else:
                print( "  * " + sample_id + " is not present in sample key - ignoring" )
        
        # rename samples using participant IDs from key
        sample_df = pd.DataFrame( columns = sample_ids )
        sample_df.rename(columns=sample_participant_lookup_s.to_dict(), inplace=True)
        sample_ids = sample_df.columns

        #print(junc_files)
        #print(sample_ids)

        # copy junctions using new IDs
        for s,j in zip(sample_ids, junc_files):
            out_path = os.path.join(junc_dir, s+'.junc') 
            shutil.copy2(j, out_path)
        
        #subprocess.check_call('gunzip -f '+os.path.join(junc_dir, '*.junc'), shell=True)
        junc_files = sorted([os.path.join(junc_dir, i+'.junc') for i in sample_ids])

        print('  * running leafcutter clustering')
        # generates ${prefix}_perind_numers.counts.gz and ${prefix}_perind.counts.gz
        with tempfile.NamedTemporaryFile(dir=args.output_dir) as tmp:
            with open(tmp.name, 'w') as f:
                f.write('\n'.join(junc_files)+'\n')
            subprocess.check_call(
                'python '+os.path.join(args.leafcutter_dir, 'leafcutter_cluster_regtools.py' \
                    +' --juncfiles '+tmp.name \
                    +' --outprefix '+args.prefix \
                    +' --checkchrom ' \
                    +' --minclureads '+args.min_clu_reads \
                    +' --mincluratio '+args.min_clu_ratio \
                    +' --maxintronlen '+args.max_intron_len), shell=True)

        print('  * compressing outputs')
        subprocess.check_call('gzip {}_pooled'.format(args.prefix), shell=True)
        subprocess.check_call('gzip {}_refined'.format(args.prefix), shell=True)

        print('  * filtering counts')
        counts_df = pd.read_csv(os.path.join(args.output_dir, args.prefix+'_perind.counts.gz'), sep='\s+').set_index('chrom')
        calculate_frac = lambda x: float(x[0])/float(x[1]) if x[1]>0 else 0
        frac_df = counts_df.applymap(lambda x: calculate_frac([int(i) for i in x.split('/')]))
        pct_zero = (frac_df==0).sum(axis=1) / frac_df.shape[1]  # for zero counts, frac is zero
        n_unique = frac_df.apply(lambda x: len(x.unique()), axis=1)
        zscore_df = ((frac_df.T-frac_df.mean(1)) / frac_df.std(1)).T

        # filter out introns with low counts or low complexity
        n = np.floor(frac_df.shape[1]*0.1)
        if n<10:
            n = 10
        mask = (pct_zero <= 0.5) & (n_unique >= n)
        # additional filter for low complexity
        ns = zscore_df.shape[1]
        mask2 = ((zscore_df.abs()<0.25).sum(1)>=ns-3) & ((zscore_df.abs()>6).sum(1)<=3)
        if np.any(mask & mask2):
            print('    ** dropping {} introns with low variation'.format(np.sum(mask & mask2)))
        mask = mask & ~mask2

        filtered_counts_df = counts_df.loc[mask]
        cluster_ids = np.unique(counts_df.index.map(lambda x: x.split(':')[-1]))
        filtered_cluster_ids = np.unique(filtered_counts_df.index.map(lambda x: x.split(':')[-1]))
        print('    ** dropping {} introns with counts in fewer than 50% of samples\n'
              '       {}/{} introns remain ({}/{} clusters)'.format(
                   counts_df.shape[0]-filtered_counts_df.shape[0], filtered_counts_df.shape[0], counts_df.shape[0], len(filtered_cluster_ids), len(cluster_ids)
            ))
        filtered_counts_file = os.path.join(args.output_dir, args.prefix+'_perind.counts.filtered.gz')
        with gzip.open(filtered_counts_file, 'wt') as f:
            filtered_counts_df.to_csv(f, sep=' ')
        # mapping clusters to genes should be done on the filtered junctions only
        print('  * mapping clusters to genes')
        subprocess.check_call(
            'Rscript' \
                +' '+os.path.abspath(os.path.join(args.leafcutter_dir, 'map_clusters_to_genes.R')) \
                +' '+os.path.join(args.output_dir, args.prefix+'_perind.counts.filtered.gz') \
                +' '+args.exons \
                +' '+args.prefix + '.leafcutter.clusters_to_genes.txt', shell=True)

        print('  * preparing phenotype table')
        subprocess.check_call(
            'python '+os.path.join(args.leafcutter_dir, 'prepare_phenotype_table.py') \
            +' '+filtered_counts_file \
            +' -p '+str(args.num_pcs), shell=True)

        print('  * concatenating BED files')
        bed_files = sorted(glob.glob(os.path.join(os.path.dirname(args.prefix), '*_perind.counts.filtered.gz.qqnorm_*')))
        print(bed_files)
        bed_df = []
        for f in bed_files:
            bed_df.append(pd.read_csv(f, sep='\t', dtype=str))
        bed_df = pd.concat(bed_df, axis=0, sort = False) # pandas gave me a warning - check whether this needs to be false or true
        print('    ** sorting')
        # leafcutter doesn't produce output for chrX --> numeric sort
        for c in ['#Chr', 'start', 'end']:
            bed_df[c] = bed_df[c].astype(int)
        bed_df = bed_df.sort_values(['#Chr', 'start', 'end'])
        print('    ** writing BED')
        bed_file = os.path.join(args.output_dir, args.prefix+'.perind.counts.filtered.qqnorm.bed.gz')
        bgzip = subprocess.Popen('bgzip -c > '+bed_file, stdin=subprocess.PIPE, shell=True, universal_newlines=True)
        stdout, stderr = bgzip.communicate(bed_df.to_csv(sep='\t', index=False))
        print('    ** indexing')
        subprocess.check_call('tabix '+bed_file, shell=True)

        print('  * converting cluster coordinates to gene coordinates')
        tss_df = gtf_to_bed(args.genes_gtf)
        cluster2gene_dict = pd.read_csv(os.path.join(args.output_dir, args.prefix + '.leafcutter.clusters_to_genes.txt'),
            sep='\t', index_col=0, squeeze=True).to_dict()

        # add 'chr' prefix
        bed_df['#Chr'] = 'chr'+bed_df['#Chr'].astype(str)
        bed_df['ID'] = 'chr'+bed_df['ID']
        
        if args.coord_mode == "cluster_middle":
            # create a dataframe of cluster starts and ends, indexed by clusterID
            
            cluster_id = bed_df.columns[3]
            cluster_df = bed_df[cluster_id].str.split(":", n = 4, expand = True)
            print(cluster_df)
            #cluster_df.to_csv("test_cluster_df.tsv", sep='\t', header = True)
            # for each junction extract start and end
            # group by clusterID and summarise min(start) and max(end)
            cluster_df.columns = ["chr", "start", "end", "clusterID"]
            cluster_df = cluster_df.astype({"chr": str, "start": int, "end": int, "clusterID": str})
            cluster_grp_df = cluster_df.groupby("clusterID").agg({'chr': 'first', 'start':min, 'end':max}) 
            # convert start and end to middle-1, middle
            cluster_grp_df['end'] = (cluster_grp_df['start'] + ( ( cluster_grp_df['end'] - cluster_grp_df['start'] ) / 2 ) ).astype(int)
            cluster_grp_df['start'] = cluster_grp_df['end'] - 1
            #cluster_grp_df.to_csv("test_cluster_grp_df.tsv", sep='\t', header = True) 
        # debugging - pickle bed_df
        #filehandler = open("bed_df.pickle", mode = 'wb')
        #pickle.dump(bed_df, filehandler)
        print('    ** assigning introns to gene mapping(s)')
        n = 0
        gene_bed_df = []
        group_s = {}
        for _,r in bed_df.iterrows():
            s = r['ID'].split(':')
            cluster_id = s[0]+':'+s[-1]
            if cluster_id in cluster2gene_dict:
                gene_ids = cluster2gene_dict[cluster_id].split(',')
                strand = cluster_id.split("_")[2]
                for g in gene_ids:
                    gi = r['ID'] + ':' + g
                    if args.coord_mode == "cluster_middle":
                         # extract clusterID and match on cluster_grp_df with cluster midpoint coords
                         clu = s[3]
                         gene_bed_df.append(cluster_grp_df.loc[clu, ['chr', 'start', 'end']].tolist() + [gi] + [g] + [strand] + r.iloc[4:].tolist())   
                    
                    if args.coord_mode == "TSS":
                        gene_bed_df.append(tss_df.loc[g, ['chr', 'start', 'end']].tolist() + [gi] + [g] + [strand] + r.iloc[4:].tolist())                    
            else:
                n += 1
        if n>0:
            print('    ** discarded {} introns without gene mapping'.format(n))
        print('  * writing FastQTL inputs')
        # set column names - columns now include "gid" and "strand"
        QTLtools_columns = list(bed_df.columns[0:4]) + ["gid"] + ["strand"] + list(bed_df.columns[4:])
        gene_bed_df = pd.DataFrame(gene_bed_df, columns=QTLtools_columns)
        print( gene_bed_df.head() ) 
        #gene_bed_df.to_csv("test_gene_bed_df.tsv", sep ="\t", header=True)
        #coordinate sort
        gene_bed_df = gene_bed_df.groupby('#Chr', sort=False, group_keys=False).apply(lambda x: x.sort_values('start'))
        # change sample IDs to participant IDs
        # this code was mangling my participant IDs
        #gene_bed_df.rename(columns={i:'-'.join(i.split('-')[:2]) for i in gene_bed_df.columns[4:]}, inplace=True)
        #write_bed(gene_bed_df, os.path.join(args.output_dir, args.prefix+'.leafcutter.bed'))
        #pd.Series(group_s).sort_values().to_csv(os.path.join(args.output_dir, args.prefix+'.leafcutter.phenotype_groups.txt'), sep='\t', header=True)
       
        gene_bed_df.rename(columns = {'#Chr': 'chr', 'ID': 'feature', 'gid': 'group'}, inplace = True)
         
        # for MMQTL - make separate matrix and metadata
        pheno_meta = gene_bed_df[list(gene_bed_df.columns[0:6]) ]  
        pheno_matrix = gene_bed_df[ ["feature"] + list(bed_df.columns[4:] ) ]
        
        print( pheno_matrix.head() )
        print( pheno_meta.head() )

        pheno_matrix.to_csv( os.path.join( args.output_dir, args.prefix + ".leafcutter.phenotype_matrix.tsv"), sep = "\t", header = True, index=False)
        pheno_meta.to_csv(  os.path.join( args.output_dir, args.prefix + ".leafcutter.phenotype_meta.tsv"), sep = "\t", header = True, index=False)
        
        

        #print('  * calculating PCs')
        #pca = PCA(n_components=args.num_pcs)
        #pca.fit(bed_df[bed_df.columns[4:]])
        #pc_df = pd.DataFrame(pca.components_, index=['PC{}'.format(i) for i in range(1,args.num_pcs + 1)],
        #    columns=bed_df.columns[4:] )
        #columns=['-'.join(i.split('-')[:2]) for i in bed_df.columns[4:]])
        # pc_df.index.name = 'ID'
        #pc_df.to_csv(args.prefix+'.leafcutter.PCs.txt', sep='\t', header = True)
        # clean up!
        shutil.rmtree(junc_dir)
        sorted_files = glob.glob("*sorted.gz")
        misc_files = glob.glob( os.path.join(os.path.dirname(args.prefix), "_perind.counts.filtered.gz.phen_*" ) )
        for fname in bed_files + sorted_files + misc_files:
            if os.path.isfile(fname):
                os.remove(fname)    
print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] done')


