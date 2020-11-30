'''
RNAseq pipeline
process fastq file into count files/matrices
'''

from ruffus import *
from cgatcore import pipeline as P
import sys

params = P.get_parameters("rnaseq_pipeline.yml")

@follows(mkdir("fastqc"))
@transform("*.fastq.gz", regex(r'(.*).fastq.gz'), r'fastqc/\1_fastqc.html')
def fastqc(infile, outfile):
    statement = "fastqc --nogroup -o fastqc %(infile)s "
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv='obds-py3')


@merge(fastqc, r'fastqc/multiqc_report.html')
def multiqc(infiles, outfile):
    statement = "multiqc -f -n %(outfile)s fastqc"
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv='obds-py3')

@follows(mkdir("bam"))
#@collate("*.fastq.gz", regex(r'(.*)_[12]_test.fastq.gz'), r'bam/\1.bam')
@collate("*.fastq.gz", regex(r'(.*)_[12](.*).fastq.gz'), r'bam/\1.bam')
def hisat2(infiles, outfile):
    read1,read2 = infiles
    statement = '''hisat2 -x %(hisat2_index)s -1 %(read1)s -2 %(read2)s %(hisat2_options)s --threads %(hisat2_threads)s --summary-file %(outfile)s.log
    | samtools sort - -o %(outfile)s
    && samtools index %(outfile)s'''
    # End up with a sorted SAM file which then ?we turn into a BAM file?
    P.run(statement, job_queue='all.q', job_threads=params["hisat2_threads"], job_memory='1G', job_condaenv='obds-py3')

@transform (hisat2, regex(r'(.*).bam'), r'\1.idxstat')
def idxstats(infile, outfile):
    statement = "samtools idxstats %(infile)s > %(outfile)s"
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv='obds-py3')

@transform (hisat2, regex(r'(.*).bam'), r'\1.flagstat')
def flagstat(infile, outfile):
    statement = "samtools flagstat %(infile)s > %(outfile)s"
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv='obds-py3')

@merge(hisat2, 'readcount.tsv')
def feature_counts(infiles, outfile):
    infiles_string = ' '.join(infiles)
    statement = '''featureCounts %(featurecounts_options)s
                                -T %(featurecounts_threads)s
                                -s %(featurecounts_strandness)s
                                -a %(featurecounts_annotation)s
                                -o %(outfile)s %(infiles_string)s'''
    P.run(statement, job_queue='all.q', job_threads=params["featurecounts_threads"], job_memory='2G', job_condaenv='obds-py3')

@merge([idxstats, flagstat, feature_counts], 'bam/bamqc.html')
def bam_qc(infiles, outfile):
    statement = "multiqc -f -n %(outfile)s bam"
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv='obds-py3')


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
