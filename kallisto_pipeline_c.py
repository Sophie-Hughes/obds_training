'''
Kallisto Pipeline
Steps:
1. Copy and paste fastq/multiqc from previous pipeline
2. Update yml to fit params
3. Do pseudoalignment with kallisto
4. Check quality of pseudoalignment
- multiqc will pick out kallisto files and produce a report
'''

from ruffus import *
from cgatcore import pipeline as P
import sys
import gzip

params = P.get_parameters("kallisto_pipeline_c.yml")

@follows(mkdir("fastqc"))
@transform("*.fastq.gz", regex(r'(.*).fastq.gz'), r'fastqc/\1_fastqc.html')
def fastqc(infile, outfile):
    statement = "fastqc --nogroup -o fastqc %(infile)s "
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv='obds-py3')

@merge(fastqc, r'fastqc/multiqc_report.html')
def multiqc(infiles, outfile):
    statement = "multiqc -f -n %(outfile)s fastqc"
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv='obds-py3')

#pseudoalignment index
@follows(fastqc, multiqc)
@follows(mkdir("kallisto"))
@transform(params["fasta_file"], regex(r'(.+).fa.gz'), r'kallisto/\1.idx')
def kallisto_index(infile, outfile):
    statement = '''kallisto index -i %(outfile)s %(infile)s'''
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv='obds-py3')


#pseudoalignment quantification
@follows(mkdir("kallisto_quant"))
@collate("*.fastq.gz", regex(r'(.*)_[12](.*).fastq.gz'), r'kallisto_quant/\1/abundance.tsv')
def kallisto_quant(infiles, outfile):
    read1,read2 = infiles
    statement = '''kallisto quant
                -t %(kallisto_quant_threads)s
                --index %(kallisto_quant_index)s
                --bootstrap = 0
                %(kallisto_quant_option)s
                -o %(outdir)s
                %(read1)s %(read2)s'''
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv='obds-py3')


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
