''' Script to split-transform-merge for pipeline '''

#Write a split-transform-merge pipeline
#Split file by chromosome
#Count the number of transcripts for each chromosome
#Read all count files, calculate the average and write to file

#@
#input fastq file
import gzip
from ruffus import *
from cgatcore import pipeline as P
import sys

@split('genes.gtf.gz', 'chr*.gtf.gz')
def split_chrom(infile, outfiles):
    with gzip.open(infile, 'rt') as inf:
        last_chrom = ''
        for line in inf:
            chrom = line.split()[0]

            if chrom == last_chrom:
                outf.write(line)
            else:
                if last_chrom != '':
                    outf.close()
                outfile = chrom+ '.gtf.gz'
                outf = gzip.open(outfile, 'wt')
                outf.write(line)
                last_chrom = chrom
                print(chrom)

#split_chrom('test.gtf.gz', "")

@transform(split_chrom, suffix('.gtf.gz'), '.count')
def count_genes(infile, outfile):
    statement = "wc -l %(infile)s > %(outfile)s"
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv='obds-py3')

@merge(count_genes, 'average')

def average(infiles, outfile):
    sum = 0
    average = 0
    for infile in infiles:
        with open(infile) as inf:
            readfile = inf.read()
            print(readfile)
            count = int(readfile.split()[0])
            sum += count
            print(sum)
    average = sum/len(infiles)
    print(average)
    with open(outfile, 'wt') as outf:
        outf.write(f'{average}')

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
