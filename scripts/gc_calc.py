import io
import os
from Bio import  SeqIO
import numpy as np
import sys
from subprocess import Popen, PIPE
from multiprocessing import Pool
class GC:
    def __init__(self, ref=False):
        self.chroms = [f"chr{i}" for i in range(1, 23)]
        self.record_holder = {}
        if ref:
            for record in SeqIO.parse("hg38.fa", 'fasta'):
                if record.id in self.chroms:
                    self.record_holder[record.id] = record.seq


    def generateReference(self, nthreads=23):

        with Pool(nthreads) as myProcess:
            myProcess.map(self.fragWindow, self.chroms)


    def fragWindow(self, chromosome):

        fraglen=400
        seq = self.record_holder[chromosome]
        chrom_array = np.full(fill_value=np.nan, shape=(len(seq) - fraglen, 2))  # prefer to use a numpy array here for speed trade-off vs memory usage
        with open('hg38.uniquely.mappable.bedgraph', 'r') as myBed: #open bed file contain mappability info
            for line in myBed:
                chrom, pos1, pos2, mapscore = line.split('\t')
                #get the bounds of each window where there is uniquely mappable sequence then get GC % in each window
                if chrom == chromosome and int(mapscore) == 1: #checks the mapscore and chrom name
                    for bp in range(int(pos1)-1, int(pos2)):
                        seq_frag = np.asarray(list(str(seq[bp:bp + fraglen])))

                        #b/c the genome fasta has large stretches of Ns that will resolve as 0% GC regions and be unmappable we instead will exclude Ns from the analysis by removing any window with an N
                        #conservative approach but keeps window size constant
                        #genome is soft masked so we are going to exclude positions that have soft-masked positions as well
                        if np.sum(seq_frag == "N") == 0 and np.sum(seq_frag == "a") == 0 and np.sum(seq_frag == "t") == 0 and np.sum(seq_frag == "g") == 0 and np.sum(seq_frag == "c") == 0:
                            gc = (np.sum(seq_frag == "G") + np.sum(seq_frag == "C")) / fraglen
                        else:
                            gc = np.nan

                        #do not need the first vector of this array anymore need to remove and re-factor
                        chrom_array[:,0][bp] = bp+1
                        chrom_array[:,1][bp] = gc

        np.save(f"hg38.{chromosome}.gc.npy", chrom_array)


    def gc_calculation(self, sample_name, bins=100):
        cwd = os.getcwd()
        gc_bins = np.asarray([ ((b+1)/bins) for b in range(bins)])
        gc_table = {((b+1)/bins):[0,0] for b in range(bins)} #a table where the value of each key is [read_depth, # of positions] so we can compute an avg depth | GC
        for chromosome in self.chroms: #iter through all chroms
            gc_input = np.load(os.path.join(cwd,"GC_UNIQ", f"hg38.{chromosome}.gc.npy"))
            chr_end = gc_input.shape[0]

            query = '{}:{}-{}'.format(chromosome, 1, chr_end)
            with Popen(['tabix', sample_name, query], stdout=PIPE, stderr=None) as process:

                #for stream in process.communicate()[0].decode("utf-8").split("\n"):
                # #stream the tabix process in line by line to save memory
                for stream in io.TextIOWrapper(process.stdout, encoding="utf-8"):
                    #stream = stream.decode('utf-8')

                    if len(stream) > 1: #sanitize empty string at end
                        chrom, pos1, pos2, depth = stream.split("\t") #newlines are ignored when converting to floats -- python magic!!!

                        pos1 = int(pos1)
                        pos2 = int(pos2)
                        pos1 = max(0, (pos1-1) )#account for 0 indexing
                        pos2 = pos2-1
                        gc = gc_input[pos1:pos2] #we will grab GC from windows, just a slice so should not use new memory assignments

                        for g in gc: #would like to speed this up but it is what it is -- this still should be o(n) even though two nested for loops
                            if ~np.isnan(g): #skip all positions with all Ns -- I think we could do better here, by removing NaNs before hand
                                fit_bin = min(gc_bins[float(g) <= gc_bins])
                                gc_table[fit_bin][0] += float(depth)
                                gc_table[fit_bin][1] += 1
                            else:
                                pass
            sys.stderr.write(f"{chromosome} complete\n")
            gc_input = None
        for element in gc_bins:
            try:
                avg_depth = gc_table[element][0]/gc_table[element][1]
            except ZeroDivisionError: #accomodate empty bins for error case
                avg_depth = 0
            print(f'{element}\t{gc_table[element][0]}\t{gc_table[element][1]}\t{avg_depth}')

if __name__ == '__main__':
    fname = sys.argv[1]


    if sys.argv[2] == "gc":
        invoke_GC = GC(ref=False)
        invoke_GC.gc_calculation(sample_name=fname)
    elif sys.argv[2] == 'ref':
        invoke_GC = GC(ref=True)
        invoke_GC.generateReference()
