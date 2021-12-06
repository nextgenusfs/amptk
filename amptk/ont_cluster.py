from bisect import bisect_left
import numpy as np
from scipy.signal import find_peaks
import pyfastx
import subprocess


def primary_clustering(fastq, outfolder, cpus=1, min_num=3):
    cmd = ['isONclust', '--ont', '--fastq', fastq, '--t', str(cpus),
           '--outfolder', outfolder]
    subprocess.call(cmd)
    cmd2 = ['isONclust', 'write_fastq', '--clusters', os.path.join(outfolder, 'final_clusters.tsv'),
            '--fastq', fastq, '--outfolder', os.path.join(outfolder, 'clust_0'), '--N', str(min_num)]
    subprocess.call(cmd2)


def subcluster(folder, outfolder, n=15, w=4, d=20):
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
    for f in os.listdir(folder):
        if f.endswith('.fastq'):
            fastq_file = os.path.join(folder, f)
            cluster_num = int(f.rstrip('.fastq'))
            reads = []
            fqindex = pyfastx.Fastq(fastq_file)
            for r in fqindex:
                reads.append((r.name, len(r)))
            lengths = [x[1] for x in reads]
            centers = get_centroids_by_length(lengths, n=n, w=w, d=d)
            clusters = {}
            if len(centers) > 0:
                for c in centers:
                    clusters[c] = []
                for read in reads:
                    closest = take_closest(centers, read[1])
                    clusters[closest].append(read[0])
            else:
                clusters[0] = [x[0] for x in reads]
            for k, v in clusters.items():
                fqoutfile = os.path.join(outfolder, '{}_{}.fq'.format(cluster_num, k))
                with open(fqoutfile, 'w') as outfile:
                    for x in v:
                        outfile.write('{}'.format(fqindex[x].raw))


def take_closest(myList, myNumber):
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return after
    else:
        return before

def get_centroids_by_length(lengths, n=15, w=4, d=20):
    s_lengths = sorted(lengths)
    r = []
    for i in range(s_lengths[0], s_lengths[-1], 1):
        r.append(s_lengths.count(i))
    x = np.array(r)
    peaks, prop = find_peaks(x, height=n, width=w, distance=d)
    len_centers = []
    if len(peaks) > 0:
        for z in peaks:
            len_centers.append(z+s_lengths[0])
    return len_centers