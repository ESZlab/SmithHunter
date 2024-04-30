import glob
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='filter out smithRNAs that has sharp peaks')
parser.add_argument('--mode', type=str, choices=['score','list','filter'], default='score', help='select mode of running, available are score(default), list and filter')
parser.add_argument('--t_five_score', type=float, default=0.40, help='set the score to be met for peaks at five end, if list mode on(default 0.40)')
parser.add_argument('--t_three_score', type=float, default=0.00, help='set the score to be met for peaks at three end, if list mode on(default 0.00)')
parser.add_argument('--penalty', type=float, default=0.10, help='set the penalty of having peaks with gap in between, default is 0.10')
parser.add_argument('--n_thre', dest='n_thre',type=float, default=0.50, help='set the value for select number of peaks, as for n50 value(default 0.50)')
parser.add_argument('--path_bedfiles', type=str, help='directory in which look for coverage file at five and three end')
parser.add_argument('--path_smith', type=str, help='directory in which look for smith file to filter')
parser.add_argument('--sex', type=str, default='NA', help='specify the sex of the animal you are testing, in order to have a column that note it')

args = parser.parse_args()

#penalty = 0.10
#n_thre = 0.5
#t_five_score = 0.40
#t_three_score = 0

def restrict_range(input):
    start_range = None
    for a in range(0,len(input)):
        if input[a] != 0:
            end_range = a
            if start_range is None:
                start_range = a
    return input[start_range:(end_range+1)]

def n_anta(input):
    decre_input = sorted(input,reverse=True)
    tmp_sum = 0
    peaks = []
    for count in range(0,len(decre_input)):
        tmp_sum = tmp_sum + decre_input[count]
        #print(tmp_sum)
        peaks.append(decre_input[count])
        if tmp_sum > (args.n_thre * sum(input)):
            break
    return count+1, peaks

def find_gappy(input, peaks):
    pos_vector = []
    for peak in peaks:
        pos_vector.append([x for x in range(0,len(input)) if input[x] == peak].pop())
    gappy = max(pos_vector) - min(pos_vector) - 1
    return gappy

def filter_fasta(string, cluster):
    if name.count('clusterid' + str(cluster + '_size')) == 1 :
        return True

#path_bedfiles = 'C:/Users/carli/Downloads/5.2_RuPh_results.clusters.bedfiles.female/5.2_RuPh_results.clusters.bedfiles.female'

all_files = glob.glob(args.path_bedfiles + '*.genomecov.3')

cluster_list = []
for i in all_files:
    cluster_list.append(i.split("/")[-1].split(".")[0])

passati = 0
non_passati = 0

for cluster in cluster_list:
    five_end = []
    three_end = []
    with open(args.path_bedfiles + cluster + ".genomecov.5") as f:
        for i in f:
            five_end.append(int(i.split("\t")[2].split("\n")[0]))
    with open(args.path_bedfiles + cluster + ".genomecov.3") as f:
        for i in f:
            three_end.append(int(i.split("\t")[2].split("\n")[0]))

    five_end_res = restrict_range(five_end)
    n_anta_out_five = list(n_anta(five_end_res))
    five_score = 1/(n_anta_out_five[0] + (args.penalty*find_gappy(five_end_res,n_anta_out_five[1])))
    three_end_res = restrict_range(three_end)
    n_anta_out_three = list(n_anta(three_end_res))
    three_score = 1 / (n_anta_out_three[0] + (args.penalty * find_gappy(three_end_res, n_anta_out_three[1])))

    if args.mode == 'list':
        if five_score > args.t_five_score and three_score > args.t_three_score:
            print(str(cluster) + '\t' + 'y' + '\t' + args.sex)
            passati += 1
        else:
            print(str(cluster) + '\t' + 'n' + '\t' + args.sex)
            non_passati += 1
    elif args.mode == 'score':
        print(str(cluster) + '\t' + str(five_score) + '\t' + str(three_score) + '\t' + args.sex)
    else:
        fasta_sequences = SeqIO.parse(open(args.path_smith),'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            if five_score > args.t_five_score and three_score > args.t_three_score:
                if filter_fasta(name,cluster):
                    with open(str("/".join(args.path_smith.split("/")[:-1]) + '/') + 'sharp_smithRNAs.fa','a') as sharp_f:
                        sharp_f.write('>' + str(name) + '\n')
                        sharp_f.write(sequence + '\n')
                        sharp_f.close()
