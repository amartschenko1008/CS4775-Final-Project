import glob
def load_file():
    basepath = '/mnt/d'

    fnames_sample = []
    fnames_random = []

    fname_sample = sorted(glob.glob(basepath + '/*1999_20_final_sample.txt', recursive=True))
    fname_random = sorted(glob.glob(basepath + '/*1999_20_final_random.txt', recursive=True))
    fnames_sample.extend(fname_sample)
    fnames_random.extend(fname_random)
    return fnames_sample, fnames_random

def inverted_analysis(fnames, output_file):
    sample_num = 0
    detect_num = 0
    with open(output_file, 'w') as out:
        for exp_id, fname in enumerate(fnames):
            print('Processing:', fname)
            found = False
            detect = False
            repeat_num = []
            cluster_num = 0
            repeat = 0
            with open(fname, 'r') as f:
                for line in f:
                    if not line.strip():
                        if found:
                            if detect:
                                continue
                            else:
                                found = False
                                detect = False
                        continue

                    if line.strip() == 'Palindromes:':
                        sample_num += 1
                        found = True
                        continue

                    if found:
                        if not line.strip():
                            if not detect:
                                found = False
                                detect = False
                        else:
                            if detect == False:
                                detect = True
                                detect_num += 1
                            out.write(line)
                            if '||' in line:
                                repeat += 1

                    if line.strip() == 'Palindromes of: sample_new.fasta' or line.strip() == 'Palindromes of: random_new.fasta':
                        found = False
                        detect = False
                        repeat_num.append(repeat)
                        repeat = 0

                for i in range(len(repeat_num)):
                    if repeat_num[i] > 20:
                        cluster_num += 1

    return sample_num, detect_num, cluster_num         


def main():
    fnames_sample, fnames_random = load_file()
    integrated_sample = '/mnt/d/integrated_sample.txt'
    integrated_random = '/mnt/d/integrated_random.txt'
    sample, inverted_sample, cluster_sam = inverted_analysis(fnames_sample, integrated_sample)
    random, inverted_random, cluster_rnd = inverted_analysis(fnames_random, integrated_random)
    print('sample_num:', sample, 'inverted_sample:', inverted_sample, 'cluster_sam:', cluster_sam)
    print('random_num:', random, 'inverted_random:', inverted_random, 'cluster_rnd:', cluster_rnd)

if __name__ == '__main__':
    main()    
