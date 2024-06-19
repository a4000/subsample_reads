from Bio import SeqIO
import random

# Get a random boolean value
def rand_bool():
    i = random.randint(0, 1)
    return bool(i)

assays = ["12S", "16S"]
nums = ["1", "2", "3", "4", "5"]
max_count = 2000 # This should be equal to the number of sequences per taxa in the simulated reads
random.seed(1)
rand_seq_count = random.randint(1, max_count)
output_sample = rand_bool()
output_set = set()

for assay in assays:
    print("Assay: " + assay)
    for num in nums:
        print("Sample number: " + num)
        fastq_file1_path = "./" + assay + "_simreads" + num + "1.fq"
        fastq_file2_path = "./" + assay + "_simreads" + num + "2.fq"

        # Specify the paths for the output files
        output_file1_path = "./subset_" + assay + "_simreads" + num + "1.fq"
        output_file2_path = "./subset_" + assay + "_simreads" + num + "2.fq"

        # Open the input FASTQ files and output files
        curr_seq_num = 0
        seq_count_i = 0
        with open(fastq_file1_path, 'r') as fastq_file1, open(fastq_file2_path, 'r') as fastq_file2, \
            open(output_file1_path, 'w') as output_file1, open(output_file2_path, 'w') as output_file2:

            for record1, record2 in zip(SeqIO.parse(fastq_file1, 'fastq'), SeqIO.parse(fastq_file2, 'fastq')):
                curr_seq_num += 1

                # There's a 50% chance of outputing the sample in a specifc readpair
                if output_sample:
                    read_name = record1.id.split("-")[0]
                    output_set.add(assay + "," + num + "," + read_name + "," + str(rand_seq_count))
                    # We are also ranomising the count of sequences we are outputing for a given sample
                    if curr_seq_num <= rand_seq_count:
                        SeqIO.write(record1, output_file1, 'fastq')
                        SeqIO.write(record2, output_file2, 'fastq')

                # We've reached the last sequence for the current sample
                if curr_seq_num == max_count:
                    curr_seq_num = 0

                    # We are doing this to make sure there are a decent number of samples with a very low sequence count
                    if rand_bool():
                        rand_seq_count = 1
                    else:
                        rand_seq_count = random.randint(1, max_count)
                    output_sample = rand_bool()

# This isn't needed, but it might be helpful to have a csv file with this information
sorted_set = sorted(output_set)
with open("output.csv", "w") as file:
    file.write("assay,file_num,read_name,sequence_count\n")
    for item in sorted_set:
        file.write(str(item) + "\n")
