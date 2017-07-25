Preprocessing Data
===================

In this exercise, we will learn how to preprocess our data for alignment. We will be doing some sequencing trimming. Trimming the left side of the reads to remove primer sequence as well as identify and trim off any adapter, or polyA/T sequence from the right edge of the read

**1\.** First, create a directory for the example in your home directory:

    cd # Go back to the home directory
    mkdir rnaseq_workshop 

---

**2\.** Next, go into that directory, create a new directory for the Raw Data and link to the directory for your raw data. This data comes from an Yeast project that we did:

    cd rnaseq_workshop
    mkdir 00-RawData
    cd 00-RawData
    ln -s path/to/data/folder .

---

**3\.** Now, take a look inside that directory.

    ls

--- 

**4\.** Pick a file and look at the contents of the file using the 'zless' command:

    zless the_file_you_choose.fastq.gz

Make sure you can identify which lines correspond to a single read and which lines are the header, sequence, and quality values. Press 'q' to exit this screen. Then, let's figure out the number of reads in this file. A simple way to do that is to count the number of lines and divide by 4 (because the record of each read uses 4 lines). In order to do this, use "zcat" to output the uncompressed file and pipe that to "wc" to count the number of lines:

    zcat the_file_you_choose.fastq.gz | wc -l

Divide this number by 4 and you have the number of reads in this file. One more thing to try is to figure out the length of the reads without counting each nucleotide. First get the first 4 lines of the file (i.e. the first record):

    zcat the_file_you_choose.fastq.gz | head -4

Then, copy and paste the sequence line into the following command (replace [sequence] with the line):

    echo -n [sequence] | wc -c

This will give you the length of the read. See if you can figure out how this command works.

---

**5\.** Now go back to your 'rnaseq_workshop' directory and create another directory called '01-Trimming':

    cd ~/rnaseq_workshop
    mkdir 01-Trimming

---

**6\.** Now, we will use software called 'bbduk' to trim the data. First we will download and install it, then run it on just one file. 

find bbdup on the internet, download, extract and move it into the 'rnaseq_workshop' folder.

Next, Looking at the usage. 

    ./bbmap/bbduk.sh -h

**7\.** In order to trim adapters we need an adapter file and the sequence file. The adapter file will depend upon which kit you used... typically you can find the adapters from the sequencing provider. In this case, Illumina TruSeq adapters were used, adapter files can be found in the resources folder.

To trim polyA/T we will need to create a file with a polyA and polyT sequence.

    printf ">polyA\nAAAAAAAAAAAAA\n>polyT\nTTTTTTTTTTTTT\n" | gzip - >  bbmap/resources/polyA.fa.gz

**8\.** Now run bbduk specifying an input file, output file, the adapters and polyA/T files, and additional parameters for trimming.

    bbmap/bbduk.sh in=00-RawData/the_file_you_choose.fastq.gz out=01-Trimmed/the_file_you_choose.fastq.gz ref=bbmap/resources/polyA.fa.gz,bbmap/resources/truseq.fa.gz k=13 ktrim=r forcetrimleft=11 useshortkmers=t mink=5 qtrim=t trimq=10 minlength=20

the parameters: k, ktrim, useshortkmers, mink, qtrim, trimq, are all parameters associated with trimming off the 

This will take a few minutes to run. When the jobs finish, you will have the file that is  trimmed.

---

**9\.** Once that is done, let's take a look at the differences between the input and output files.

---

**10\.** We have run through adapter & quality trimming for one sample, but now we need to do it for all the samples. 
    
We will need to generate a file called 'samples.txt' that contains all the sample IDs. We will extract this information from the filenames using 'cut'. First, get a listing of all the R1 files:

    ls -1 *R1*.fastq.gz

We will pipe this output to 'cut' to get the fields we want. Give cut the options "-d" with an underscore (usually above the minus sign on a keyboard) as the parameter, and the "-f" option with "1,2" as the parameter in order to get the first and second fields:

    ls -1 *R1*.fastq.gz | cut -d_ -f1,2

This gives us the all the sample IDs. Now we just need to redirect that output to a file:

    ls -1 *R1*.fastq.gz | cut -d_ -f1,2 > samples.txt

Use 'cat' to view the contents of the file to make sure it looks right:

    cat samples.txt

---

**14\.** Take a look at the other script:

    cat qa_for_loop.sh

This script has similar commands, but instead of using a task array, it is using a for loop. So this will loop through all the IDs in samples.txt and assign a new ID on every iteration of the loop. You should use this script if you will be running jobs NOT on a cluster, but on a single machine.

