<h1>Demo for transcripts_with_artefacts_v4_20_buckets</h1>

<h2>Instructions</h2>

1. Download the source code: transcripts_with_artefacts_v4_20_buckets.c

2. Compile the source code using a C compiler such as gcc

   ```
   gcc transcripts_with_artefacts_v4_20_buckets.c -o transcripts_with_artefacts_v4_20_buckets
   ```

3. Download the two input files "in_features.txt" and "input.sam" that are provided in this demo folder and place them in the same folder as the compiled program.

4. Run the program using the command

   ```
   ./transcripts_with_artefacts_v4_20_buckets input.sam
   ```
   
5. The program will ask you whether to load all features into memory, and reply 'y'

6. The program will run in circa 1 second and terminate, creating an output file "output_input.txt"





<h2>Expected Output</h2>

The expected output file "output_input.txt" should be identical to the one provided in this demo folder. It contains multiple columns. First column - a list of all features listed in the "in_features.txt" file. Columns 2 to 11 - the number of reads in the input.sam file that land within a position range along the transcript, divided up into 10 decile buckets. Column 12 - if more than 80% of all reads for a transcript land in only 2 buckets, the transcript gets flagged as having a PCR artefact. Otherwise, the 12th column is blank.





<h2>Expected run time for the demo</h2>

1 minute





<h2>Notes</h2>

The input SAM file can vary in format depending on the mapping program that was used to map the  reads, and the sequences that the reads were mapped against. In the example provided here, "Query template name" (the identifier of each read in the first column of the SAM file) is made of two strings separated by a space, rather than one string. Furthermore the "Reference sequence name" (transcript ID in the 3rd column) is also made of two strings separated by a space, because the original FASTA file with transcript sequences had a version number for each transcript, and the mapper that was used to generate this SAM file separated out that version number. The key code line that determines how count_features_nt_range_v2 reads the .sam input file is line 188. Please make sure this matches the format of your SAM file and adjust it as necessary.

