#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


int main(int argc, char* argv[]){

	FILE *in_SAM, *in_features, *out, *outDetails;
	int num_features;
	int max_feature_len;
	int i, j, n;
	char c;
	int *counts;
	int found;
	char name[50]="";
	char read[150];
	int pos, posNorm;
	int linesProcessed;
	int start, end, length;
	int start_search_pos;
	int read_length;
	int already_sorted;
	int transcriptLength;


	struct transcript{
		char name[50];
		int totalLength;
		int startCodonPos;
		int stopCodonPos;
	};
	
	struct transcript *f;
	
	struct index{
		int pos;
		char name[50];
	};
	
	struct index *fi;
	
	/* OPEN FILES, INITIALIZE*/
	
	if(argc != 3){
		printf("Usage: metagene_plot_5UTR_CDS_3UTR [sam_filename] [gene to search e.g. NM_000014]\n");
		exit(1);
	}

	if(NULL == (in_SAM = fopen(argv[1],"r")) ){
    	printf("Error opening input file: %s\n", argv[1]);
    	exit(1);
	}
	if(NULL == (in_features = fopen("in_features.txt","r")) ){
    	printf("Error opening input file: \"in_features.txt\"! \n");
    	exit(1);
	}

	if(NULL == (out = fopen("output.txt","w")) ){
    	printf("Error opening output file: \"output.txt\"! \n");
    	exit(1);
	}
	if(NULL == (outDetails = fopen("output-reads.txt","w")) ){
    	printf("Error opening output file: \"outDetails.txt\"! \n");
    	exit(1);
	}
	

	
	
	// count how many features
	num_features=0;
	while(fscanf(in_features,"%*s %*i %*i %*i") != EOF){
		num_features += 1;
	}
	rewind(in_features);
	printf("[INFO] Found %i features.\n", num_features);

	
	// input feature data into memory
	printf("Input feature data into memory (requires %lu MB of memory) [y/n]? ", sizeof(struct transcript)*num_features/1048576);
	scanf("%c", &c);
	if(c != 'y'){
		printf("[ERROR] Terminating due to user request!\n");
		fclose(in_SAM);
		fclose(in_features);
		fclose(out);
		exit(1);
	}
	
	f=malloc(sizeof(struct transcript) * num_features);
	if(f==NULL){
		printf("[ERROR] Could not allocate memory for features !\n Terminating !\n");
		fclose(in_SAM);
		fclose(in_features);
		fclose(out);
		exit(1);
	}
	n = 0;
	while(fscanf(in_features,"%49s\t%i\t%i\t%i", f[n].name, &f[n].totalLength, &f[n].startCodonPos, &f[n].stopCodonPos) != EOF){
		n++;
	}



	// find the transcript length for the gene to be counted
	n = 0;
	while(n<num_features && strcmp(argv[2], f[n].name)) n++;
	if(strcmp(argv[2], f[n].name)){
		printf("[ERROR] Could not find %s in feature list.\n Terminating\'n", argv[2]);
		exit(1);
	}
	transcriptLength = f[n].totalLength;
	printf("[INFO] %s has a length of %i nt.\n", argv[2], transcriptLength);
	
	// allocate memory for counts and set all to zero
	counts=malloc(sizeof(int)*100*transcriptLength);
	if(counts==NULL){
		printf("[ERROR] Could not allocate memory for counts !\nTerminating!\n");
		exit(1);
	}

	// footprint length varies from 1 to 100
	// index i = transcriptLength*(footprint length - 1) + position on transcript
	for(i=0; i<(100 * transcriptLength); i++){
		counts[i]=0;
	}


	// Go to end of header rows
	printf("[INFO] Reading input SAM file.\n");
	printf("[INFO] Looking for end of header rows...\n");
	while( (c=getc(in_SAM)) == '@' ){
		while( (c=getc(in_SAM)) != 10 ); // go to end of line
	}
	printf("[INFO] Done.\n");


	// Count !
	printf("[INFO] Counting...\n");
	printf("Ten-thousands of lines processed: ");
	fflush(stdout);
	linesProcessed=0; //use this to provide an update on progress
	while(fscanf(in_SAM,"%*s %*s %*i %49s %i %*s %*s %*s %*s %*s %149s", name, &i, read)!=EOF){
		while( (c=getc(in_SAM)) != EOF && c!= 10 ); // go to end of line
		
		if(strcmp(name,"*")){  // name=* means it did not map (hence excluded); read=* is a secondary mappings, which have a flag of 256 and * as read sequence
			//calculate read length
			if(strcmp(read,"*")){ //primary mapping - the read sequence is displayed
				read_length = strlen(read);
				if(read_length < 1) read_length = 1;
				if(read_length > 100) read_length = 100;
			} 
			// else {nothing}.. This is a secondary mapping - use the same read_length as the previous read
			
			if(!strcmp(name, argv[2])){ //read is in the transcript of interest
				pos = i-1;
				counts[transcriptLength*(read_length - 1) + pos] += 1;
				fprintf(outDetails,"%s len=%i pos=%i\n", read, read_length,i);
			}

		}
		linesProcessed += 1;
		if(linesProcessed>10000){
			printf(".");
			fflush(stdout);
			linesProcessed = 0;
		}
	}
	printf("[INFO] Done.\n");


	// Output results
	printf("[INFO] Writing output to file.\n");

	// print output file header
	fprintf(out,"position-footprint length:");
	for(i=1; i<=100; i++){
		fprintf(out,"\t%i", i);
	}
	fprintf(out,"\n");
	for(pos=1; pos<=transcriptLength; pos++){
		fprintf(out,"%i",pos);
		for(j=1; j<101; j++){
			fprintf(out, "\t%i", counts[transcriptLength*(j - 1) + pos-1]);
		}
		fprintf(out,"\n");
	}

	

	/* CLOSE AND EXIT */

	fclose(in_SAM);
	fclose(in_features);
	fclose(out);
	fclose(outDetails);
}

