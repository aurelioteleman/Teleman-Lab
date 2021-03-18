#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


int main(int argc, char* argv[]){

	FILE *in_SAM, *in_features, *out;
	int num_features;
	int max_feature_len;
	int i, j, n;
	char c;
	int *counts, *countsNorm;
	int found;
	char name[50]="";
	char read[150];
	int pos, posNorm;
	int linesProcessed;
	int start, end, length;
	int start_search_pos;
	int read_length;
	int already_sorted;

	struct feature{
		char name[50];
		int pos;
	};
	
	struct feature *f;
	
	struct index{
		int pos;
		char name[50];
	};
	
	struct index *fi;
	
	/* OPEN FILES, INITIALIZE*/
	

	if(NULL == (in_SAM = fopen("input.sam","r")) ){
    	printf("Error opening input file: \"input.sam\"! \n");
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
	

	
	
	// count how many features
	num_features=0;
	while(fscanf(in_features,"%*s\t%*i") != EOF){
		num_features += 1;
	}
	rewind(in_features);
	printf("[INFO] Found %i features.\n", num_features);
	
	
	// input feature data into memory
	
	f=malloc(sizeof(struct feature) * num_features);
	if(f==NULL){
		printf("[ERROR] Could not allocate memory for features !\n Terminating !\n");
		fclose(in_SAM);
		fclose(in_features);
		fclose(out);
		exit(1);
	}
	n = 0;
	while(fscanf(in_features,"%49s\t%i", f[n].name, &f[n].pos) != EOF){
		n++;
	}




	// sort features alphabetically and index
	if(num_features<100){
		printf("[ERROR] This version of the program uses indexes and requires at least 100 features!\n");
		printf("[ERROR] Terminating !");
		exit(1);
	}
	printf("[INFO] Sorting features alphabetically...");
	fflush(stdout);
	
	already_sorted = 1; // check if already sorted
	i = 0;
	while(i<(num_features-1) && already_sorted){
		if(strcmp(f[i].name, f[i+1].name)>0) {
			already_sorted = 0;
			}
		i = i+1;
	}

	if(already_sorted)	{
		printf("[INFO] Features already alphabetically sorted.\n");
	} else {
		printf("[INFO] Features not already alphabetically sorted.\n");	
		for(i=0; i<(num_features-1); i++){
			for(j=i+1; j<num_features; j++){
				if(strcmp(f[i].name, f[j].name)>0){ //feature i should come after j
					strcpy(name, f[i].name);
					pos = f[i].pos;
				
					strcpy(f[i].name, f[j].name);
					f[i].pos = f[j].pos;
				
					strcpy(f[j].name, name);
					f[j].pos = pos;
				}
			}
		}
	}
	printf("Done.\n");
	
	
	printf("[INFO] Creating feature index...");
	fi = malloc(sizeof(struct index)*101);
	if(fi==NULL){
		printf("[ERROR] Could not allocate memory for feature index!\n");
		printf("[ERROR] Terminating !\n");
		exit(1);
	}
	for(i=0; i<100; i++){
		fi[i].pos = i*num_features/100;
		strcpy(fi[i].name, f[i*num_features/100].name);
	}
	fi[100].pos = num_features;
	strcpy(fi[100].name, "");
	printf("Done.\n");
	
	
	
	// allocate memory for counts and set all to zero
	counts=malloc(sizeof(int) * 200 * 100); // separately count footprints of lengths 1-100 nt
	if(counts==NULL){
		printf("[ERROR] Could not allocate memory for counts !\n Terminating !\n");
		fclose(in_SAM);
		fclose(in_features);
		fclose(out);
		exit(1);
	}

	// start codon will be at i=99 (ie subtract 99 from i to obtain position)
	// index i = (position relative to ATG + 99) + (200*(footprint length-1)) 
	for(i=0; i<200*100; i++){ 
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
	printf("Hundred-thousand lines processed: ");
	fflush(stdout);
	linesProcessed=0; //use this to provide an update on progress
//	while(fscanf(in_SAM,"%*s %*s %*i %49s %*i %i %*s %*s %*s %*s %*s %149s", name, &i, read)!=EOF){
	while(fscanf(in_SAM,"%*s %*s %*i %49s %i %*s %*s %*s %*s %*s %149s", name, &i, read)!=EOF){
		while( (c=getc(in_SAM)) != 10  && c!=EOF); // go to end of line
		
		if(strcmp(name,"*") && strcmp(read,"*")){  // excludes unmapped reads and secondary mappings, which have a flag of 256 and * as read sequence
			//calculate read length
			read_length = strlen(read);
			if(read_length < 1) read_length = 1;
			if(read_length > 100) read_length = 100;

			//find position relative to CDS start
			found = 1; //not yet found
		
			// search the index to see where to start searching
			n=0;
			while(n<99 && strcmp(fi[n+1].name, name)<=0){
				n++;
			}
		
			for(j=fi[n].pos; j<fi[n+1].pos && found; j++){
				if(!strcmp(f[j].name, name)){ //found it !
					found = 0;
				}
			}
			if(found == 0){
				j = j-1;
				while(!strcmp(f[j].name, name)){	
					if(i-f[j].pos > -100 && i-f[j].pos < 101){ //read is within the feature
						pos = i - f[j].pos + 99;
						counts[pos + 200*(read_length-1)] += 1;
					}
					j=j+1;
				}
			}
		}
		linesProcessed += 1;
		if(linesProcessed>100000){
			printf(".");
			fflush(stdout);
			linesProcessed = 0;
		}
	}
	printf("\n[INFO] Done.\n");


	// Output results
	printf("[INFO] Writing output to file.\n");
	fprintf(out,"position-footprint length:");
	for(i=1; i<=100; i++){
		fprintf(out,"\t%i", i);
	}
	fprintf(out,"\n");
	for(i=0; i<200; i++){
		fprintf(out,"%i",i-99);
		for(j=1; j<101; j++){
			fprintf(out, "\t%i", counts[i+200*(j-1)]);
		}
		fprintf(out,"\n");
	}
	

	/* CLOSE AND EXIT */

	fclose(in_SAM);
	fclose(in_features);
	fclose(out);
}

