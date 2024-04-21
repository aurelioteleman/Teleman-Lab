/*  ******************
 *	v2 - improved the feature search function by using indexing
 *	v3 - improved the indexing
 *
 *  ******************* */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define BUCKETS 20


int main(int argc, char* argv[]){

	FILE *in_SAM, *in_features, *out;
	int num_features;
	int i, j, k, n;
	char c;
	int found;
	char name[50]="";
	int pos, posNorm;
	int linesProcessed;
	int start, end, length;
	int start_search_pos;
	int reads[BUCKETS];
	int totalReads, eightyPercent;
	char outputFilename[100];

	struct feature{
		char name[50];
		int start;
		int end;
		int reads[BUCKETS];
	};
	
	struct feature *f;
	
	
	struct index{
		int pos;
		char name[50];
	};
	
	struct index *fi;

	
	/* OPEN FILES, INITIALIZE*/

	if(argc != 2){
		printf("[ERROR] Need to supply input SAM file name as first arguement.\n");
		printf("[ERROR] Terminating.\n");
		exit(1);
	}
	if(NULL == (in_SAM = fopen(argv[1],"r")) ){
   	printf("Error opening input file: %s! \n", argv[1]);
    	exit(1);
	}
	
	strcpy(outputFilename, "output_");
	strcat(outputFilename, argv[1]);
	outputFilename[strlen(outputFilename)-4]='\0';
	strcat(outputFilename,".txt");
	
	if(NULL == (in_features = fopen("in_features.txt","r")) ){
    	printf("Error opening input file: \"in_features.txt\"! \n");
    	exit(1);
	}

	if(NULL == (out = fopen(outputFilename,"w")) ){
    	printf("Error opening output file: %s! \n", outputFilename);
    	exit(1);
	}
	

	
	
	// count how many features
	num_features=0;
	while(fscanf(in_features,"%*s\t%i\t%i", &i, &j) != EOF){
		num_features += 1;
	}
	rewind(in_features);
	printf("[INFO] Found %i features.\n", num_features);
	
	
	// input feature data into memory
	printf("Input feature data into memory (requires %lu MB of memory) [y/n]? ", sizeof(struct feature)*num_features/1048576);
	scanf("%c", &c);
	if(c != 'y'){
		printf("[ERROR] Terminating due to user request!\n");
		fclose(in_SAM);
		fclose(in_features);
		fclose(out);
		exit(1);
	}
	
	f=malloc(sizeof(struct feature) * num_features);
	if(f==NULL){
		printf("[ERROR] Could not allocate memory for features !\n Terminating !\n");
		fclose(in_SAM);
		fclose(in_features);
		fclose(out);
		exit(1);
	}
	n = 0;
	while(fscanf(in_features,"%49s\t%i\t%i", f[n].name, &f[n].start, &f[n].end) != EOF){
		for(i=0; i<BUCKETS; i++){
			f[n].reads[i]=0;
		}
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
	for(i=0; i<(num_features-1); i++){
		for(j=i+1; j<num_features; j++){
			if(strcmp(f[i].name, f[j].name)>0){ //feature i should come after j
				strcpy(name, f[i].name);
				start = f[i].start;
				end = f[i].end;
				
				strcpy(f[i].name, f[j].name);
				f[i].start = f[j].start;
				f[i].end = f[j].end;
				
				strcpy(f[j].name, name);
				f[j].start = start;
				f[j].end = end;
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
//	while(fscanf(in_SAM,"%*s %*s %*i %49s %*i %i", name, &i)!=EOF){ //NM version
	while(fscanf(in_SAM,"%*s %*s %*i %49s %i", name, &i)!=EOF){ //ENST version
		while( (c=getc(in_SAM)) != 10 && c!=13 && c!=EOF); // go to end of line
		//find position relative to feature start
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
			if(i>=f[j].start && i<f[j].end){
				pos = (i-f[j].start) * BUCKETS / (f[j].end-f[j].start);
				f[j].reads[pos] += 1;
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
	printf("[INFO] Figuring out which transcripts have footprint artefacts and writing output to file...\n");
	fprintf(out, "List of transcripts containing footprint artefacts:\n");
	for(i=0; i<num_features; i++){
		fprintf(out,"%s", f[i].name);
		for(j=0; j<BUCKETS; j++){
			fprintf(out, "\t%i", f[i].reads[j]);
		}
		totalReads = 0;
		for(j=0; j<BUCKETS; j++) totalReads += f[i].reads[j];
		eightyPercent = totalReads*8/10;
		fprintf(out, "\t%i", totalReads);
		
		//find top two buckets
		j=0;
		k=0;
		for(n=1; n<BUCKETS; n++){
			if(f[i].reads[n]>f[i].reads[j]){
				k=j;
				j=n;
			}
		}
		if((f[i].reads[j] + f[i].reads[k]) > eightyPercent){
			fprintf(out,"\t%s\n", f[i].name);
		} else {
			fprintf(out,"\n");
		}
	}
	

	/* CLOSE AND EXIT */

	fclose(in_SAM);
	fclose(in_features);
	fclose(out);
}

