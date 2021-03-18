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
	int found;
	char name[100]="";
	int pos, posNorm;
	int linesProcessed;
	int start, end, length;
	int start_search_pos;
	char filename[100];
	int already_sorted;

	struct feature{
		char name[100];
		int start;
		int end;
		int count;
	};
	
	struct feature *f;
	
	struct index{
		int pos;
		char name[50];
	};
	
	struct index *fi;
	
	/* OPEN FILES, INITIALIZE*/
	
	if(argc != 2){
		printf("USAGE: count_features [name of SAM input file]\n");
		printf("Requires one input SAM file as arguement !\nTerminating.\n");
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

	strcpy(filename,"output_");
	strcat(filename, argv[1]);
	filename[strlen(filename)-3]='\0';
	strcat(filename,"txt");
	if(NULL == (out = fopen(filename,"w")) ){
    	printf("Error opening output file: %s\n", filename);
    	printf("Terminating\n");
    	exit(1);
	}
	

	
	
	// count how many features
	num_features=0;
	while(fscanf(in_features,"%*s %*i %*i\n") != EOF){
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
    for(i=0; i<num_features; i++){
    	fscanf(in_features,"%99s\t%i\t%i", f[i].name, &f[i].start, &f[i].end);
		f[i].count = 0;
	}



	// sort features alphabetically and index
	if(num_features<100){
		printf("[ERROR] This version of the program uses indexes and requires at least 100 features!\n");
		printf("[ERROR] Terminating !");
		exit(1);
	}
	printf("[INFO] Sorting features alphabetically...\n");
	fflush(stdout);
	
	already_sorted = 1; // check if already sorted
	i = 0;
	while(i<(num_features-1) && already_sorted){
		if(strcmp(f[i].name, f[i+1].name)>0) {
			already_sorted = 0;
			printf("%s - %s\n",f[i].name, f[i+1].name);
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
	while(fscanf(in_SAM,"%*s %*s %*i %49s %i", name, &i)!=EOF){
		while( (c=getc(in_SAM)) != 10 )
			; // go to end of line
        
        if(strcmp(name,"*")){ // * means it didn't align to anything
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
                if(i>=f[j].start && i<=f[j].end){ //read is within the feature
					f[j].count += 1;
				}          
            }
        }
		linesProcessed += 1;
		if(linesProcessed>10000){
			printf(".");
			fflush(stdout);
			linesProcessed = 0;
		}
	}
	printf("\n[INFO] Done.\n");


	// Output results
	printf("[INFO] Writing output to file.\n");
	for(i=0; i<num_features; i++){
		fprintf(out, "%s\t%i\n", f[i].name, f[i].count);
	}
	

	/* CLOSE AND EXIT */

	fclose(in_SAM);
	fclose(in_features);
	fclose(out);
}

