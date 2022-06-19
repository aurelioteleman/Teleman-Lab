#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define BREAKPOINT_START 5
#define BREAKPOINT_END 93

int main(int argc, char* argv[]){

	FILE *in, *out;
	char transcriptName[200];
	float reads[100];
	float filteredReads[98];
	int i;
	int breakpoint, minBreakpoint;
	float average1, average2, sumOfSquares[100], sos, minSos;
	int maxPos1, maxPos2;
	

	
	/* OPEN FILES, INITIALIZE*/

	if(NULL == (in = fopen("input.txt","r")) ){
   	printf("Error opening input file: input.txt \n");
    	exit(1);
	}
	
	if(NULL == (out = fopen("output.txt","w")) ){
    	printf("Error opening output file: output.txt! \n");
    	exit(1);
	}
	
	fprintf(out,"Transcript");
	for(i=1; i<101; i++){
		fprintf(out,"\tReads_%i", i);
	}
	fprintf(out,"\tBreakpoint\tAverage 1\tAverage 2\tDifference\tSum of Squares\n");
	

	while(fscanf(in, "%s", transcriptName) !=EOF){
		for(i=0; i<100 && fscanf(in, "%f", &reads[i])!=EOF; i++){
			;
		}
		
		//FIND TWO MAXIMUM DATAPOINTS AND REMOVE (OFTEN OUTLIERS)
		if(reads[0]>reads[1]){
			maxPos1 = 0;
			maxPos2 = 1;
		} else {
			maxPos1 = 1;
			maxPos2 = 0;
		}
		for(i=2; i<100; i++){
			if(reads[i]>reads[maxPos1]){		
				maxPos2 = maxPos1;
				maxPos1 = i;
			} else if(reads[i]>reads[maxPos2]){
				maxPos2 = i;
			}
		}
		if(maxPos1>maxPos2){
			i = maxPos1;
			maxPos1 = maxPos2;
			maxPos2 = i;
		}
		for(i=0; i<maxPos1; i++){
			filteredReads[i]=reads[i];
		}
		for(i=maxPos1+1; i<maxPos2; i++){
			filteredReads[i-1]=reads[i];
		}
		for(i=maxPos2+1; i<100; i++){
			filteredReads[i-2]=reads[i];
		}
		
		
		for(breakpoint=BREAKPOINT_START; breakpoint<BREAKPOINT_END; breakpoint++){
			//CALCULATE AVERAGES BEFORE AND AFTER THE BREAKPOINT
			average1 = 0;
			average2 = 0;
			for(i=0; i<breakpoint; i++){
				average1 += filteredReads[i]; // sum reads from 0 to breakpoint-1
			}
			average1 = average1/((float) breakpoint);
			for(i=breakpoint; i<98; i++){
				average2 += filteredReads[i]; // sum reads from breakpoint to end (99)
			}
			average2 = average2/(98-breakpoint);
			
			//CALCULATE SUM OF SQUARE
			sos=0;
			for(i=0; i<breakpoint; i++){
				sos += pow((filteredReads[i]-average1),2);
			}
			for(i=breakpoint; i<98; i++){
				sos += pow((filteredReads[i]-average2),2);
			}
			sumOfSquares[breakpoint] = sos;
		}
		
		//FIND THE BREAKPOINT WITH MINIMUM SUM OF SQUARES
		minBreakpoint = BREAKPOINT_START;
		for(breakpoint=BREAKPOINT_START; breakpoint<BREAKPOINT_END; breakpoint++){
			if(sumOfSquares[breakpoint] < sumOfSquares[minBreakpoint]){
				minBreakpoint = breakpoint;
			}
		}
		
		//RECALCULATE AVERAGE1 AND AVERAGE2 FOR MINBREAKPOINT
		average1 = 0;
		average2 = 0;
		for(i=0; i<minBreakpoint; i++){
			average1 += filteredReads[i]; // sum reads from 0 to breakpoint-1
		}
		average1 = average1/(minBreakpoint);
		for(i=minBreakpoint; i<98; i++){
			average2 += filteredReads[i]; // sum reads from breakpoint to end (99)
		}
		average2 = average2/(98-minBreakpoint);


		//RECALCULATE SUM OF SQUARES, BUT NORMALIZING FOR AVERAGE
		sos=0;
		for(i=0; i<minBreakpoint; i++){
			sos += pow((filteredReads[i]-average1)/average1,2);
		}
		for(i=minBreakpoint; i<98; i++){
			sos += pow((filteredReads[i]-average2)/average2,2);
		}
		
		//to indicate breakpoint relative to entire data array, need to add +1 or +2 again
		if(minBreakpoint > maxPos1 && minBreakpoint > maxPos1){
			minBreakpoint += 2;
		} else if(minBreakpoint < maxPos1 && minBreakpoint < maxPos1){
			;
		} else {
			minBreakpoint += 1;
		}

		//OUTPUT RESULTS
		fprintf(out,"%s", transcriptName);
		for(i=0; i<100; i++){
			fprintf(out,"\t%f", reads[i]);
		}
		fprintf(out,"\t%i\t%f\t%f\t%f\t%f\n", minBreakpoint, average1, average2, average2-average1, sos);
	}
	

	

	/* CLOSE AND EXIT */

	fclose(in);
	fclose(out);
}

