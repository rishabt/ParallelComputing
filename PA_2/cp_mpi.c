#include "cp_mpi.h"

int rank, np;
int rows, columns;

int buffer_size, lines, lines_per_buffer, overflow, size_of_line;
double *cp, *m;

MPI_Status status;

int main(int argc, char *argv[]){
	char *fileName;
	int charsPerLine;
	
	struct timeval start, end; //struct used to compute execution time

	//check user input
	if(argc != 2){
		printf("usage: mpirun -np # ./cp_mpi matrix.txt.\n");
		return -1;
	}

	//initialize MPI
	if(MPI_SUCCESS != MPI_Init(&argc, &argv)){
		printf("error initializing MPI.\n");
		exit(-1);	
	}
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	fileName = argv[1];
/*	int z;
	for(z=0;z<100;z++){
	if(rank == 0){
		gettimeofday(&start, NULL);
	}*/

	//get size of matrix
	get_matrix_size(fileName);
	//get size of each line (in bytes)
	get_size_of_line(&size_of_line);
	//read user matrix
		/*
			note, to ensure that the matrix will always be located in a contiguous
			portion of memory, only a single dimensional array will be used to hold the 
			matrix. 
		*/
	double *matrix = user_matrix_read(fileName, size_of_line);
//	print_matrix(matrix, lines, columns);
	//create clicking probabilities vector

	if((m = malloc((rows * columns) * sizeof(double))) == NULL || (cp = malloc(lines * sizeof(double))) == NULL){
			printf("error allocating matrix.\n");
			exit(-1);
	}	
	//computation//
	calculate_REF(matrix, lines, columns);	
	calculate_RREF_p(matrix, lines);
	
//	print_matrix(matrix, lines, columns);
	divide_by_max(matrix, lines, find_global_max(find_local_max(matrix, lines)));

	MPI_Gather(matrix,columns*lines,MPI_DOUBLE,m,columns*lines,MPI_DOUBLE,0,MPI_COMM_WORLD);

	write_clicking_probabilities(m, rows);
	print_best_acceptance_threshold(cp, lines);
//	if(rank == 0){
//		gettimeofday(&end, NULL);
//		printf(/*"\n\nAlgorithm's computational part duration :*/"%ld\n", \
					((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));
//	}	

//	}
	//terminate MPI
	if(MPI_SUCCESS != MPI_Finalize()){
		printf("error terminating MPI.\n");
		exit(-1);
	}
	return 0; 
}


int get_proc_with_row(int row){
	int i, acc=0;
	for(i=0;i<np;i++){
		acc += get_limit_of_row(i);
		if(row < acc){
			return i;
		}
	}
	return -1;
}

void calculate_REF(double *matrix, int row_limit, int columns){
	
	int row, j, i;
	double *tmp;
	if((tmp = malloc(columns*row_limit*sizeof(double))) == NULL){
		printf("error mallocing.\n");
		exit(-1);
	}	

	for(i=0;i<rank*row_limit;i++){
		MPI_Bcast(tmp, columns, MPI_DOUBLE, i/row_limit, MPI_COMM_WORLD);
		//adjust rows based on tmp
		for(row=0;row<row_limit;row++){
			for(j=i+1;j<columns;j++){
				matrix[row*columns + j] = matrix[row*columns + j] - matrix[row*columns + i] * tmp[j];
			}	
			matrix[row*columns + i] = 0;
		}
	}

	for(row=0;row<row_limit;row++){
		for(j=rank*row_limit+row+1;j<columns;j++){
			matrix[row*columns+j] = matrix[row*columns+j] / matrix[row*columns+rank*row_limit+row];
		}
		matrix[row*columns+rank*row_limit+row] = 1;

		for(i=0;i<columns;i++){
			tmp[i] = matrix[row*columns+i];
		}
		MPI_Bcast(tmp, columns, MPI_DOUBLE, rank, MPI_COMM_WORLD);

		for(i=row+1;i<row_limit;i++){
			for(j=rank*row_limit+row+1;j<columns;j++){
				matrix[i*columns+j] = matrix[i*columns+j] - matrix[i*columns+row+rank*row_limit]*tmp[j];
			}
			matrix[i*columns+row+rank*row_limit] = 0;
		}
	}
	//TODO: fix
	for(i=(rank+1)*get_limit_of_row(rank+1);i<columns-1;i++){
		MPI_Bcast(tmp, columns, MPI_DOUBLE, i/get_limit_of_row(rank+1), MPI_COMM_WORLD);
	}
/*	for(i=(rank+1)*row_limit;i<columns-1;i++){
		printf("%d\n", i/row_limit);
		MPI_Bcast(tmp, columns, MPI_DOUBLE, i/row_limit, MPI_COMM_WORLD);
	}*/

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(matrix,columns*row_limit,MPI_DOUBLE,m,columns*row_limit,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

//returns the number of rows in the inputed processor
int get_limit_of_row(int proc){
	return get_size_of_buffer(proc)/size_of_line;
}

void calculate_RREF(double *matrix, int row_limit){
	int row, col, i;
	for(col=columns-2, i=rows-1;col>=0;col--, i--){
		for(row=0;row<row_limit;row++){
			if(col == row)
				break;
			matrix[row*columns+columns-1] -= matrix[row*columns+col]*matrix[i*columns+columns-1];
			matrix[row*columns+col] = 0;
		}
	}
}

void calculate_RREF_p(double *matrix, int row_limit){
	int row, col, i, j;
	double b=1;
	for(col=columns-2, i=rows-1, j=1;col>=0;col--, i--, j++){
		
		if(rank == get_proc_with_row(i)){
			if(j%row_limit == 0){
				MPI_Bcast(&matrix[columns-1], 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
				b = matrix[columns-1];
			}else{
				MPI_Bcast(&matrix[(row_limit-j%row_limit)*columns+columns-1], 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
				b = matrix[(row_limit-j%row_limit)*columns+columns-1];
			}	
		}else{
			MPI_Bcast(&b, 1, MPI_DOUBLE, get_proc_with_row(i), MPI_COMM_WORLD);
		}

		for(row=0;row<row_limit;row++){
			if(col == row || matrix[row*columns+col] == 1)
				break;
			matrix[row*columns+columns-1] -= matrix[row*columns+col]*b;
			matrix[row*columns+col] = 0;
		}
	}
	
}

void write_clicking_probabilities(double *matrix, int row_limit){
	if(rank == 0){
	FILE *file;

	if((file = fopen("clicking_probabilities.txt", "w")) == 0){
		printf("error creating file.\n");
		exit(-1);
	}


	int i;
	for(i=0;i<row_limit;i++){
		fprintf(file, "%lf\n", *(matrix+i*columns+columns-1));
	}

	fclose(file);
	}
/*	MPI_File file;

	char *fileName = "clicking_probabilities.txt";
	MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_CREATE, MPI_INFO_NULL, &file);

	int i;
	for(i=0;i<row_limit;i++){
		printf("%lf\n", *(matrix+i*columns+columns-1));
		MPI_File_write_shared(file, matrix+(i*columns+columns-1), 1, MPI_DOUBLE, &status);
	}

	MPI_Barrier(MPI_COMM_WORLD); //wait until all processes finished reading
	MPI_File_close(&file);*/

}

void divide_by_max(double *matrix, int row_limit, double max){
	int i;
	for(i=0;i<row_limit;i++){
		matrix[i*columns+columns-1] = fabs(matrix[i*columns+columns-1] / max);
		*(cp+i) = *(matrix+i*columns+columns-1);
	}
}

double find_global_max(double localMax){
	double globalMax = localMax;
	if(np == 1){
		return globalMax;
	}
	//send all local maxi to the root to process
	if(rank != 0){
		MPI_Send(&localMax, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MPI_Bcast(&globalMax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}else{ //root

		int i;
		double *maxi;
		maxi = malloc(sizeof(double) * np-1);
		for(i=1;i<np;i++){
			MPI_Recv((maxi+i-1), 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
		}
	
		for(i=0;i<np-1;i++){
			if(globalMax < maxi[i]){
				globalMax = maxi[i];
			}
		}
		MPI_Bcast(&globalMax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	return globalMax;
}

double find_local_max(double *matrix, int row_limit){
	
	int i;
	double localMax = 0;
	for(i=0;i<row_limit;i++){
		if(fabs(matrix[i*columns+columns-1]) > localMax)
			localMax = fabs(matrix[i*columns+columns-1]);
	}
	return localMax;
}


double * allocate_matrix(int lines){
	double *matrix;
	if((matrix = malloc((lines*columns) * sizeof(double))) == NULL){
			printf("error allocating matrix.\n");
			exit(-1);
	}	
	return matrix;
}

void print(double *matrix){
	if(rank == 0){
	int i, j;
	for(i=0;i<rows;i++){
		for(j=0;j<columns;j++){
			printf("%lf ", matrix[i*columns+j]);
		}
		printf("\n");
	}
	}
}

void print_matrix(double *matrix, int i, int j) 
{
	if(rank == 0){
		int row, column;
		for (row = 0; row < i; row++) {
	 		for (column = 0; column < j; column++) {
	 			printf("%lf ",matrix[row*columns+column]);
		 	}
		 	printf("\n");
		}	

		//print from other process
		int p;
		for(p=1;p<np;p++){
			double *buffer = malloc(9*4*sizeof(double));
			if(buffer == NULL){
				printf("error mallocing.\n");
				exit(-1);
			}
			MPI_Recv(buffer, 9*4, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &status);
			for(row =0;row<i;row++){
				for (column = 0; column < j; column++) {
					printf("%lf ", buffer[row*columns + column]);
				}
				printf("\n");
			}
		}
	}else{
		MPI_Send(matrix, 9*4, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
}

//loads contents in the buffer into a matrix
double * user_matrix_load(char *buffer, int size_of_line){
	lines = buffer_size / size_of_line;
	
	double *matrix = allocate_matrix(lines);
	
	char temp[8];
	int i, j, k;
	j = i = k = 0;
	while(*buffer != '\0'){
		if(*buffer == '\n'){
			i++;
			j=0;
		}else if(*buffer == ' '){
			j++;
		}else{
			temp[k++] = *buffer;
		}
		if(k == 8){
			k = 0;
			matrix[i*columns+j] = atof(temp);
		}
		*buffer++;
	}
	return matrix;
}

int get_size_of_buffer(int proc){
	return (proc < overflow)? (lines_per_buffer+1) * size_of_line : lines_per_buffer * size_of_line; 
}

//reads matrix from file and stores into buffers
double * user_matrix_read(char *fileName, int size_of_line){
	MPI_File file;

	//every process will get a buffer to read its section
	//lets divide up our file..

	//how many lines per buffer?
	lines_per_buffer = rows / np;
	//don't forget the overflow!
	overflow = rows % np; //overflow: [0- np-1]

	//so what's the size of each buffer? 
		//overflow check
	buffer_size = get_size_of_buffer(rank);//(rank < overflow)? (lines_per_buffer+1) * size_of_line : lines_per_buffer * size_of_line; 

	//lets create the buffer
	char *buffer;
	if((buffer = malloc(buffer_size)) == NULL){
		printf("error creating buffer.\n");
	}

	//lastly, each process needs to know where to seek the file pointer
		//because of the overflow, it's not just a simple multiplication
	int seek;
	if(rank == 0){
		seek = 0;
	}else if(overflow == 0){
		seek = ((lines_per_buffer)*size_of_line * rank);
	}else if(rank <= overflow){
		seek = ((lines_per_buffer+1)*size_of_line * rank);
	}else{
		seek = ((lines_per_buffer+1)*size_of_line * (rank-1)) + ((lines_per_buffer)*size_of_line * (rank-overflow));
	}

	MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
	MPI_File_seek(file, seek, MPI_SEEK_SET);
	MPI_File_read(file, buffer, buffer_size, MPI_CHAR, &status);


	MPI_Barrier(MPI_COMM_WORLD); //wait until all processes finished reading
	MPI_File_close(&file);

	return user_matrix_load(buffer, size_of_line);
}

//determine size of each line in bytes
void get_size_of_line(int *size_of_line){
	//			each element in the row   +   the spaces       + the carriage return					 
	*size_of_line = sizeof(double) * columns + sizeof(char) * (columns-1) + sizeof(char);
}

//determine size of matrix
void get_matrix_size(char *fileName){
	FILE *file;
	//open & check if valid file
	if((file = fopen(fileName, "r")) == 0){
		printf("error openening file %s.\n", fileName);
		exit(-1);
	}

	rows = 1;
	columns = 1;
	char c;
	int columns_known = 0;
	while(!feof(file)){
		c = fgetc(file);
		if(c == ' '){
			if(!columns_known)
				(columns)++;
		}
		if(c == '\n'){
			(rows)++;
			columns_known = 1;
			continue;
		}
	}
	fclose(file);
}

void print_best_acceptance_threshold(double *cp, int row_limit) {
	
	double threshold, local_profit;
	double global_profit, max_profit = 0, t;
	int i;

	for(threshold=0.2;threshold<=1.0;threshold+=0.2){
		local_profit = 0;
		for(i=0;i<row_limit;i++){
			if(*(cp+i) > threshold){
				local_profit += *(cp+i) - 2*(1-*(cp+i));
			}
		}
		MPI_Reduce(&local_profit, &global_profit, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if(rank == 0){
			if(global_profit > max_profit){
				max_profit = global_profit;
				t = threshold;
			}
		}
	}

	if(rank == 0)
		printf("Acceptance threshold to max profit: %lf -> $%lf\n", t, max_profit);

}

