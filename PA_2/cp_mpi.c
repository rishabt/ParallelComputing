#include "cp_mpi.h"

int rank;
int np;
int number_rows;
int number_columns;

int buffer_size;
int lines;
int lines_per_buffer; 
int overflow;
int size_of_line;
double* cp;
double* o;

MPI_Status status;

int main(int argc, char *argv[])
{
	char *file_name;
	int charsPerLine;
	
	// Check time
	struct timeval start, end;

	// Verify user input
	if(argc != 2)
	{
		printf("usage: mpirun -np # ./cp_mpi user_matrix.txt\n");
		return -1;
	}

	// MPI initialization
	if(MPI_Init(&argc, &argv) != MPI_SUCCESS)
	{
		printf("MPI initialization error\n");
		exit(-1);	
	}

    // Rank of the calling MPI process within the specified communicator
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Number of MPI processes
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	file_name = argv[1];

	// Get matrix size
	get_matrix_size(file_name);

	// Get line size (in bytes)
	get_size_of_line(&size_of_line);

	// Read the user matrix
	double *user_matrix = user_matrix_read(file_name, size_of_line);

	// Clicking probabilities vector creation
	if((o = malloc((number_rows * number_columns) * sizeof(double))) == NULL || 
			(cp = malloc(lines * sizeof(double))) == NULL)
	{
			printf("Matrix memory allocation error\n");
			exit(-1);
	}

	// Calculate REF and RREF_p
	calculate_REF(user_matrix, lines, number_columns);	
	calculate_RREF_p(user_matrix, lines);
	
	divide_by_max(user_matrix, lines, get_global_maximum(get_local_maximum(user_matrix, lines)));

	MPI_Gather(user_matrix, (number_columns * lines), MPI_DOUBLE, o, 
					number_columns * lines, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	write_file(o, number_rows);
	
	print_best_acceptance_threshold(cp, lines);

	((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));	

	// MPI termination
	if(MPI_SUCCESS != MPI_Finalize())
	{
		printf("Error terminating MPI.\n");
		exit(-1);
	}
	return 0; 
}


int get_row_proc(int row)
{
	int i;
	int a = 0;

	for(i = 0; i < np; i++)
	{
		a += get_limit_of_row(i);

		if(row < a)
		{
			return i;
		}
	}

	return -1;
}

void calculate_REF(double* mat, int lines, int number_columns){
	
	int r, b, a;

	double* temp = malloc(number_columns * lines * sizeof(double));

	if(temp == NULL)
	{
		printf("Malloc error\n");
		exit(-1);
	}	

	for(a = 0; a < rank * lines; a++)
	{
		MPI_Bcast(temp, number_columns, MPI_DOUBLE, (a/lines), MPI_COMM_WORLD);

		for(r = 0; r < lines; r++)
		{
			for(b = a + 1; b < number_columns; b++)
			{
				mat[r * number_columns + b] = mat[r * number_columns + b] - mat[r * number_columns + a] * temp[b];
			}	
			mat[r * number_columns + a] = 0;
		}
	}

	for(r = 0; r < lines; r++)
	{
		for(b = rank * lines + r + 1; b < number_columns; b++)
		{
			mat[r * number_columns + b] = mat[r * number_columns + b]/mat[r * number_columns + rank * lines + r];
		}
		mat[r * number_columns + rank * lines + r] = 1;

		for(a = 0; a < number_columns; a++)
		{
			temp[a] = mat[r*number_columns+a];
		}

		MPI_Bcast(temp, number_columns, MPI_DOUBLE, rank, MPI_COMM_WORLD);

		for(a = r + 1; a < lines; a++)
		{
			for(b = rank * lines + r + 1; b < number_columns; b++)
			{
				mat[a * number_columns + b] = mat[a * number_columns + b] - 
										mat[a * number_columns + r + rank * lines] * temp[b];
			}
			mat[a * number_columns + r + rank * lines] = 0;
		}
	}

	for(a = ((rank + 1) * get_limit_of_row(rank + 1)); a < number_columns - 1; a++)
	{
		MPI_Bcast(temp, number_columns, MPI_DOUBLE, (a/get_limit_of_row(rank + 1)), MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(mat, number_columns * lines, MPI_DOUBLE, o, number_columns * lines, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

int get_limit_of_row(int processor)
{
	int result = get_size_of_buffer(processor)/size_of_line;

	return result;
}

void calculate_RREF(double* mat, int lines)
{
	int r, c, a;
	for(c = number_columns - 2, a = number_rows - 1; c >= 0; c--, a--)
	{
		for(r = 0; r < lines; r++)
		{
			if(c == r)
			{
				break;
			}

			mat[r * number_columns + number_columns - 1] -= mat[r * number_columns + c] * 
													mat[a * number_columns + number_columns - 1];
			mat[r * number_columns + c] = 0;
		}
	}
}

void calculate_RREF_p(double* mat, int lines)
{
	int r, c, a, b;
	double x = 1;

	for(c = number_columns - 2, a = number_rows - 1, b = 1; c >= 0; c--, a--, b++)
	{	
		if(rank == get_row_proc(a))
		{
			if(b%lines == 0)
			{
				MPI_Bcast(&mat[number_columns - 1], 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
				x = mat[number_columns-1];
			}
			else
			{
				MPI_Bcast(&mat[(lines - b%lines) * number_columns + number_columns-1], 
														1, MPI_DOUBLE, rank, MPI_COMM_WORLD);

				x = mat[(lines-b%lines)*number_columns+number_columns-1];
			}	
		}
		else
		{
			MPI_Bcast(&x, 1, MPI_DOUBLE, get_row_proc(a), MPI_COMM_WORLD);
		}

		for(r = 0;r < lines; r++)
		{
			if(c == r || mat[r*number_columns+c] == 1)
			{
				break;
			}

			mat[r * number_columns + number_columns - 1] -= mat[r * number_columns + c] * x;

			mat[r * number_columns + c] = 0;
		}
	}
	
}

void write_file(double* mat, int lines)
{
	if(rank == 0){
	FILE *file;

	if((file = fopen("clicking_probabilities.txt", "w")) == 0)
	{
		printf("File creation error\n");
		exit(-1);
	}


	int a;
	for(a = 0; a < lines; a++)
	{
		fprintf(file, "%lf\n", *(mat + a * number_columns + number_columns - 1));
	}

	fclose(file);
	}
}

void divide_by_max(double *mat, int lines, double maximum)
{
	int a;
	for(a = 0; a < lines; a++)
	{
		mat[a * number_columns + number_columns - 1] = fabs(mat[a * number_columns + number_columns - 1]/maximum);

		*(cp+a) = *(mat+a*number_columns+number_columns-1);
	}
}

double get_global_maximum(double local_maximum)
{
	double global_maximum = local_maximum;

	if(np == 1)			// If number of processes is 1
	{
		return global_maximum;
	}
	
	if(rank != 0)		// Send to root
	{
		MPI_Send(&local_maximum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MPI_Bcast(&global_maximum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	else
	{
		int a;
		double *max = malloc(sizeof(double) * np - 1);

		for(a = 1; a < np; a++)
		{
			MPI_Recv((max + a - 1), 1, MPI_DOUBLE, a, 0, MPI_COMM_WORLD, &status);
		}
	
		for(a = 0; a < np - 1; a++)
		{
			if(global_maximum < max[a])
			{
				global_maximum = max[a];
			}
		}

		MPI_Bcast(&global_maximum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	return global_maximum;
}

double get_local_maximum(double* mat, int lines)
{	
	int a;
	double local_maximum = 0;

	for(a = 0; a < lines; a++)
	{
		if(fabs(mat[a * number_columns + number_columns - 1]) > local_maximum)
		{
			local_maximum = fabs(mat[a * number_columns + number_columns - 1]);
		}
	}

	return local_maximum;
}


double* aloc_matrix(int lines)
{
	double* mat = malloc((lines*number_columns) * sizeof(double));

	if(mat == NULL)
	{
		printf("Matrix allocation error\n");
		exit(-1);
	}

	return mat;
}

void print(double *mat)
{
	if(rank == 0)
	{
		int a, b;

		for(a = 0; a < number_rows; a++)
		{
			for(b = 0; b < number_columns; b++)
			{
				printf("%lf ", mat[a * number_columns + b]);
			}

			printf("\n");
		}
	}
}

void print_matrix(double *mat, int i, int j) 
{
	if(rank == 0)
	{
		int r, c;
		for (r = 0; r < i; r++) 
		{
	 		for (c = 0; c < j; c++) 
	 		{
	 			printf("%lf ",mat[r * number_columns + c]);
		 	}

		 	printf("\n");
		}

		int pr;
		for(pr = 1; pr < np; pr++)
		{
			double *buf = malloc(9 * 4 * sizeof(double));

			if(buf == NULL)
			{
				printf("Malloc error\n");
				exit(-1);
			}

			MPI_Recv(buf, 9 * 4, MPI_DOUBLE, pr, 0, MPI_COMM_WORLD, &status);

			for(r = 0; r < i; r++)
			{
				for (c = 0; c < j; c++) 
				{
					printf("%lf ", buf[r * number_columns + c]);
				}

				printf("\n");
			}
		}
	}
	else
	{
		MPI_Send(mat, 9 * 4, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
}

double* user_matrix_load(char* buf, int size_of_line)
{
	lines = buffer_size/size_of_line;
	double* mat = aloc_matrix(lines);
	
	char tmp[8];

	int a = 0;
	int b = 0;
	int k = 0;

	while(*buf != '\0')
	{
		if(*buf == '\n')
		{
			a++;
			b = 0;
		}
		else if(*buf == ' ')
		{
			b++;
		}
		else
		{
			tmp[k++] = *buf;
		}

		if(k == 8)
		{
			k = 0;
			mat[a * number_columns + b] = atof(tmp);
		}

		*buf++;
	}

	return mat;
}


int get_size_of_buffer(int pr)
{
	int result;

	if(pr < overflow)
	{
		result = (lines_per_buffer+1) * size_of_line;
	}
	else
	{
		result = lines_per_buffer * size_of_line;
	}

	return result;
}

double * user_matrix_read(char* file_name, int size_of_line)
{
	MPI_File mp_file;

	lines_per_buffer = number_rows/ np;
	overflow = number_rows % np;

	buffer_size = get_size_of_buffer(rank);

	char *buf;
	if((buf = malloc(buffer_size)) == NULL){
		printf("Buffer creation error\n");
	}

	int s;

	if(rank == 0)
	{
		s = 0;
	}
	else if(overflow == 0)
	{
		s = ((lines_per_buffer) * size_of_line * rank);
	}
	else if(rank <= overflow)
	{
		s = ((lines_per_buffer+1)*size_of_line * rank);
	}
	else
	{
		s = ((lines_per_buffer+1)*size_of_line * (rank-1)) + 
				((lines_per_buffer)*size_of_line * (rank - overflow));
	}

	MPI_File_open(MPI_COMM_WORLD, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &mp_file);
	MPI_File_seek(mp_file, s, MPI_SEEK_SET);
	MPI_File_read(mp_file, buf, buffer_size, MPI_CHAR, &status);


	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&mp_file);

	return user_matrix_load(buf, size_of_line);
}

void get_size_of_line(int* size_of_line)
{
	*size_of_line = sizeof(double) * number_columns + 
						sizeof(char) * (number_columns-1) + sizeof(char);
}

void get_matrix_size(char *file_name)
{
	FILE *f;

	if((f = fopen(file_name, "r")) == 0)
	{
		printf("Couldn't open file:  %s.\n", file_name);
		exit(-1);
	}

	number_rows = 1;
	number_columns = 1;

	char c;
	int known_columns = 0;

	while(!feof(f))
	{
		c = fgetc(f);

		if(c == ' ')
		{
			if(!known_columns)
				(number_columns)++;
		}

		if(c == '\n')
		{
			(number_rows)++;
			known_columns = 1;
			continue;
		}
	}
	fclose(f);
}

void print_best_acceptance_threshold(double* cp, int lines) 
{	
	double threshold; 
	double l_profit;
	double g_profit;
	double profit_max = 0;
	double i;

	int a;

	for(threshold = 0.2; threshold <= 1.0; threshold += 0.2)
	{
		l_profit = 0;

		for(a = 0; a < lines; a++)
		{
			if(*(cp+a) > threshold)
			{
				l_profit += *(cp+a) - 2*(1-*(cp+a));
			}
		}

		MPI_Reduce(&l_profit, &g_profit, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		if(rank == 0)
		{
			if(g_profit > profit_max)
			{
				profit_max = g_profit;
				i = threshold;
			}
		}
	}

	if(rank == 0)
	{
		printf("Acceptance threshold to max profit: %lf -> $%lf\n", i, profit_max);
	}

}

