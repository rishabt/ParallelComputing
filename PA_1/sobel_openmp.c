#include "lodepng.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>

void sobelize(char* input_filename, char* output_filename, int thread_count)
{
  unsigned error;
  unsigned char *image, *new_image;
  unsigned width, height;

  error = lodepng_decode32_file(&image, &width, &height, input_filename);
  if(error) printf("error %u: %s\n", error, lodepng_error_text(error));
  new_image = malloc(width * height * 4 * sizeof(unsigned char));
  unsigned char value;

  struct timeval start, end; // struct used to compute execution time
  gettimeofday(&start, NULL);  // set starting point

  /* TODO: put your OpenMP parallel block here */

  gettimeofday(&end, NULL);
  printf("\n\nAlgorithm's computational part duration : %ld\n", \
               ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));


  lodepng_encode32_file(output_filename, new_image, width, height);

  free(image);
  free(new_image);
}

int main(int argc, char *argv[])
{
  char* input_filename = argv[1];
  char* output_filename = argv[2];
  int thread_count = atoi(argv[3]);

  sobelize(input_filename, output_filename, thread_count);
  return 0;
}