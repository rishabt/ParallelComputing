#include "lodepng.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>

/* TODO: it may help to put some global variables here 
for your threads to use */

unsigned part;
unsigned char *image, *new_image;

int partitions;
unsigned total_width;
unsigned total_height;

void *worker_thread(void *arg) {
  /* TODO: put image processing code here */

  int thread_id = (int) arg;

  int start = thread_id * partitions;
  int end = (thread_id + 1) * partitions;

  if(start == 0){
    start += 1;
  }

  if(start >= total_width){
    start -= 1;
  }

  if(end == total_width){
    end -= 1;
  }

  unsigned char value;
  for (int j = 1; j < total_height - 1; j++) {
    for (int i = start; i < end; i++) {

      value = abs((image[4*total_width*(i-1) + 4*(j-1)] + 2 * image[4*total_width*(i-1) + 4*j] +
                image[4*total_width*(i-1) + 4*(j+1)]) - (image[4*total_width*(i+1) + 4*(j-1)] +
                2 * image[4*total_width*(i+1) + 4*j] + image[4*total_width*(i+1) + 4*(j+1)])) +
              abs((image[4*total_width*(i-1) + 4*(j+1)] + 2 * image[4*total_width*i + 4*(j+1)] +
                image[4*total_width*(i+1) + 4*(j+1)]) - (image[4*total_width*(i-1) + 4*(j-1)] +
                2 * image[4*total_width*i + 4*(j-1)] + image[4*total_width*(i+1) + 4*(j-1)]));

      new_image[4*total_width*i + 4*j] = value;
      new_image[4*total_width*i + 4*j + 1] = value;
      new_image[4*total_width*i + 4*j + 2] = value;
      new_image[4*total_width*i + 4*j + 3] = 255;
    }
  }

  pthread_exit(NULL);
}

void sobelize(char* input_filename, char* output_filename, int thread_count)
{
  unsigned error;
  unsigned width, height;

  // load image from PNG into C array
  error = lodepng_decode32_file(&image, &width, &height, input_filename);
  if(error) printf("error %u: %s\n", error, lodepng_error_text(error));
  new_image = malloc(width * height * 4 * sizeof(unsigned char));

  struct timeval start, end; // struct used to compute execution time
  gettimeofday(&start, NULL);  // set starting point

  /* TODO: create your thread team here and send each thread an argument 
  telling it which part of "image" to process 

  remember to join all threads!
  */
  total_width = width;
  total_height = height;

  partitions = width/thread_count;

  part = (width * height)/thread_count;

  int threads[thread_count];

  int i;
  for (i = 0; i < thread_count; i++) {
    pthread_t t;

    threads[i] = pthread_create(&t, NULL, worker_thread, (void*) i);
  }

  int j;
  for (j = 0; j < thread_count; j++){
    pthread_join(threads[i], NULL);
  }


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
