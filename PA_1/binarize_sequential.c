#include "lodepng.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#define THRESHOLD 200

void binarize(char* input_filename, char* output_filename)
{
  unsigned error;
  unsigned char *image, *new_image;
  unsigned width, height;

  error = lodepng_decode32_file(&image, &width, &height, input_filename);
  if(error) printf("error %u: %s\n", error, lodepng_error_text(error));
  new_image = malloc(width * height * 4 * sizeof(unsigned char));

  struct timeval start, end; // struct used to compute execution time
  gettimeofday(&start, NULL);  // set starting point

  /* TODO: Loop through the "image" and generate the "new_image". */
  unsigned char value;
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      /* TODO use thresholding on image[4*width*i + 4*j], put
      the result into the "value" variable (which then gets put into
      the new image's pixel)
      */

      if(image[4*width*i + 4*j] > THRESHOLD)
      {
          value = 255;
      }
      else
      {
          value = 0;
      }

      new_image[4*width*i + 4*j] = value;
      new_image[4*width*i + 4*j + 1] = value;
      new_image[4*width*i + 4*j + 2] = value;
      new_image[4*width*i + 4*j + 3] = 255;
    }
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

  binarize(input_filename, output_filename);

  return 0;
}
