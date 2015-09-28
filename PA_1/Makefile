all: binarize_sequential binarize_pthreads binarize_openmp sobel_sequential sobel_openmp sobel_pthreads

binarize_sequential: binarize_sequential.o lodepng.o
		gcc -o binarize_sequential binarize_sequential.o lodepng.o

binarize_pthreads: binarize_pthreads.o lodepng.o
		gcc -o binarize_pthreads binarize_pthreads.o lodepng.o

binarize_openmp: binarize_openmp.o lodepng.o
		gcc -fopenmp -o binarize_openmp binarize_openmp.o lodepng.o

binarize_sequential.o: binarize_sequential.c lodepng.h
		gcc -c binarize_sequential.c

binarize_pthreads.o: binarize_pthreads.c lodepng.h
		gcc -c -pthread binarize_pthreads.c

binarize_openmp.o: binarize_openmp.c lodepng.h
		gcc -c binarize_openmp.c

sobel_sequential: sobel_sequential.o lodepng.o
		gcc -o sobel_sequential sobel_sequential.o lodepng.o

sobel_pthreads: sobel_pthreads.o lodepng.o
		gcc -o sobel_pthreads sobel_pthreads.o lodepng.o

sobel_openmp: sobel_openmp.o lodepng.o
		gcc -fopenmp -o sobel_openmp sobel_openmp.o lodepng.o

sobel_sequential.o: sobel_sequential.c lodepng.h
		gcc -c sobel_sequential.c

sobel_pthreads.o: sobel_pthreads.c lodepng.h
		gcc -c -pthread sobel_pthreads.c

sobel_openmp.o: sobel_openmp.c lodepng.h
		gcc -c sobel_openmp.c

lodepng.o: lodepng.c lodepng.h
		gcc -c lodepng.c

clean:
		rm *.o binarize_sequential binarize_openmp binarize_pthreads sobel_sequential sobel_openmp sobel_pthreads 
