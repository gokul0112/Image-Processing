/********************************************************
***IMPORTANT NOTE***
If you have problems with the provided sample code,
part of the reason might be due to the function "fopen".
Please try changing "r/w" to "rb/wb" or the other way
when you use this function.
*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <malloc.h>  
#include <memory.h>

#define max(x, y) ((x>y) ? (x):(y))
#define min(x, y) ((x<y) ? (x):(y))

#define SIZE 10
#define SIGMA 1
#define SIGMA_SPACE 25.0
#define SIGMA_RANGE 25.0

int xdim;
int ydim;
int maxraw;
double sum_of_values=0;
double normalize=0;
unsigned char *image;
unsigned char *image2, *image3;
void ReadPGM(FILE*);
double gaussian(int x, int y);
void gaussian_kernel(double kernel[SIZE][SIZE]);
void bilateral_filter(int xdim, int ydim, FILE* fp, char **argv);
void gaussian_filter(FILE* fp,double kernel[SIZE][SIZE],char **argv);
void WritePGM(FILE*,int xdim, int ydim);
void WritePGM2(FILE* fp,int xdim,int ydim);
void kernel2(double values[SIZE*SIZE],double kernel[SIZE][SIZE],int median);
int i, j;

int main(int argc, char **argv)
{
  char str[20];
  double kernel[SIZE][SIZE];
  FILE *fp;
  printf(" value of argc: %d", argc);
  if (argc != 4){
    printf("Usage: MyProgram <input_ppm> <output_ppm> <output2_ppm>\n");
    printf("       <input_ppm>: PGM file \n");
    printf("       <output_ppm>: PGM file \n");
    printf("       <output2_ppm>: PGM file \n");
    exit(0);              
  }

  /* begin reading PGM.... */
  printf("begin reading PGM.... \n");
  if ((fp=fopen(argv[1], "rb"))==NULL){
    printf("read error...\n");
    exit(0);
  }
  ReadPGM(fp);
  gaussian_kernel(kernel);
  printf("Applying gaussian filter...\n");
  gaussian_filter(fp, kernel,argv);
  printf("Applying bilateral filter...\n");
  bilateral_filter(xdim, ydim,fp,argv);
  free(image);
  free(image2);
  return (1);
}
void gaussian_filter(FILE* fp,double kernel[SIZE][SIZE],char **argv) {
    int i, j, k, l;
    int sum;
    image2 = (unsigned char*)malloc(sizeof(unsigned char)*xdim*ydim);
    // Iterate over each pixel in the image
    for(j = 0; j < ydim; j++) {//height
        for(i = 0; i < xdim; i++) {//width
            sum = 0;
            // Apply the kernel to the neighborhood of the current pixel
            for(k = 0; k < SIZE; k++) {
              for(l = 0; l < SIZE; l++) {
                    sum += kernel[k][l] * image[(j+k)* xdim +(i+l)];
                }
            }
            image2[j * xdim + i] = sum;
        }
    }
    printf("Begin writing PGM.... \n");
    if ((fp=fopen(argv[2], "wb")) == NULL){
     	printf("write pgm error....\n");
     	exit(0);
    }
    WritePGM(fp,xdim,ydim);
}
double gaussian(int x, int y) { //used for gaussian filter
    return 1 / (2 * M_PI * SIGMA * SIGMA) * exp((-(x * x + y * y) / (2 * SIGMA * SIGMA)));
}
double gaussian2(double x, double sigma) { // used for bilateral filter
    return exp((-(x*x) / (2 * sigma * sigma)));
}
void gaussian_kernel(double kernel[SIZE][SIZE]) { // global mask
    int center = SIZE / 2;
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            kernel[i][j] = gaussian((i - center),(j - center));
            sum_of_values+=kernel[i][j];
        }
    }
    //printf("sum_of_values ->%f--\n",sum_of_values);
    for (int i = 0; i < SIZE; i++) {		//normalizing
        for (int j = 0; j < SIZE; j++) {
            kernel[i][j] = kernel[i][j]/sum_of_values;
            normalize+=kernel[i][j];
        }
    }
    //printf("normalize ->%f--\n",normalize); //uncomment to check if values sum to 1.
}

void bilateral_filter(int xdim, int ydim, FILE* fp, char **argv) {
    int x, y, i, j,k,l;
    double intensity_difference, spatial_weight, intensity_weight, total_weight;
    double gaussian_space, gaussian_intensity;
    double values[SIZE*SIZE];
    image3 = (unsigned char*)malloc(sizeof(unsigned char)*xdim*ydim);
    for (x = 0; x < ydim; x++) {
        for (y = 0; y < xdim; y++) {
        double kernel[SIZE][SIZE];
            total_weight = 0;
            for (i = 0; i < SIZE; i++) {
                for (j = 0; j < SIZE; j++) {
                    
                    {
                        values[i*SIZE+j] = image[(x+i)*xdim+(y+j)];
                
                    }
                }
            }
            
            int center = SIZE / 2;
            sum_of_values = 0;
            normalize = 0;
            for (int i = 0; i < SIZE; i++) { 				//creating mask for each pixel. 
            	for (int j = 0; j < SIZE; j++) {
        		int intensity_val = values[i*SIZE+j] - values[center];
            		kernel[i][j] = gaussian2((i-center),SIGMA_SPACE)*gaussian2(intensity_val,SIGMA_RANGE);
            		sum_of_values+=kernel[i][j];
            	}
            }
            for (int i = 0; i < SIZE; i++) {                        //normalizing.
        	for (int j = 0; j < SIZE; j++) {
            		kernel[i][j] = kernel[i][j]/sum_of_values;
            		normalize+=kernel[i][j];
        	}
    	    } 
    	    //printf("normalize ->%f--\n",normalize); //uncomment to check if values sum to 1.
            double sum=0;
            for(k = 0; k < SIZE; k++) {
            	for(l = 0; l < SIZE; l++) {
                    	sum += kernel[k][l] * image[(x+k)* xdim +(y+l)];
                }
            }
            image3[x * xdim + y] = sum;
        }
    }
    printf("Begin writing PGM.... \n");
    if ((fp=fopen(argv[3], "wb")) == NULL){
     	printf("write pgm error....\n");
     	exit(0);
    }
    WritePGM2(fp,xdim,ydim); 
}
void ReadPGM(FILE* fp)
{
    int c;
    int i,j;
    int val;
    unsigned char *line;
    char buf[1024];


    while ((c=fgetc(fp)) == '#')
        fgets(buf, 1024, fp);
     ungetc(c, fp);
     if (fscanf(fp, "P%d\n", &c) != 1) {
       printf ("read error ....");
       exit(0);
     }
     if (c != 5 && c != 2) {
       printf ("read error ....");
       exit(0);
     }

     if (c==5) {
       while ((c=fgetc(fp)) == '#')
         fgets(buf, 1024, fp);
       ungetc(c, fp);
       if (fscanf(fp, "%d%d%d",&xdim, &ydim, &maxraw) != 3) {
         printf("failed to read width/height/max\n");
         exit(0);
       }
       printf("Width=%d, Height=%d \nMaximum=%d\n",xdim,ydim,maxraw);

       image = (unsigned char*)malloc(sizeof(unsigned char)*xdim*ydim);
       getc(fp);

       line = (unsigned char *)malloc(sizeof(unsigned char)*xdim);
       for (j=0; j<ydim; j++) {
          fread(line, 1, xdim, fp);
          for (i=0; i<xdim; i++) {
            image[j*xdim+i] = line[i];
         }
       }
       free(line);

     }

     else if (c==2) {
       while ((c=fgetc(fp)) == '#')
         fgets(buf, 1024, fp);
       ungetc(c, fp);
       if (fscanf(fp, "%d%d%d", &xdim, &ydim, &maxraw) != 3) {
         printf("failed to read width/height/max\n");
         exit(0);
       }
       printf("Width=%d, Height=%d \nMaximum=%d,\n",xdim,ydim,maxraw);

       image = (unsigned char*)malloc(sizeof(unsigned char)*xdim*ydim);
       getc(fp);

       for (j=0; j<ydim; j++)
         for (i=0; i<xdim; i++) {
            fscanf(fp, "%d", &val);
            image[j*xdim+i] = val;
         }

     }

     fclose(fp);
}


void WritePGM(FILE* fp,int xdim,int ydim)
{
  int i,j;
  fprintf(fp, "P5\n%d %d\n%d\n", xdim, ydim, 255);
  for (j=0; j<ydim; j++)
    for (i=0; i<xdim; i++) {
      fputc(image2[j*xdim+i], fp);
    }

  fclose(fp);
  
}
void WritePGM2(FILE* fp,int xdim,int ydim)
{
  int i,j;
  fprintf(fp, "P5\n%d %d\n%d\n", xdim, ydim, 255);
  for (j=0; j<ydim; j++)
    for (i=0; i<xdim; i++) {
      fputc(image3[j*xdim+i], fp);
    }

  fclose(fp);
  
}
