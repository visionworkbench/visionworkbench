#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include <vw/FileIO/JP2.h>

using namespace vw;


int main(int argc, char** argv)
{
  uint8* d;
  FILE* fp;
  uint64 nbytes;
  uint64 i;
  int c;

  if(argc < 2)
  {
    fprintf(stderr, "USAGE: %s jp2file\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  if(!(fp = fopen(argv[1], "r")))
  {
    fprintf(stderr, "ERROR: Unable to open file '%s' for reading\n", argv[1]);
    exit(EXIT_FAILURE);
  }

  for(nbytes = 0; fgetc(fp) != EOF; nbytes++);
  rewind(fp);

  d = new uint8[nbytes];

  for(i = 0; i < nbytes && (c = fgetc(fp)) != EOF; d[i] = (uint8)c, i++);

  if(!(i == nbytes && fgetc(fp) == EOF))
  {
    fprintf(stderr, "ERROR: File '%s' has changed size\n", argv[1]);
    exit(EXIT_FAILURE);
  }

  fclose(fp);

  JP2File f(d, nbytes);
  f.print();

  delete[] d;

  return 0;
}
