#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "input.h"
#include "util.h"
#include <math.h>
#include <sys/time.h>

parameterBlock::parameterBlock( void )
{
        defaults_set = 0;
}

void printParamBlock( parameterBlock *block );

void setDefaults( parameterBlock *block )
{
  block->defaults_set = 1;

  block->do_planar = 0;

  block->T = 0.592;
  block->r = 10;
  block->kc = 20;
  block->c0 = 0.0;
  block->chi = 0.0;
  block->eta = 0.128;
  block->Ap = 0.0;
  block->D0 = 0.0;

  if(!block->do_planar){
    block->mode_max = 60;
  }
  if(block->do_planar){
    block->mode_max = 4.0;
  }

  block->rate = 0.01;
  block->Nsteps = 7000;
}


#define FILE_PASS               0
#define COMMAND_LINE_PASS       1
#define END_PASS                2

void getInput( const char **argv, int argc, parameterBlock *block)
{
  printf("Input line:\n");
  for( int x = 0; x < argc; x++ )
    {
      printf("%s ", argv[x] );
    }
  printf("\n");

  char *word1 = (char *)malloc( sizeof(char) * 4096 );
  char *word2 = (char *)malloc( sizeof(char) * 4096 );

  int is_file = 1;

  if( argc > 1 )
    {
      for( int t = 0; t < strlen(argv[1]); t++ )
	{
	  if( argv[1][t] == '=' )
	    is_file = 0;
	}
    }

  FILE *theFile = NULL;
  int pass = FILE_PASS;

  int c = 2;
  if( argc <= 1 )
    {
      printf("Using defaults.\n");
      pass = END_PASS;
    }
  else if( !is_file )
    {
      printf("Interpreting the first argument as setting a parameter, not an input file.\n");
      pass = COMMAND_LINE_PASS;
      c = 1;
    }
  else
    {
      theFile = fopen(argv[1],"r");
      if( !theFile )
	{
	  printf("Couldn't open input file '%s'\n", argv[1] );
          exit(1);
        }
      printf("Input file:\n");
      printf("#####################################\n");
    }
  char *buffer = (char *)malloc( sizeof(char) * 100000 );
  if( !block->defaults_set )
    setDefaults(block);

  int ERROR = 0;

  while( pass != END_PASS )
    {
      const char *tbuf;
      if( pass == FILE_PASS )
	{
	  getLine( theFile, buffer );
	  if( feof(theFile) )
	    {
	      printf("#####################################\n");
              pass = COMMAND_LINE_PASS;
              continue;
	    }
	  printf("%s\n", buffer );
          if( strlen(buffer) > 4095 )
	    {
	      printf("Line %s too long.\n", buffer );
              ERROR = 1;
	    }
	  if( buffer[0] == '\0' )
	    continue;
	  tbuf = buffer;
          const char *p = tbuf;
          while( *p == '\t' || *p == ' ' ) p += 1;

          if( *p == '#' )
	    {
	      continue;
	    }

	  int nr = sscanf( buffer, "%s %s", word1, word2 );

	  if( nr != 2 )
	    {
	      printf("Could not interpret input line '%s'.\n", buffer );
              ERROR = 1;
	    }
	}
      else if( pass == COMMAND_LINE_PASS )
	{
	  if( c >= argc )
	    {
	      pass = END_PASS;
              continue;
	    }
	  word1[0] = '\0';
          word2[0] = '\0';
          int x;
          for( x = 0; x < strlen(argv[c]); x++ )
	    {
	      if( argv[c][x] == '=' )
		break;
	      word1[x] = argv[c][x];
              word1[x+1] = '\0';
	    }
	  x++;
          for( int y = 0; x < strlen(argv[c]); x++, y++ )
	    {
	      word2[y] = argv[c][x];
              word2[y+1] = '\0';
	    }

	  if( strlen(word1) < 1 || strlen(word2) < 1 )
	    {
	      printf("Couldn't parse command line argument '%s'.\n", argv[c] );
	      ERROR = 1;
	    }

	  tbuf = argv[c];
          c++;
	}
      else break;

      if( !strcasecmp( word1, "do_planar" ) )
	{
          if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
            block->do_planar = 1;
          else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
            block->do_planar = 0;
          else
            {
              printf("Could not interpret input line '%s'.\n", tbuf );
	      ERROR = 1;
	    }
        }

      else if( !strcasecmp( word1, "T" ) )
	block->T = atof(word2);
      else if( !strcasecmp( word1, "r" ) )
	block->r = atof(word2);
      else if( !strcasecmp( word1, "kc" ) )
	block->kc = atof(word2);
      else if( !strcasecmp( word1, "chi" ) )
        block->chi = atof(word2);
      else if( !strcasecmp( word1, "c0" ) )
        block->c0 = atof(word2);
      else if( !strcasecmp( word1, "eta" ) )
        block->eta = atof(word2);
      else if( !strcasecmp( word1, "Ap" ) )
        block->Ap = atof(word2);
      else if( !strcasecmp( word1, "D0" ) )
        block->D0 = atof(word2);
      else if( !strcasecmp( word1, "mode_max" ) )
        block->mode_max = atof(word2);
      else if( !strcasecmp( word1, "rate" ) )
        block->rate = atof(word2);
      else if( !strcasecmp( word1, "Nsteps" ) )
        block->Nsteps = atof(word2);
      else
	{
	  printf("Could not interpret input line '%s'.\n", tbuf );
          ERROR = 1;
	}
    }
  if ( ERROR )
    {
      exit(1);
    }

  printParamBlock(block);
  
  free(word1);
  free(word2);
  free(buffer);

}

void printParamBlock( parameterBlock *block )
{
  printf("Parameters:\n");
  printf("/****************************/\n");
  if( !block->do_planar )
    {
      printf("Simulation Type:      Spherical\n");
    }
  if( block->do_planar )
    {
      printf("Simulation Type:      Planar\n");
    }
  printf("T:      %f kcal/mol\n", block->T);
  printf("radius:      %f microns\n", block->r);
  printf("Bending Modulus:      %f kcal/mol\n", block->kc);
  printf("Coverage:      %f\n", block->chi);
  printf("Spontaneous Curvature:      %f microns^-1\n", block->c0);
  printf("Area per lipid:      %e microns^2\n", block->Ap);
  printf("Diffusion Costant:      %f microns^2/s\n", block->D0);
  printf("Water Viscosity at 25C:      %f kcal s/mol microns^3\n", block->eta);
  printf("Mode Max:      %f\n", block->mode_max);
  printf("Nstesp:      %d\n", block->Nsteps);
  printf("rate:      %f\n", block->rate);
  printf("/****************************/\n");
}
