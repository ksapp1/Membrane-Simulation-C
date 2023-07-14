#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "util.h"
#include <stdarg.h>

void getLine( FILE *theFile, char *theBuffer )
{
        int i = 0;

        while( !feof(theFile) )
        {
                char tc = fgetc(theFile);

                if( tc != '\n' && i < 50000 )
                {
                        theBuffer[i++] = tc;
                }
                else if( tc != '\n' && i >= 50000 )
                {
                }
                else
                        break;
        }

        theBuffer[i] = '\0';
}
