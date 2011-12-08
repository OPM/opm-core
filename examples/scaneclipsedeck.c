/* scaneclipsedec finds technically valid Eclipse keywords in an ascii file.
 * Copyright (c) 2010 Jostein R. Natvig <jostein.natvig@gmail.com>
 *
 * The MIT License
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>



static char* 
read_keyword(FILE *fp, char *buf)
{
    int i, j, c;
    
    /* Clear buf */
    for (i=0; i<9; ++i) {
        buf[i] = '\0';
    }

    /* Read first character and check if it is uppercase*/
    buf[0] = fgetc(fp);
    if ( !isupper( buf[0] ) ) {        
        return NULL;          /* NOT VALID CHARACTER */
        ungetc(buf[0], fp);
    }
    
    
    /* Scan as much as possible possible keyword, 8 characters long */
    i = 1;
    while ( (c = fgetc(fp)) && 
            (c != EOF     ) &&  
            (c != '\n'    ) && 
            (c != '/'     ) && 
            (i < 8        )) {
        buf[i++] = c;
    }

    /* Skip rest of line */
    if (c != '\n'){
        while ( (c = fgetc(fp)) && 
                (c != EOF     ) &&  
                (c != '\n'    )) {
            ;
        }
    }

    if(c == '\n') {
        ungetc(c, fp);
    }
    
    /* Find first non-uppercase or non-digit character */ 
    for (i=0; i<8; ++i) {
        if ( !(isupper(buf[i]) || isdigit(buf[i])) ) {                    
            break;
        }
    }
    
    /* Check if remaining characters are blank */
    for (j = i; j<8; ++j) {
        if(!isspace(buf[j]) && buf[j] != '\0') {
            return NULL; /* CHARACTER AFTER SPACE OR INVALID CHARACTER */
        }
        buf[j] = '\0';
    }
    return buf;
}


int main(int argc, char *argv[])
{
    FILE         *fp;
    char buf[10];
    
    if (argc != 2)
    {
        fprintf(stderr, "Usage: <app> filename.grdecl\n");
        exit(EXIT_FAILURE);
    }
    fp = fopen(argv[1], "ra");
    if (fp == NULL)
    {
        fprintf(stderr, "No such file...\n");
        exit(EXIT_FAILURE);
    }

    int c, lineno = 0, nkw = 0;
    
    if (read_keyword(fp, buf) != NULL) {
        ++nkw;
        fprintf(stderr, "%s\n", buf);
    }
                


    while ((c = getc(fp)) != EOF) {    /* Eat large chunks */
        if ( c == '\n') {
            ++lineno;
        
            if (read_keyword(fp, buf) != NULL) {
                fprintf(stderr, "%s\n", buf);
                ++nkw;
            }
        }
    }
    
    fprintf(stderr, "Scanned %d lines, found %d keywords.\n", lineno, nkw);

    return 0;
}
/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
