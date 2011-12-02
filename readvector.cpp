/*===========================================================================
//
// Author: Jostein R. Natvig     <Jostein.R.Natvig@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.
*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include <vector>

#include "readvector.hpp"

struct ParserState {
    int             error;
};

template <typename T>
static void  parse_next_token  (FILE *fp,  struct ParserState *s, 
				T *v, int *count);
static char *get_token_string  (FILE *fp, struct ParserState *s, char *buf);
static void  skip_line         (FILE *fp, struct ParserState *s);
static char *convert_string    (int *value, char *str);
static char *convert_string    (double *value, char *str);
static void  split_token_string(char *str, char *rstr, char *vstr);
static int   get_repeat_count  (char *rstr, struct ParserState *s);
template <typename T>
static void get_value          (char *vstr, struct ParserState *s, T *v);

/* ------------------------------------------------------------------ */

static const int bufsz = 1024;

/* ------------------------------------------------------------------ */

template <typename T>
static void readvector(FILE *fp, struct ParserState *s, std::vector<T>& v)
{
    T   value;
    int count, i;

    count = 1;
    s->error = 0;
    while (count > 0){

        parse_next_token(fp, s, &value, &count);
        for (i=0; i<count; ++i)
        {           
	    v.push_back(value);
        }
    }

    if (s->error) {
        assert(0);
    }
}


/* ------------------------------------------------------------------ */
template <typename T>
static void parse_next_token(FILE *fp,  struct ParserState *s, T *v, int *r)
{
    char str[bufsz];
    char rstr[bufsz];
    char vstr[bufsz];

    str[0]='\0';
    while (strlen(str)==0)
    {
        get_token_string(fp, s, str);       
    }

    if (str[0] == '/')
    {
        /* end-of-vector marker */
        *r = 0;
    }
    else
    {
        split_token_string(str, rstr, vstr);
        *r = get_repeat_count(rstr, s);
	get_value(vstr, s, v);
        
        if (s->error) 
        {
            /* signal abort to caller */
            *r = 0;
        }
    }
}

/* ------------------------------------------------------------------ */

/* split string 'rrrr*vvvv' in strings 'rrrr' and 'vvvv' */
static void split_token_string(char *str, char *rstr, char *vstr)
{
    char *ptr;

    ptr=strchr(str, '*');
    if ((ptr != NULL) && (ptr!=str)) {
        while(str!=ptr) 
        { 
            *rstr++=*str++; 
        }
        *rstr='\0';
        
        str++;
    }
    else
    {
        rstr[0]='\0';
    }
    
    while((*vstr++=*str++) )
    { 
        continue ; 
    }
    *vstr='\0';
}

/* ------------------------------------------------------------------ */

static int get_repeat_count(char *rstr, struct ParserState *s)
{
    int   r;
    char *ptr;

    if (strlen(rstr)>0)
    {
        ptr=convert_string(&r, rstr);
        
        if ((ptr==rstr) || (strlen(ptr)>0) || (r<1))
        {
            s->error = 1;
        }
    }
   else
   {
       r = 1;
   }

   return r;
}

/* ------------------------------------------------------------------ */
template <typename T>
static void get_value(char *vstr, struct ParserState *s, T *v)
{
    char *ptr;

    if (strlen(vstr)>0)
    {
        /* Convert FORTRAN style floats 3.13D-6 to IEEE */
        for(ptr=vstr; *ptr; ++ptr)
        {
            if (*ptr=='D')
            {
                *ptr='E';
            }
        }
        
        ptr = convert_string(v, vstr); 

        if ((ptr==vstr) || (strlen(ptr)>0))
        {
            s->error = 1;
        }
    }
    else
    {
        s->error = 1;
    }
}


/* ------------------------------------------------------------------ */

static char *get_token_string(FILE *fp, struct ParserState *s, char *str)
{
    char *ptr = str;
    int   c;
    
    /* Skip leading blanks */
    while((c = fgetc(fp)) != EOF && isspace(c)) 
    {
	;
    }

    *ptr++ = c;

    /* Catch end marker */
    if (c == '/') {
        skip_line(fp, s);
        *ptr++ = '\0';
        return str;
    }

    while((c = fgetc(fp)) != EOF && !isspace(c)) 
    {
        *ptr++ = c;
        
        /* Break and skip rest of line if '--' if encountered */
        if (*(ptr-1) == '-' && ptr-str>1 && *(ptr-2) == '-'){
            ptr = ptr - 2;
            skip_line(fp, s);
            break;
        }
        
        /* If end marker is encontered, push character back onto stream. */
        if (c=='/') {
            ungetc(c, fp);
            ptr--;
            break;
        }
        
        assert(ptr-str < bufsz);
    }
    
    *ptr='\0';
    return str;
}

/* ------------------------------------------------------------------ */

static void skip_line(FILE *fp, struct ParserState *s)
{
    static_cast<void>(s);
    int c;
    while((c = fgetc(fp))!=EOF && c != '\n') {
        ;
    }
}

/* ------------------------------------------------------------------ */

// int version
static char *convert_string(int *value, char *str)
{
    char *q;
    *value = strtol(str, &q, 10);
    return q;
}

/* ------------------------------------------------------------------ */

// double version
static char *convert_string(double *value, char *str)
{
    char *q;
    *value = strtod(str, &q);
    return q;
}

/* ------------------------------------------------------------------ */
template <typename T>
void read_vector_from_file(const std::string& fn, std::vector<T>& v)
{
    FILE *fp = fopen(fn.c_str(), "r");
    struct ParserState s = { 0 };
    readvector(fp, &s, v);
    fclose(fp);
}


void read_vector_from_file(const std::string& fn, std::vector<int>& v)
{
    read_vector_from_file<int>(fn, v);
}

void read_vector_from_file(const std::string& fn, std::vector<double>& v)
{
    read_vector_from_file<double>(fn, v);
}
