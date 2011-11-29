/* Copyright 2010 (c) Jostein R. Natvig <Wbfgrva.angivt@tznvy.pbz> */

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <cassert>

#include <vector>
#include <cmath>
/* ------------------------------------------------------------------ */
/*      Interface for parsing records of '/'-separated strings        */
/* ------------------------------------------------------------------ */
#define BUFSIZE 1024

struct State {
    int error;
    int lineno;
};

static int   parse_next_token   (FILE *fp, struct State *s, char **t, int *c);
static char *get_token_string   (FILE *fp, struct State *s, char *buf);
static void  split_token_string (char *str, char *rstr, char *vstr);
static int   get_default_count  (char *rstr, struct State *s);
static void  skip_line          (FILE *fp, struct State *s);
static char *string_to_int      (int *value, char *str);


static int 
readtokens(FILE *fp, char **tokens, int numtokens, struct State *s)
{
    char  *token;
    int    defaultcount;
    int    i;
    int    pos;
    int finish;

    defaultcount   = 0;
    s->error = 0;
    pos     = 0;
    
    while ((finish=parse_next_token(fp, s, &token, &defaultcount)) && pos < numtokens)
    {
        if (defaultcount == 0)
        {
            tokens[pos++] = token;
        }
        else
        {
            if (defaultcount > numtokens - pos)
            {
                fprintf(stderr, "Too many default values '%s' on line %d.  ", token, s->lineno);
                break;
            }

            for (i=0; i<defaultcount; ++i)
            {
                tokens[pos++] = NULL; /* strdup("DEFAULT"); */
            }

            free(token);
        }
    }

    for (; pos < numtokens; ++pos)
    {
        tokens[pos] = NULL; /* strdup("DEFAULT"); */
    }

    if (finish != 0)
    {
        fprintf(stderr, "Too many tokens?\n");
    }
    return finish == 0;
}

/* ------------------------------------------------------------------ */

static int 
parse_next_token(FILE *fp,  struct State *s, char **tok, 
                            int *count)
{
    char str[BUFSIZE];
    char rstr[BUFSIZE];
    char vstr[BUFSIZE];
    int  ok = 1;
    
    str[0]='\0';
    while (strlen(str)==0)
    {
        get_token_string(fp, s, str);       
    }

    if (str[0] == '/')
    {
        ok = 0;
    }
    else
    {
        *tok = strdup(str);
        *count = 0;

        split_token_string(str, rstr, vstr);

        if ((strlen(rstr)>0) && (strlen(vstr)==0))
        {
            *count = get_default_count(rstr, s);
            
            if (s->error)
            {
                fprintf(stderr, "Warning: string token ends with a '*'\n");

                /* attempt to proceed */
                s->error = 0;
                *count = 0;
            }
        }
    }
    return ok;
}

/* ------------------------------------------------------------------ */

static void 
split_token_string(char *str, char *rstr, char *vstr)
{
    /* str:'rrrr*vvvv' -> rstr:'rrrr' vstr:'vvvv' */
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
    
    while((*vstr++=*str++))
    { 
        ; 
    }
    *vstr='\0';
}

/* ------------------------------------------------------------------ */

static int 
get_default_count(char *rstr, struct State *s)
{
    int   r;
    char *ptr;

    if (strlen(rstr)>0)
    {
        ptr=string_to_int(&r, rstr);
        
        /* no conversion, extra characters, ... */
        if ((ptr==rstr) || (strlen(ptr)>0) || (r<1))
        {
            s->error = 1;
        }
    }
   else
   {
       r = 0;
   }

   return r;
}

/* ------------------------------------------------------------------ */
/* hint: always read '/' as a separate token */

static char *
get_token_string(FILE *fp, struct State *s, char *str)
{
    char *ptr = str;
    int   c;
    int   quote=0;
    
    /* skip leading blanks */
    while((c = fgetc(fp)) != EOF && isspace(c)) {if (c=='\n') s->lineno++;}

    
    /* catch end marker */
    if (c == '/') {
        skip_line(fp, s);
        *ptr++ = c;
        *ptr++ = '\0';
        return str;
    }
    

    /* otherwise, read the wole token in the next block*/
    ungetc(c, fp);

    while((c = fgetc(fp)) != EOF && (quote || !isspace(c)))
    {
        *ptr++ = c;
        if (c=='\''){  
            quote = quote ? 0 : 1;
        }

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
        
        assert(ptr-str < BUFSIZE);
    }
    
    if (c=='\n') 
    {
        s->lineno++;
    }
    
    *ptr='\0';
    return str;
}

/* ------------------------------------------------------------------ */

static void 
skip_line(FILE *fp, struct State *s)
{
    int c;
    while((c = fgetc(fp))!=EOF && c != '\n') {
        ;
    }
    ++s->lineno;
 
}

/* ------------------------------------------------------------------ */

static char *
string_to_int(int *value, char *str)
{
    char *q;
    *value = strtol(str, &q, 10);
    return q;
}

/* ------------------------------------------------------------------ */

static int 
all_defaults(char **t, int n)
{
    int i;
    for (i=0; i<n; ++i)
    {
        if (t[i] != NULL)
        {
            return 0;
        }
    }
    return 1;
}

/* ------------------------------------------------------------------ */
// Convert this to std::vector<std::string>
struct Record {
    char **fields;
    size_t numfields;

    Record(size_t n)
    {
	numfields = n;
	fields = new char*[numfields];
	for(size_t i=0; i<numfields; ++i)
	{
	    fields[i] = NULL;
	}
    }
    ~Record() {
	for (size_t i=0; i<numfields; ++i)
        {
            free(fields[i]);
            fields[i] = NULL;
        }
	delete fields;
    }
    Record(const Record& R)
    {
	numfields = R.numfields;
	fields = new char*[numfields];
	for(size_t i=0; i<numfields; ++i)
	{
	    if (R.fields[i] != NULL)
	    {
		fields[i] = strdup(R.fields[i]);
	    }
	    else
	    {
		fields[i] = NULL;
	    }
	}
    }

    void display(FILE *stream)
    {
	fprintf(stream, "[");
	for(size_t i=0; i<numfields; ++i)
	{
	    fprintf(stream, "%s, ", fields[i]);
	}
        fprintf(stream, "\b\b]\n");
    }
};

/* ------------------------------------------------------------------ */

void read_records(const char *fn, size_t record_width, std::vector<Record> &v)
{
    Record  r(record_width);
    struct State   s = {0, 1};
    
    FILE *fp = fopen(fn, "r");
    int   ok;
    
    /* read until am empty record is found (which is read as 'only
     * default values') */
    while ( (1 == (ok=readtokens (fp, r.fields, r.numfields, &s))) && 
            (0 == all_defaults (r.fields, r.numfields))   )
    {
	v.push_back(r);
    }

    fclose(fp);
}

void 
display_records(std::vector<Record> v, FILE *stream)
{
    for (size_t i=0; i<v.size(); ++i)
    {
	v[i].display(stream);
    }
}


/* ------------------------------------------------------------------ */

void
read_compdat(const char *fn, std::vector<Record> &v)
{
    size_t sz = 14;
    read_records(fn, sz, v);
}
void
read_welspecs(const char *fn, std::vector<Record> &v)
{
    size_t sz = 16;
    read_records(fn, sz, v);
}
void
read_wconinje(const char *fn, std::vector<Record> &v)
{
    size_t sz = 14;
    read_records(fn, sz, v);
}

static char *trimspace (char *str)
{
    char *ibuf, *obuf;
   
    if (str) {
        for (ibuf = obuf = str; *ibuf; ) {
            while (*ibuf && (isspace (*ibuf))) {
                ibuf++;
            }
            if (*ibuf && (obuf != str)) {
                *(obuf++) = ' ';
            }
            while (*ibuf && (!isspace (*ibuf))) {
                *(obuf++) = *(ibuf++);
            }
        }
        *obuf = '\0';
    }
    return (str);
}
static char* remove_quotes(char *str)
{
    int len = strlen(str);
    if ((str[0]=='\'') && (str[strlen(str)-1]=='\''))
    {
        str[0] = ' ';
        str[len-1] = '\0';
        trimspace(str);
    }
    return str;
}
static void get_wellnames(std::vector<Record> welspecs, Record &wnames)
{
    int n = welspecs.size();

    for (int i=0; i<n; ++i)
    {
	assert(welspecs[i].fields[0] != NULL);
	
	wnames.fields[i] = strdup(remove_quotes(welspecs[i].fields[0]));
    }
}
static int get_well_number(Record wnames, char *name)
{
    size_t i;
    for(i=0; i<wnames.numfields; ++i)
    {
	if (strcmp(wnames.fields[i], remove_quotes(name))!=0)
	{
	    break;
	}
    }
    return i;
}
static double effective_pressure_radius(const Record &completion, 
                                        double dx, double dy, double dz,
                                        double kx, double ky, double kz)
{
    double k1,k2,dx1,dx2;
    double r0;
    
//    assert( completion != NULL);

    if (completion.fields[13] != NULL)
    {
        sscanf(completion.fields[13], "%lf", &r0);
    }
    else
    {
        /* if a field is NULL, use the default value */
        if ((completion.fields[12] == NULL) || 
            (strcmp(completion.fields[12], "Z")==0))
        {
            dx1=dx;dx2=dy;k1=kx;k2=ky;                
        }
        else if (strcmp(completion.fields[12], "X")==0)
        {
            dx1=dy;dx2=dz;k1=ky;k2=kz; 
        }
        else if(strcmp(completion.fields[12], "Y")==0)
        {
            dx1=dx;dx2=dz;k1=kx;k2=kz;
        }
        else
        {
            fprintf(stderr, "ERROR: unsupported completion direction (item 13)\n");
            exit(EXIT_FAILURE);
        }
        
        /* Effective pressure radius used in Peaceman formula for
         * rectangular grid block and permeability aligned with axes
         * of grid cell */
        r0  = 0.28;
        r0 *= sqrt(sqrt(k2/k1)*dx1*dx1 + sqrt(k1/k2)*dx2*dx2);
        r0 /= sqrt(sqrt(k1/k2)) + sqrt(sqrt(k2/k1));
    }

    return r0;
}                                        
static double compute_transmissibility(const Record &completion)
{
    double r; 
    double r0;
    double kh;
    double skin = 0;
    double twopi = 6.283185307179586;

    double trans = -1.0;

    if (completion.fields[7] != NULL )
    {
                sscanf(completion.fields[7], "%lf", &trans);
    } 

    
    if (completion.fields[7] == NULL || trans == 0)
    {

        if(completion.fields[8] !=NULL)
        {
            sscanf(completion.fields[8], "%lf", &r);
        }
        else
        {
            r = 0.3048;
        }

        if(completion.fields[9] !=NULL)
        {
            sscanf(completion.fields[9], "%lf", &kh);

        }
        else
        {
            kh = 1;
        }

        if(completion.fields[10]!=NULL)
        {
            sscanf(completion.fields[10], "%lf", &skin);
        }
        else
        {
            skin = 0;
        }

        r0 = effective_pressure_radius(completion, 100,100,1,10,10,1);     

        trans = twopi*kh/(logf(r0/r)+skin);

    }
    return trans;
}
static void get_well_completions(const Record &completion, 
				 std::vector<Record> welspecs, 
				 int dims[3],
				 std::vector<std::vector<int> > &completions)
{
    int k;
    int i,j,k1,k2;
    int pcell;
    Record wnames(welspecs.size());
    
    get_wellnames(welspecs, wnames);
    int wnum = get_well_number(wnames, completion.fields[0]);

    assert(completion.fields[1]!=NULL);
    assert(completion.fields[2]!=NULL);
    assert(completion.fields[3]!=NULL);
    assert(completion.fields[4]!=NULL);

    sscanf(completion.fields[1], "%d", &i);
    sscanf(completion.fields[2], "%d", &j);
    sscanf(completion.fields[3], "%d", &k1);
    sscanf(completion.fields[4], "%d", &k2);
    
    

    compute_transmissibility(completion);
    
    fprintf(stderr, "%d %d %d %d: ", i,j,k1,k2);
    
    i = i==0 ? 0 : i-1;
    j = j==0 ? 0 : j-1;
    k1 = k1==0 ? 0 : k1-1;
    k2 = k2==0 ? 0 : k2-1;

    for (k = k1; k<=k2; ++k)
    {
        pcell = i + dims[0]*(j + dims[1]*k);
	completions[wnum].push_back(pcell);
    }
    fprintf(stderr, "\n");

}


int main()
{
    std::vector<Record> compdat;
    read_compdat("compdat.txt", compdat);
    display_records(compdat, stderr);
    fprintf(stderr, "\n\n");
    
    std::vector<Record> welspecs;
    read_welspecs("welspecs.txt", welspecs);
    display_records(welspecs, stderr);
    fprintf(stderr, "\n\n");
   
    std::vector<Record> wconinje;
    read_wconinje("wconinje.txt", wconinje);
    display_records(wconinje, stderr);
 
    Record wnames(welspecs.size());
    get_wellnames(welspecs, wnames);
    wnames.display(stderr);
    
    int dims[3] = {10,10,5};
    std::vector<std::vector<int> > completions(2);
    for (size_t i=0; i<compdat.size(); ++i)
    {
	get_well_completions(compdat[i], welspecs, dims, completions);
    }

    fprintf(stderr, "here:\n");
    for (size_t i=0; i<completions.size(); ++i)
    {
	fprintf(stderr, "%d: ", (int)i);
	for (size_t j=0; j<completions[i].size(); ++j)
	{
	    fprintf(stderr, "%d ", completions[i][j]);
	}
	fprintf(stderr, "\n");
    }
    return 0;
}
