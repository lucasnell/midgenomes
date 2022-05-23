/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------- NPStat ---------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/* Code to extract pool statistics (theta, neutrality tests)
from pileup data and a fasta file representing the outgroup sequence */

/* Compile with
gcc -O3 -o npstat NPStat-vXX.c -lgsl -lgslcblas -lm
substituting XX with version number
*/

/* Arguments:
run "npstat" to see the arguments
*/

/* Output stats:
window number, bases (after filtering),
 bases with known outgroup allele (after filtering),
 average read depth, number of segregating sites, Watterson theta, Tajima's Pi,
 Tajima's D, Fay and Wu's H, divergence per base (from outgroup), HKA etc...
*/



/* Include libraries */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>


/* Define substitutions */

#define DEB(x) //x

// this fixes the minimum acceptable coverage in terms of reads aligned:
#define minimum_coverage 2
// this fixes the maximum acceptable coverage in terms of reads aligned:
#define maximum_coverage 1000
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

/* Declare structures */

struct tests
{
    unsigned long cov;
    unsigned long l;
    unsigned long l_out;
    unsigned long s;
    double num_t;
    double num_p;
    double num_hl;
    double num_hq;
    double den_t;
    double den_p;
    double den_hl;
    double den_hq;
};

struct fst_calc
{
    unsigned long l;
    double gen_diff;
    double c_s;
};

struct combinatorial
{
    double *d_t;
    double *d_p;
    double *d_hl;
    double *d_hq;
};

struct combinatorial_fst
{
    double **c_s;
};

/* Declare functions */

int generate_pool_covariance_matrix(double *covmat, double *covmatpool,
                                    int na, int nb, int n);
int generate_covariance_matrix(double *covmat, int n);




int base_to_num(char base);

int read_line_pileup(FILE * bam_file, unsigned long min_qual,
                     unsigned long min_mqual, unsigned long * pos_base,
                     unsigned long * n_ref, unsigned long * n_alt_allele,
                     unsigned long * n_tot_allele, unsigned long * n_alt,
                     int * ref_base, int *  alt_base) ;

int extract_outgroup_base(FILE * fasta_out, unsigned long pos,
                          unsigned long oldpos, int fasta_length);

void extract_stats(struct tests * test, struct combinatorial * comb, int n0,
                   unsigned long n_ref, unsigned long n_alt_allele,
                   unsigned long rd, unsigned long * n_alt, int ref_base,
                   int alt_base, int out_base, int mb);

void extract_fst(struct fst_calc * fst, struct combinatorial_fst * combfst,
                 int n01, unsigned long n_ref1, unsigned long n_alt_allele1,
                 unsigned long rd1, unsigned long * n_alt_1, int ref_base1,
                 int alt_base1, int n02, unsigned long n_ref2,
                 unsigned long n_alt_allele2, unsigned long rd2,
                 unsigned long * n_alt_2, int ref_base2, int alt_base2,
                 int out_base, int mb);

//SNPS
void extract_snps(unsigned long pos, FILE * output_snps,
                  unsigned long * n_alt_1, unsigned long * n_alt_2, int mb);

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*-------------------------- Functions -------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/


int generate_covariance_matrix(double *covmat, int n)
{
    double s,*a,*b,bb,bz,zz,a2,*xi_null2,*xi_b2;
    int i,j,n2,ca,cb;
    double theta;
    gsl_matrix *c,*invc,*sigma,*ct,*invct;

    if(n<2) return -10000;

    n2=ceil((double)n/2);
    n=(int)n2*2;

    a = (double *)calloc((unsigned long int)n+1,sizeof(double));
    b = (double *)calloc((unsigned long int)n,sizeof(double));

    a[0]=0;
    a2=0;
    for(i=2;i<=n+1;i++){
        a[i-1]=a[i-2]+1/(double)(i-1);
    }
    for(i=1;i<n;i++){
        a2+=1/(double)(i*i);
        b[i-1]=2*(double)n*(a[n]-a[i-1])/(double)((n-i+1)*(n-i))-2/(double)(n-i);
    }
    b[n-1]=0;
    sigma=gsl_matrix_alloc(n-1,n-1);
    gsl_matrix_set_zero(sigma);
    for (i=1;i<n;i++){
        if (2*i<n) {
            gsl_matrix_set(sigma,i-1,i-1,b[i]);
        }  else {
            if (2*i>n) {
                gsl_matrix_set(sigma,i-1,i-1,b[i-1]-1/gsl_pow_2((double)i));
            }  else {
                gsl_matrix_set(sigma,i-1,i-1,(a[n-1]-a[i-1])*2/(double)(n-i)-1/
                    gsl_pow_2((double)i));
            }
        }
    }
    for (i=1;i<n;i++){
        for (j=1;j<i;j++){
            if (i+j<n) {
                gsl_matrix_set(sigma,i-1,j-1,(b[i]-b[i-1])/2);
                gsl_matrix_set(sigma,j-1,i-1,(b[i]-b[i-1])/2);
            }  else {
                if (i+j>n) {
                    gsl_matrix_set(sigma,i-1,j-1,(b[j-1]-b[j])/2-1/(double)(i*j));
                    gsl_matrix_set(sigma,j-1,i-1,(b[j-1]-b[j])/2-1/(double)(i*j));
                }  else {
                    gsl_matrix_set(sigma,i-1,j-1,(a[n-1]-a[i-1])/(double)(n-i)+
                        (a[n-1]-a[j-1])/(double)(n-j)-(b[i-1]+b[j])/2-1/(double)(i*j));
                    gsl_matrix_set(sigma,j-1,i-1,(a[n-1]-a[i-1])/(double)(n-i)+
                        (a[n-1]-a[j-1])/(double)(n-j)-(b[i-1]+b[j])/2-1/(double)(i*j));
                }
            }
        }
    }

    for(ca=1;ca<n;ca++){
        for(cb=1;cb<n;cb++){
            covmat[(cb-1)*(n-1)+ca-1]=gsl_matrix_get(sigma,ca-1,cb-1);
        }
    }


    free(a);
    free(b);

    gsl_matrix_free(sigma);
    return 1;
}



int base_to_num(char base)
{
    switch(base)
    {
    case 'A': return 1; break;
    case 'a': return 1; break;
    case 'C': return 2; break;
    case 'c': return 2; break;
    case 'G': return 3; break;
    case 'g': return 3; break;
    case 'T': return 4; break;
    case 't': return 4; break;
    default: return 0;
    };
};


/*--------------------------------------------------------------*/

int read_line_pileup(FILE * bam_file, unsigned long min_qual,
                     unsigned long min_mqual, unsigned long * pos_base,
                     unsigned long * n_ref, unsigned long * n_alt_allele,
                     unsigned long * n_tot_allele, unsigned long * n_alt,
                     int * ref_base, int * alt_base) {

    int count_i;
    char *cline, *cchrom, *cpileup, *cqual, *cmqual, crefbase;
    unsigned long count_j, n_ins;
    unsigned long nseq;
    size_t nline;

    DEB(printf("entering read routine\n")); //debug

    cchrom=(char *) malloc(100);
    cpileup=(char *) malloc(40);
    cqual=(char *) malloc(20);
    cmqual=(char *) malloc(20);
    cline=(char *) malloc(1);

    nline=1;

    getdelim(&cline,&nline,9,bam_file);
    sscanf(cline,"%s\t",cchrom);
    getdelim(&cline,&nline,9,bam_file);
    sscanf(cline,"%lu\t",pos_base);
    getdelim(&cline,&nline,9,bam_file);
    sscanf(cline,"%c\t",&crefbase);
    getdelim(&cline,&nline,9,bam_file);
    sscanf(cline,"%lu\t",&nseq);
    getdelim(&cline,&nline,9,bam_file);
    if (strlen(cline)>=40) cpileup=(char *)realloc(cpileup,strlen(cline)+1);
    sscanf(cline,"%s\t",cpileup);
    getdelim(&cline,&nline,10,bam_file);
    if (strlen(cline)>=20) cqual=(char *)realloc(cqual,strlen(cline)+1);
    sscanf(cline,"%s\t",cqual);

    DEB(printf("read data %s\t%lu\t%c\t%lu\t%s\t%s\t%s\n",cchrom,*pos_base,
               crefbase,nseq,cpileup,cqual,cmqual)); //debug

    for(count_i=0;count_i<5;count_i++){
        n_alt[count_i]=0;
    };
    count_j=0;
    *n_ref=0;
    *n_alt_allele=0;
    *n_tot_allele=0;

    if((nseq>=minimum_coverage)&&(nseq<=maximum_coverage)) {

        for(count_i=0;count_i<strlen(cpileup);count_i++){

            switch (cpileup[count_i]) {
            case '^': count_i++; break;
            case '$': break;
            case '*': count_j++; break;
            case '+': {
                count_i++;
                for(n_ins=0;isdigit(cpileup[count_i])!=0;count_i++){
                    n_ins=n_ins*10+(cpileup[count_i]-48);
                };
                for(n_ins=n_ins-1;n_ins>0;n_ins--) {
                    count_i++;
                };
            }; break;
            case '-': {
                count_i++;
                for(n_ins=0;isdigit(cpileup[count_i])!=0;count_i++){
                    n_ins=n_ins*10+(cpileup[count_i]-48);
                };
                for(n_ins=n_ins-1;n_ins>0;n_ins=n_ins-1) {
                    count_i++;
                };
            }; break;
            case 'N': count_j++; break;
            case 'n': count_j++; break;
            case '.': {if((cqual[count_j]>=min_qual+33)) (n_alt[0])++; count_j++;}; break;
            case ',': {if((cqual[count_j]>=min_qual+33)) (n_alt[0])++; count_j++;}; break;
            case 'A': {if((cqual[count_j]>=min_qual+33)) (n_alt[1])++; count_j++;}; break;
            case 'C': {if((cqual[count_j]>=min_qual+33)) (n_alt[2])++; count_j++;}; break;
            case 'G': {if((cqual[count_j]>=min_qual+33)) (n_alt[3])++; count_j++;}; break;
            case 'T': {if((cqual[count_j]>=min_qual+33)) (n_alt[4])++; count_j++;}; break;
            case 'a': {if((cqual[count_j]>=min_qual+33)) (n_alt[1])++; count_j++;}; break;
            case 'c': {if((cqual[count_j]>=min_qual+33)) (n_alt[2])++; count_j++;}; break;
            case 'g': {if((cqual[count_j]>=min_qual+33)) (n_alt[3])++; count_j++;}; break;
            case 't': {if((cqual[count_j]>=min_qual+33)) (n_alt[4])++; count_j++;}; break;
            };


        };

        *n_alt_allele = 0;
        *n_ref=0; //n_alt[0];
        for(count_i=0;count_i<5;count_i++){
            if (*n_ref < n_alt[count_i]) {
                *n_ref = n_alt[count_i]; *ref_base=count_i;
            };
        };
        for(count_i=0;count_i<5;count_i++){
            if ((*n_alt_allele < n_alt[count_i])&&(count_i!=*ref_base)) {
                *n_alt_allele = n_alt[count_i]; *alt_base=count_i;
            };
        };
        if (*ref_base==0) { *ref_base=base_to_num(crefbase); };
        if (*alt_base==0) { *alt_base=base_to_num(crefbase); };
        *n_tot_allele = *n_ref + *n_alt_allele;


    };

    DEB(printf("using data %lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%u\t%u\n",
               *n_ref, *n_alt_allele, *n_tot_allele, n_alt[0], n_alt[1],
               n_alt[2], n_alt[3], n_alt[4], *ref_base, *alt_base)); //debug

    DEB(printf("exit read routine\n")); //debug
    free(cchrom);
    free(cpileup);
    free(cqual);
    free(cmqual);
    free(cline);

};








/*--------------------------------------------------------------*/


int extract_outgroup_base(FILE * fasta_out, unsigned long pos,
                          unsigned long oldpos, int fasta_length)
{
    unsigned long diff;
    ldiv_t pos_newlines, oldpos_newlines;
    pos_newlines=ldiv(pos-1,(unsigned long)fasta_length);
    oldpos_newlines=ldiv(oldpos-1,(unsigned long)fasta_length);
    diff=pos-oldpos+pos_newlines.quot-oldpos_newlines.quot;
    fseek(fasta_out,diff-1,SEEK_CUR);
    return base_to_num(fgetc(fasta_out));
};





/*--------------------------------------------------------------*/


void extract_stats(struct tests * test, struct combinatorial * comb, int n0,
                   unsigned long n_ref, unsigned long n_alt_allele,
                   unsigned long rd, unsigned long * n_alt, int ref_base,
                   int alt_base, int out_base, int mb)
{
    char is_out=0;

    if (out_base>0)
    {
        if (ref_base==out_base) { is_out=1; }
        else if (n_alt[out_base]==n_alt_allele) { is_out=2; };
    };

    test->cov+=rd;
    test->l+=1;
    test->den_t+=comb->d_t[rd-1];
    test->den_p+=comb->d_p[rd-1];
    if (is_out>0)
    {
        test->l_out+=1;
        test->den_hl+=comb->d_hl[rd-1];
        test->den_hq+=comb->d_hq[rd-1];
    };

    if ((n_ref>mb)&&(n_ref<rd-mb))
    {
        test->s+=1;
        test->num_t+=1;
        test->num_p+=(double)(2*n_alt_allele*n_ref)/(double)(rd*(rd-1));
        if (is_out==1){
            test->num_hl+=(double)(n_alt_allele)/(double)(rd-1);
            test->num_hq+=(double)(n_alt_allele*n_alt_allele)/(double)(rd*(rd-1));
        };
        if (is_out==2)
        {
            test->num_hl+=(double)(n_ref)/(double)(rd-1);
            test->num_hq+=(double)(n_ref*n_ref)/(double)(rd*(rd-1));
        };
    };


    DEB(printf("using data %u\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%u\t%u\t%u\n",
               n0, n_ref, n_alt_allele, rd, n_alt[0], n_alt[1],n_alt[2],
               n_alt[3],n_alt[4], ref_base, alt_base, out_base));

};



/*--------------------------------------------------------------*/

unsigned long extract_pos_snpinput( FILE * snpinput){
    char *line;
    int c,i,j;
    size_t l_line;
    unsigned long pos;
    line = (char *) malloc (100*sizeof(char));
    l_line=100;
    char test; test=fgetc(snpinput); if (test!=EOF) {
        ungetc(test, snpinput);
        getline(&line,&l_line,snpinput);
        sscanf(line,"%lu",&pos);
    } else { pos=2000000000; };
    free(line);
    return pos;
};




//SNPS
void extract_snps(unsigned long pos, FILE * output_snps,
                  unsigned long * n_alt_1, unsigned long * n_alt_2,
                  int mb)
{
    int i, c1, c2;

    c2=0;
    c1=0;
    for(i=1;i<5;i++)
    {
        if ((n_alt_1[i]+n_alt_2[i])>(n_alt_1[c1]+n_alt_2[c1]))
        {
            c1=i;
        };
    };
    if (c1==0) c2=1;
    for(i=c2+1;i<5;i++)
    {
        if (((n_alt_1[i]+n_alt_2[i])>(n_alt_1[c2]+n_alt_2[c2]))&&(i!=c1))
        {
            c2=i;
        };
    };
    if((n_alt_1[c2]+n_alt_2[c2])>1)
    {
        fprintf(output_snps, "%lu\t%lu\t%lu\t%lu\t%lu\n", pos, n_alt_1[c2],
                n_alt_2[c2], n_alt_1[c1], n_alt_2[c1]);
    };
};

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*------------------------ Main program ------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/



int main(int argc, char *argv[])
{
    /*
     arguments: bam file 1, bam file 2, fasta outgroup file, window size,
     haploid sample size, minimum coverage, maximum coverage,
     minimum base quality, minimum mapping quality,
     low frequency alleles removed
     */

    /* Input flags */
    double input_coverage1;
    /* Variables */

    int count_i,count_j;
    FILE *bam_file, *bam_file1, *bam_file2, *fasta_ref, *fasta_out;
    char ct, ct1, ct2;
    int lt;
    unsigned long line, line1, line2, pos, oldpos, pos_base;
    unsigned long pos_base1, pos_base2, nseq;
    unsigned long n_ref1, n_alt_allele1, n_tot_allele1;
    unsigned long n_alt_1[5];
    unsigned long n_ref2, n_alt_allele2, n_tot_allele2;
    unsigned long n_alt_2[5];
    unsigned long max_cov, min_cov, min_qual, min_mqual;
    int m_bar, m1t, m2t;
    unsigned long window_size;
    int n0max, n01, n02;
    unsigned long rd, rd1, rd2;
    FILE *output_stat1, *output_stat2, *output_fst;
    char *temp_file_name1, *temp_file_name2;
    //unsigned long
    int fasta_length;
    int ref_base, ref_base1, ref_base2, alt_base1, alt_base2, out_base;
    /* Gipo's variables */
    long int n_pos1;
    int array_counts1[2000] = {0};
    unsigned long int int_positions1[1000] = {0};
    int n_variation1;
    n_variation1=1;
    double error_rate1 = 0.001;

    // Random number generator
    const gsl_rng_type * type;
    gsl_rng * r;
    gsl_rng_env_setup();
    type = gsl_rng_default;
    r = gsl_rng_alloc (type);
    gsl_rng_env_setup();

    unsigned long n_window;
    unsigned long div;

    unsigned long *vec_rd;
    double *vec_s, *vec_p, *vec_h, *vec0_s, *vec0_d, *vec0_h;
    double *covmat;

    struct combinatorial comb1, comb2;
    struct combinatorial_fst combfst;
    struct tests test1, test2;
    struct fst_calc fst;

    FILE * list_snps;
    unsigned long pos_snp;

    FILE * gff;
    char if_gff;
    unsigned long cds_start, cds_end;
    int phase_cds, frame;
    char strand, feature[256], *gff_file_name;
    unsigned long psyn, pnon, dsyn, dnon;



    /* Main part of the program */

    /* Read arguments */
    window_size=0;
    n01=0;
    n02=0;
    min_cov=4;
    max_cov=100;
    min_qual=10;
    min_mqual=10;
    m_bar=1;

    char from_stdin, outgroup_available, compute_fst, ext_snps, map_qual_pileup;
    char *pileup_file_name, *pileup2_file_name, *outgroup_file_name;
    char *snp_file_name;
    int arg_i, col_chrom, col_pos;
    from_stdin=2;
    outgroup_available=0;
    compute_fst=0;
    ext_snps=0;
    map_qual_pileup=0;
    if_gff=0;

    for(arg_i=1;arg_i<argc-1;arg_i++)
    {
        if (strcmp(argv[arg_i], "-n") == 0) {
            arg_i++;
            sscanf(argv[arg_i], "%u", &n01);
        } else if (strcmp(argv[arg_i], "-l") == 0) {
            arg_i++;
            sscanf(argv[arg_i], "%lu", &window_size);
        } else if (strcmp(argv[arg_i], "-cov") == 0) {
            arg_i++;
            sscanf(argv[arg_i], "%lf", &input_coverage1);
        } else if (strcmp(argv[arg_i], "-mincov") == 0) {
            arg_i++;
            sscanf(argv[arg_i], "%lu", &min_cov);
        } else if (strcmp(argv[arg_i], "-maxcov") == 0) {
            arg_i++;
            sscanf(argv[arg_i], "%lu", &max_cov);
        } else if (strcmp(argv[arg_i], "-minqual") == 0) {
            arg_i++;
            sscanf(argv[arg_i], "%lu", &min_qual);
        } else if (strcmp(argv[arg_i], "-nolowfreq") == 0) {
            arg_i++;
            sscanf(argv[arg_i], "%u", &m_bar);
        } else if (strcmp(argv[arg_i], "-outgroup") == 0) {
            outgroup_available=1; arg_i++;
            outgroup_file_name=argv[arg_i];
        } else if (strcmp(argv[arg_i], "-snpfile") == 0) {
            ext_snps=1;
            arg_i++;
            snp_file_name=argv[arg_i];
        } else if (strcmp(argv[arg_i], "-annot") == 0) {
            if_gff=1;
            arg_i++;
            gff_file_name=argv[arg_i];
        }
    }
    if (arg_i==argc-1) {
        if (strcmp(argv[arg_i], "-") == 0) {
            from_stdin=1;
        } else {
            from_stdin=0;
            pileup_file_name=argv[arg_i];
        }
    };

    if ((n01==0)||(window_size==0)||(from_stdin==2))
    {
        fprintf(stderr,"Missing values in command line!\n");
        fprintf(stderr,"  Command:\n");
        fprintf(stderr,"    npstat [options] file.pileup\n");
        fprintf(stderr,"  or to read from standard input:\n");
        fprintf(stderr,"    npstat [options] -\n");
        fprintf(stderr,"  Options:\n   -n samplesize : haploid sample size\n");
        fprintf(stderr,"   -l windowlength : window length\n");
        fprintf(stderr,"   -mincov minimum_coverage : filter on minimum coverage (default 4)\n");
        fprintf(stderr,"   -maxcov maximum_coverage : filter on maximum coverage (default 100)\n");
        fprintf(stderr,"   -minqual minimum_base_quality : filter on base quality (default 10)\n");
        fprintf(stderr,"   -nolowfreq m : filter on minimum allele count mac>m\n");
        fprintf(stderr,"   -outgroup file.fa : outgroup file in FASTA\n");
        fprintf(stderr,"   -annot file.gff3 : annotation file in GFF3\n");
        fprintf(stderr,"   -snpfile file.snp : consider SNPs only if present in file.snp\n");
        return(-1);
    };


    /* Initialize and open file */

    temp_file_name1=malloc(strlen(pileup_file_name)+strlen(pileup_file_name)+10);
    temp_file_name2=malloc(strlen(pileup_file_name)+strlen(pileup_file_name)+10);

    if (from_stdin==0) {
        bam_file1=fopen(pileup_file_name,"r");
        if (bam_file1 == NULL){
            fprintf(stderr,"Error: the pileup file cannot be opened!");
            return(-1);
        };
    } else {
        bam_file1=stdin;
        pileup_file_name="NPStat_input";
    };

    if(ext_snps==1){
        list_snps=fopen(snp_file_name,"r");
        if (list_snps == NULL){
            fprintf(stderr,"Error: the SNP list file cannot be opened!");
            return(-1);
        };
    };

    if(if_gff==1){
        gff=fopen(gff_file_name,"r");
        if (gff == NULL){
            fprintf(stderr,"Error: the GFF3 file cannot be opened!");
            return(-1);
        };
    };

    strcpy(temp_file_name1,pileup_file_name);
    strcpy(temp_file_name2,".stats");
    strcat(temp_file_name1,temp_file_name2);
    output_stat1=fopen(temp_file_name1,"w");
    if (output_stat1 == NULL){
        fprintf(stderr,"Error: the output file cannot be written!");
        return(-1);
    };

    //SNPS
    fasta_length=0;
    if (outgroup_available==1)
    { fasta_length=0;
        fasta_out=fopen(outgroup_file_name,"r");
        for (;iscntrl(fgetc(fasta_out))==0;)
        {
        };
        for (count_i=0;iscntrl(fgetc(fasta_out))==0;count_i++)
        {
        };
        fasta_length=count_i;
        rewind(fasta_out);
        for (;iscntrl(fgetc(fasta_out))==0;)
        {
        };
    };

    /* Load combinatorics */
    printf("Initializing combinatorics...\n");

    vec_rd=(unsigned long *)malloc(max_cov*sizeof(unsigned long));
    vec_s=(double *)malloc(max_cov*(n01-1)*sizeof(double));
    vec_p=(double *)malloc(max_cov*(n01-1)*sizeof(double));
    vec_h=(double *)malloc(max_cov*(n01-1)*sizeof(double));
    vec0_s=(double *)malloc(max_cov*sizeof(double));
    vec0_d=(double *)malloc(max_cov*sizeof(double));
    vec0_h=(double *)malloc(max_cov*sizeof(double));

    comb1.d_t=(double *)malloc(max_cov*sizeof(double));
    comb1.d_p=(double *)malloc(max_cov*sizeof(double));
    comb1.d_hl=(double *)malloc(max_cov*sizeof(double));
    comb1.d_hq=(double *)malloc(max_cov*sizeof(double));

    for(count_i=0;count_i<max_cov;count_i++)
    {
        comb1.d_t[count_i]=0;
        comb1.d_p[count_i]=0;
        comb1.d_hl[count_i]=0;
        comb1.d_hq[count_i]=0;
    };
    for(rd=2*m_bar+1;rd<=max_cov;rd++)
    {
        int i,j,k;
        comb1.d_p[rd-1]+=((double)n01-1)/(double)n01;
        comb1.d_hl[rd-1]+=(double)rd*(n01-1)/(double)(n01*(rd-1));
        comb1.d_hq[rd-1]+=(double)rd*((double)((n01-1)*(rd+1))/(double)(2*rd))/
            (double)(n01*(rd-1));
        for(k=1;k<n01;k++)
        {
            comb1.d_t[rd-1]+=(1-gsl_pow_int((double)k/(double)n01,rd)-
                gsl_pow_int(1-(double)k/(double)n01,rd))/(double)k;
            for(lt=1;lt<=m_bar;lt++)
            {
                comb1.d_t[rd-1]+=-gsl_sf_choose(rd,lt)*
                    gsl_pow_int((double)k/(double)n01,lt-1)*
                    gsl_pow_int(1-(double)k/(double)n01,rd-lt-1)/(double)n01;
                comb1.d_p[rd-1]+=-2*gsl_sf_choose(rd-2,lt-1)*
                    gsl_pow_int((double)k/(double)n01,lt-1)*
                    gsl_pow_int(1-(double)k/(double)n01,rd-lt-1)/(double)n01;
                comb1.d_hq[rd-1]+=-1/(double)(n01*(rd-1))*
                    ((double)lt*gsl_sf_choose(rd-1,lt-1)*
                    gsl_pow_int((double)k/(double)n01,lt-1)*
                    gsl_pow_int(1-(double)k/(double)n01,rd-lt) +
                    (double)(rd-lt)*gsl_sf_choose(rd-1,rd-lt-1)*
                    gsl_pow_int((double)k/(double)n01,rd-lt-1)*
                    gsl_pow_int(1-(double)k/(double)n01,lt));
            };
            comb1.d_hl[rd-1]+=((double)rd/(double)n01)*((1-2/((double)rd-1)) *
                (double)k/(double)n01-1) *
                gsl_pow_int((double)k/(double)n01,rd-2);
            comb1.d_hq[rd-1]+=-(double)rd/(double)(n01*(rd-1))*
                (gsl_pow_int((double)k/(double)n01,rd-1));
        };
    };

    for(rd=1;rd<=max_cov*(n01-1); rd++){
        vec_s[rd-1]=0;
        vec_p[rd-1]=0;
        vec_h[rd-1]=0;
    };

    for(rd=2*m_bar+1;rd<=max_cov; rd++){
        int k,covt;
        vec0_s[rd-1]=0;
        vec0_d[rd-1]=0;
        vec0_h[rd-1]=0;
        for(k=1;k<n01;k++){
            vec_s[(k-1)+(n01-1)*(rd-1)]=(1-gsl_pow_int((double)k/(double)n01,rd)-
                gsl_pow_int(1-(double)k/(double)n01,rd));
            vec_p[(k-1)+(n01-1)*(rd-1)]=(double)(2*k*(n01-k))/(double)(n01*n01);
            vec_h[(k-1)+(n01-1)*(rd-1)]=(double)(k*k)/(double)(n01*n01)+
                (double)(k)/(double)(n01*(rd-1))-(double)rd/(double)(rd-1)*
                gsl_pow_int((double)(k)/(double)(n01),rd);
            for(covt=1;covt<=m_bar;covt++){
                vec_s[(k-1)+(n01-1)*(rd-1)]+=-(
                    gsl_ran_binomial_pdf(covt,(double)k/(double)n01,rd)+
                    gsl_ran_binomial_pdf(rd-covt,(double)k/(double)n01,rd));
                vec_p[(k-1)+(n01-1)*(rd-1)]+=-(
                    2*(double)(covt*(rd-covt))/(double)(rd*(rd-1))*
                    (gsl_ran_binomial_pdf(covt,	(double)k/(double)n01,rd)+
                    gsl_ran_binomial_pdf(rd-covt,(double)k/(double)n01,rd)));
                vec_h[(k-1)+(n01-1)*(rd-1)]+=-(
                    (double)(covt*covt)/(double)(rd*(rd-1))*
                        gsl_ran_binomial_pdf(covt,(double)k/(double)n01,rd)+
                        (double)(rd-covt)*(double)(rd-covt)/(double)(rd*(rd-1))*
                        gsl_ran_binomial_pdf(rd-covt,(double)k/(double)n01,rd));
            };
            for(covt=m_bar+1;covt<=rd-m_bar-1;covt++){
                vec0_s[rd-1]+=gsl_ran_binomial_pdf(covt,(double)k/(double)n01,rd)/k;
                vec0_d[rd-1]+=gsl_pow_int(2*(double)(covt*(rd-covt))/
                    (double)(rd*(rd-1))/comb1.d_p[rd-1]-1/comb1.d_t[rd-1],2)*
                        gsl_ran_binomial_pdf(covt,(double)k/(double)n01,rd)/k;
                vec0_h[rd-1]+=gsl_pow_int((double)(covt*covt)/(double)(rd*(rd-1))/
                    comb1.d_hq[rd-1]-2*(double)(covt*(rd-covt))/(double)(rd*(rd-1))/
                        comb1.d_p[rd-1],2)*
                            gsl_ran_binomial_pdf(covt,(double)k/(double)n01,rd)/k;
            };
        };
    };

    covmat=(double *)malloc((n01-1)*(n01-1)*sizeof(double));
    generate_covariance_matrix(covmat,n01);


    /* Initialize variables */
    pos=0;
    oldpos=0;
    pos_base1=0;
    pos_base2=0;
    n_window=1;
    cds_start=0;
    cds_end=0;


    fprintf(output_stat1, "window\tlength\tlength_outgroup\tread_depth\tS\t");
    fprintf(output_stat1, "Watterson\tPi\tTajima_D\tvar_S\tvar_Watterson\t");
    fprintf(output_stat1, "unnorm_FayWu_H\tFayWu_H\tdiv\tnonsyn_pol\tsyn_pol\t");
    fprintf(output_stat1, "nonsyn_div\tsyn_div\talpha\n");

    printf("Computing statistics for the window:");
    /* Run across all bases */
    for(pos=1;(ct1!=EOF); pos++)
    {
        DEB(printf("new line\n"));

        if (pos==(n_window-1)*window_size+1)
        {
            // INITIALIZE EVERYTHING HERE!!!!!!!!!!!!!!

            printf(" %lu ", n_window);
            if ( n_window % 10 == 0 ) printf("\n");

            test1.cov=0;
            test1.l=0;
            test1.l_out=0;
            test1.s=0;
            test1.num_t=0;
            test1.num_p=0;
            test1.num_hl=0;
            test1.num_hq=0;
            test1.den_t=0;
            test1.den_p=0;
            test1.den_hl=0;
            test1.den_hq=0;
            div=0;
            for(rd=1;rd<=max_cov; rd++){
                vec_rd[rd-1]=0;
            };
            psyn=0;
            dsyn=0;
            pnon=0;
            dnon=0;

        };
        /* Read data */
        if(ext_snps==1){
            while (pos_snp<pos) pos_snp=extract_pos_snpinput(list_snps);
        };

        if(if_gff==1){
            while (cds_end<pos) {
                char *line_gff;
                size_t n_line_gff;
                n_line_gff=1;
                line_gff=malloc(sizeof(char));
                if(getline(&line_gff,&n_line_gff,gff)!=-1){
                    while ((line_gff[0]=='#')||(line_gff[0]=='\0')) {
                        getline(&line_gff,&n_line_gff,gff);
                    }; sscanf(line_gff,"%*s\t%*s\t%s\t%lu\t%lu\t%*s\t%c\t%u\t",
                              feature,&cds_start,&cds_end,&strand,&phase_cds);
                    if((strcmp(feature,"CDS")!=0)&&(strcmp(feature,"cds")!=0)&&
                       (strcmp(feature,"SO:0000316")!=0)) {
                        cds_end=0;
                    };
                } else {cds_end=-1;};
            };
        };

        if (pos_base1<pos)
        {
            DEB(printf("reading new base from file 1\n")); //debug
            read_line_pileup(bam_file1, min_qual, min_mqual, &pos_base1,
                             &n_ref1, &n_alt_allele1, &rd1, n_alt_1,
                             &ref_base1, &alt_base1);
        };

        /* Extract statistics */
        if ((pos_base1==pos)&&(rd1>=min_cov)&&(rd1<=max_cov))
        {
            if(outgroup_available==1){
                out_base=extract_outgroup_base(fasta_out,pos,oldpos,fasta_length); //1;
            } else {
                out_base=0;
            };
            if((ext_snps==1)&&(pos_snp!=pos)) {
                if (n_ref1>=n_alt_allele1){
                    n_ref1+=n_alt_allele1;
                    n_alt_allele1=0;
                } else {
                    n_alt_allele1+=n_ref1;
                    n_ref1=0;
                };
            };
            if (rd1>=max(min_cov,2*m_bar+2))
            {
                if (if_gff==1){
                    if (cds_start<=pos){
                        if(strand=='+'){
                            frame=((pos-cds_start+3-phase_cds)%3)+1;
                        } else {
                            if(strand=='-'){
                                frame=((cds_end+3-phase_cds-pos)%3)+1;
                            } else frame=0;
                        };
                    } else frame=0;
                };
                extract_stats(&test1, &comb1, n01, n_ref1, n_alt_allele1, rd1,
                              n_alt_1, ref_base1, alt_base1, out_base, m_bar);
                if(out_base!=0){
                    if (((n_alt_allele1>=rd1-m_bar)&&(alt_base1!=out_base))||
                        ((n_ref1>=rd1-m_bar)&&(ref_base1!=out_base))){
                        div++;
                        if (if_gff==1){
                            if (frame==3) dsyn++;
                            if ((frame==1)||(frame==2)) dnon++;
                        };
                    };
                    if ((if_gff==1)&&(n_alt_allele1>m_bar)&&(n_ref1>m_bar)){
                        if (frame==3) psyn++;
                        if ((frame==1)||(frame==2)) pnon++;
                    };
                };
                vec_rd[rd1-1]++;
            };


            oldpos=pos;

        };
        /* Print output */
        ct1=fgetc(bam_file1); ungetc(ct1,bam_file1);

        if ((pos==(n_window*window_size))||(ct1==EOF))
        {
            double theta1_val, pi1_val, d1_val, h1_val, theta2_val, pi2_val;
            double d2_val, h2_val, pia_val, fst_val, cov1_val, cov2_val;
            double div_val, var_h, var_d, var_s, var0_s, var0_d, var0_h;
            double vk_s[n01-1], vk_d[n01-1], vk_h[n01-1];
            int k;
            DEB(printf("printing output\n")); //debug

            var_s=0;
            var_d=0;
            var_h=0;
            var0_s=0;
            var0_d=0;
            var0_h=0;
            if ((test1.den_t>0)&&(test1.den_p>0)) {
                if (test1.den_hq>0) {
                    for(k=1;k<n01;k++){
                        vk_s[k-1]=0;
                        vk_d[k-1]=0;
                        vk_h[k-1]=0;
                        for(rd=max(min_cov,2*m_bar+2);rd<=max_cov; rd++){
                            vk_s[k-1]+=vec_rd[rd-1]*vec_s[(k-1)+(n01-1)*(rd-1)];
                            vk_d[k-1]+=vec_rd[rd-1]*(vec_p[(k-1)+(n01-1)*(rd-1)]/
                                test1.den_p-vec_s[(k-1)+(n01-1)*(rd-1)]/test1.den_t);
                            vk_h[k-1]+=vec_rd[rd-1]*(vec_h[(k-1)+(n01-1)*(rd-1)]/
                                test1.den_hq-vec_p[(k-1)+(n01-1)*(rd-1)]/test1.den_p);
                        };
                    };
                    for(rd=max(min_cov,2*m_bar+2);rd<=max_cov; rd++){
                        var0_s+=vec_rd[rd-1]*vec0_s[rd-1];
                        var0_d+=vec_rd[rd-1]*vec0_d[rd-1];
                        var0_h+=vec_rd[rd-1]*vec0_h[rd-1];
                    };
                    for(k=1;k<n01;k++){
                        for(lt=1;lt<n01;lt++){
                            var_s+=vk_s[k-1]*vk_s[lt-1]*covmat[(n01-1)*(lt-1)+k-1];
                            var_d+=vk_d[k-1]*vk_d[lt-1]*covmat[(n01-1)*(lt-1)+k-1];
                            var_h+=vk_h[k-1]*vk_h[lt-1]*covmat[(n01-1)*(lt-1)+k-1];
                        };
                    };
                    var0_d=var0_d/(double)(test1.l*test1.l);
                    var0_h=var0_h/(double)(test1.l_out*test1.l_out);
                } else {
                    for(k=1;k<n01;k++){
                        vk_s[k-1]=0;
                        vk_d[k-1]=0;
                        for(rd=max(min_cov,2*m_bar+2);rd<=max_cov; rd++){
                            vk_s[k-1]+=vec_rd[rd-1]*vec_s[(k-1)+(n01-1)*(rd-1)];
                            vk_d[k-1]+=vec_rd[rd-1]*(vec_p[(k-1)+(n01-1)*(rd-1)]/
                                test1.den_p-vec_s[(k-1)+(n01-1)*(rd-1)]/test1.den_t);
                        };
                    };
                    for(rd=max(min_cov,2*m_bar+2);rd<=max_cov; rd++){
                        var0_s+=vec_rd[rd-1]*vec0_s[rd-1];
                        var0_d+=vec_rd[rd-1]*vec0_d[rd-1];
                    };
                    for(k=1;k<n01;k++){
                        for(lt=1;lt<n01;lt++){
                            var_s+=vk_s[k-1]*vk_s[lt-1]*covmat[(n01-1)*(lt-1)+k-1];
                            var_d+=vk_d[k-1]*vk_d[lt-1]*covmat[(n01-1)*(lt-1)+k-1];
                        };
                    };
                    var0_d=var0_d/(double)(test1.l*test1.l);
                };
            };

            if(test1.l>0) {
                cov1_val=(double)(test1.cov)/(double)(test1.l);
            } else { cov1_val=-1; };
            if(test1.den_t>0) {
                theta1_val=test1.num_t/test1.den_t;
            } else { theta1_val=-1; };
            if(test1.den_p>0) {
                pi1_val=test1.num_p/test1.den_p;
            } else { pi1_val=-1; };
            if((test1.den_t>0)&&(test1.den_p>0)) {
                d1_val=pi1_val-theta1_val;
            } else { d1_val=-1; };
            if((test1.den_p>0)&&(test1.den_hq>0)) {
                h1_val=test1.num_hq/test1.den_hq-pi1_val;
            } else { h1_val=-1; };
            if(test1.l_out>0) {
                div_val=(double)(div)/(double)(test1.l_out);
            } else { div_val=-1; };

            var0_s=var0_s*theta1_val;
            var0_d=var0_d*theta1_val;
            var0_h=var0_h*theta1_val;
            var_s=var_s*theta1_val*theta1_val;
            var_d=var_d*theta1_val*theta1_val;
            var_h=var_h*theta1_val*theta1_val;

            fprintf(output_stat1, "%lu\t%lu\t%lu",n_window, test1.l, test1.l_out);
            if (test1.l>0) {
                fprintf(output_stat1, "\t%f\t%lu\t%f\t%f\t%f\t%f\t%f", cov1_val,
                        test1.s, theta1_val, pi1_val, d1_val/sqrt(var0_d+var_d),
                        var0_s+var_s, (var0_s+var_s)/(test1.den_t*test1.den_t));
            } else {
                fprintf(output_stat1, "\tNA\t0\tNA\tNA\tNA\tNA\tNA");
            };
            if (test1.l_out>0) {
                fprintf(output_stat1, "\t%f\t%f\t%f", h1_val,
                        h1_val/sqrt(var0_h+var_h), div_val);
            } else {
                fprintf(output_stat1, "\tNA\tNA\tNA");
            };
            if (if_gff==1) {
                fprintf(output_stat1, "\t%lu\t%lu\t%lu\t%lu",
                        pnon, psyn, dnon, dsyn);
                if(dsyn*pnon==0){
                    if(psyn*dnon==0){
                        fprintf(output_stat1, "\tNA");
                    } else {
                        fprintf(output_stat1, "\tInf");
                    };
                } else {
                    fprintf(output_stat1, "\t%f",
                            1-(double)(dsyn*pnon)/(double)(psyn*dnon));
                };
            } else {
                fprintf(output_stat1, "\tNA\tNA\tNA\tNA\tNA");
            };
            DEB(fprintf(output_stat1, "\tvars\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",
                        var0_s,var_s,var0_s/(test1.den_t*test1.den_t+0.000001),
                        var_s/(test1.den_t*test1.den_t+0.000001),
                        var0_d/(theta1_val),var_d/(theta1_val*theta1_val),
                        var0_h/(theta1_val),var_h/(theta1_val*theta1_val));)
                fprintf(output_stat1, "\n");

            n_window++;
        };

    };
    /* Close files, free gsl random number generator */
    fclose(bam_file1);
    if (outgroup_available==1) { fclose(fasta_out); };
    if (ext_snps==1) { fclose(list_snps); };
    if (if_gff==1) { fclose(gff); };
    fclose(output_stat1);
    gsl_rng_free(r);
    printf("\n");
    printf("Computation of statistics completed.\n");
}
