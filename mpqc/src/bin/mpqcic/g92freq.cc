#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <strings.h>

#include <util/keyval/keyval.h>
#include <math/scmat/matrix.h>
#include <math/scmat/local.h>
#include <chemistry/molecule/molecule.h>

#include "g92_int.h"

int
g92_freq_driver(char *name_in, RefMolecule& mole,
                const RefKeyVal &g92_keyval,
                char *method, double &energy,
                RefSCVector& gradient, RefSCVector &frequencies,
                RefSCVector &normalmodes, int& nmodes, int& nimag)
{
    // Read necessary pieces from KeyVal
    // Read g92_method iff methd parameter not set
    char *g92_method;
    if (!method)
        if (g92_keyval->exists("g92_method"))
            g92_method=g92_keyval->pcharvalue("g92_method");
        else
        {
            fprintf(stderr,"Could not find g92_method in input\n");
            return -1;
        }
    else
        g92_method=method;

    int  charge=0;
    if (g92_keyval->exists("charge"))
        charge=g92_keyval->intvalue("charge");
    else
    {
        //make sure that this molecule has a closed shell neutral configuration
        int net_charge=0;
        for (int i=0; i<mole->natom(); i++)
            net_charge+=mole->atom(i).element().number();
        if (net_charge%2 !=0)
        {
            fprintf(stderr," Neutral closed shell molecule not possible\n");
            fprintf(stderr," You must add a g92:charge keyval\n");
            return -1;
        }
    }
    int  multiplicity=1;
    if (g92_keyval->exists("multiplicity"))
        multiplicity=g92_keyval->intvalue("multiplicity");
    int  memory=4000000;
    if (g92_keyval->exists("memory"))
        memory=g92_keyval->intvalue("memory");
    int  keep_g92_log=1;
    if (g92_keyval->exists("keep_g92_log"))
        keep_g92_log=g92_keyval->intvalue("keep_g92_log");

    // Directories needed for setting environmental variables
    char *scratch_dir=NULL;
    if (g92_keyval->exists("scratch_dir"))
        scratch_dir=g92_keyval->pcharvalue("scratch_dir");
    char *g92_dir=NULL;
    if (g92_keyval->exists("g92_dir"))
        g92_dir=g92_keyval->pcharvalue("g92_dir");

    // Make sure we have a legitimate run name
    char *name;
    if (!name_in)
    {
        name = new char[8];
        strcpy(name,"g92_tmp");
    }
    else
    {
        name = new char[strlen(name_in)+1];
        strcpy(name, name_in);
    }

    // Convert the runtype name to an integer
    int runtype=0;
    while (runtype < n_g92_calc_types &&
           strcasecmp(g92_method, g92_calc[runtype].name)) runtype++;

    if (runtype == n_g92_calc_types)   
    {
        fprintf(stderr,"Run_g92 does not recognize the calculation type: %s\n",
                g92_method);
        return -1;
    }

    // if this runtype needs a basis set, get it from the keyval
    // and check if a checkpoint file can be used for a guess
    char *basis=NULL;
    int use_checkpoint_guess=0;
    if (g92_calc[runtype].basis_req)
    {
        if (g92_keyval->exists("basis"))
            basis=g92_keyval->pcharvalue("basis");
        else
        {
            fprintf(stderr,"G92 method requires basis set to be specified in"
                    " input\n");
            return -1;
        }            
        if (g92_keyval->exists("use_checkpoint_guess"))
            use_checkpoint_guess=g92_keyval->intvalue("use_checkpoint_guess");
    }
        
    // Run g92 calculations and then parse the output 
    run_g92_freq(name, runtype, basis, memory,
                 use_checkpoint_guess, scratch_dir, g92_dir, mole,
                 charge, multiplicity);

    if (parse_g92_freq(name, g92_calc[runtype].parse_string, mole,
                       energy, gradient, frequencies, normalmodes,
                       nmodes, nimag))
    {
        fprintf(stderr,"Error parsing G92 Force calculation\n");
        fprintf(stderr,"Check output in %s.out\n",name);
        return -1;
    }

    // Clean up output files
    char commandstr[100];
    if (keep_g92_log)
    {
        sprintf(commandstr,"cat %s.g92freq.out >> %s.freqlog",name,name);
        system(commandstr);
    }
    else if (!use_checkpoint_guess)
    {
        sprintf(commandstr,"%s.chk",name);
        unlink(commandstr);
    }
    sprintf(commandstr,"%s.g92freq.out",name);
    unlink(commandstr);

    return 0;
}    

int
run_g92_freq(char *prefix, int runtype, char *basis, int memory,
             int chk_guess,char *scratch, char *g92_dir, RefMolecule &mole,
             int charge, int multiplicity)
{
    FILE *fp_g92_input;
    char *infilename=(char*) malloc(strlen(prefix)+5);
    char *outfilename=(char*) malloc(strlen(prefix)+13);
    char commandstr[1000];
    
    strcat(strcpy(infilename,prefix),".com");
    strcat(strcpy(outfilename,prefix),".g92freq.out");
    //printf("G92 input filename: %s\n",infilename);  
    //printf("G92 output filename: %s\n",outfilename);
    
    // First assemble the input file 
    fp_g92_input=fopen(infilename, "w");
    
    // Write out required headers 
    fprintf(fp_g92_input,"%%chk=%s\n",prefix);
    fprintf(fp_g92_input,"%%mem=%d\n",memory);
    fprintf(fp_g92_input,"#p units=au Freq %s", g92_calc[runtype].command); 
    if (basis)
        fprintf(fp_g92_input,"/%s\n",basis);
    else
        fprintf(fp_g92_input,"\n");

    // Request to use guess from checkpoint file
    if (chk_guess)
    {
        char *tmp=(char*)malloc(sizeof(char)*(5+strlen(prefix)));
        sprintf(tmp,"%s.chk",prefix);
        if (FILE* fp_chk=fopen(tmp,"r"))
        {
            fprintf(fp_g92_input,"guess=check\n");
            fclose(fp_chk);
        }
        free (tmp);
    }

    fprintf(fp_g92_input,"\nG92 input generated by MPQC opt\n");
    fprintf(fp_g92_input,"\n %d %d\n",charge, multiplicity);
    for (int i=0; i<mole->natom(); i++)
        fprintf(fp_g92_input,"%s %lf %lf %lf\n",
                mole->atom(i).element().symbol(),
                mole->atom(i).operator[](0), mole->atom(i).operator[](1),
                mole->atom(i).operator[](2));
    fprintf(fp_g92_input,"\n");
    fclose(fp_g92_input);

    // Set environmental variable necessary for g92 run
    // Do this only once, or else the execution shell loses it
    static env_var_set=0;
    if (!env_var_set)
    {
        env_var_set=1;
        static char *g92_env_str;
        if (g92_dir)
        {
            g92_env_str=(char*) malloc(sizeof(char)*(14+strlen(g92_dir)));
            sprintf(g92_env_str,"GAUSS_EXEDIR=%s",g92_dir);
            putenv(g92_env_str);
        }
        static char *g92_scr_str;
        if (scratch)
        {
            g92_scr_str=(char*) malloc(sizeof(char)*(16+strlen(scratch)));
            sprintf(g92_scr_str,"GAUSS_SCRDIR=%s",scratch);
        }
        else
            g92_scr_str="GAUSS_SCRDIR=./";
        putenv(g92_scr_str);
    }

    // assemble and execute g92 command 
    sprintf(commandstr,"%s/g92 < %s > %s",g92_dir,infilename, outfilename);
    int ret = system(commandstr);
    
    // Free filesnames
    free(outfilename);
    free(infilename);

    return ret;
}

int
parse_g92_freq(char *prefix,char *parse_string, RefMolecule mole,
               double &energy, RefSCVector& gradient,
               RefSCVector& frequencies, RefSCVector &normalmodes,
               int& nmodes, int& nimag)
{
    /* read and parse output file */
    FILE *fp_g92_output;
    char *outfilename=(char*) malloc(strlen(prefix)+13);
    char *word;
    char line[120];
    
    int natoms=mole->natom();
    strcat(strcpy(outfilename,prefix),".g92freq.out");
    fp_g92_output=fopen(outfilename,"r");
    
    if (!fp_g92_output)
    {
        fprintf(stderr,"Trouble in parse_g92, no output file:\n",
                outfilename);
        return -1;
    }
    
    /* Read through the output file, finding the nuclear repulsion
        energy, converged scf energy, and archive entry */
    while(fgets(line,120,fp_g92_output))
    {
        if (strstr(line,
                   "Harmonic frequencies (cm**-1), IR intensities (KM/Mole),"
                   ))
            goto found_freq;
        if (strstr(line,
                   "imaginary frequencies (negative signs)"
                   ))
        {               
            word = strtok(line,"* ");
            nimag=strtol(word,NULL,10);
        }
    }
    fprintf(stderr,"Could not find forces in G92 scf output\n");
    return -1;
    
  found_freq:
    // Skip next two lines
    fgets(line,120,fp_g92_output); fgets(line,120,fp_g92_output);

    // Figure out how many degrees of freedom
    nmodes=3*natoms-6;
    if (natoms==2) nmodes=1;
    int nsets=nmodes/5;
    if (nmodes%5) nsets++;
    int count1=0, count2=0, total_count=0;

    // Now read in natoms lines of forces
    for (int i=0; i<nsets; i++)
    {
        int nrows=(i==nsets-1 && (nmodes%5))?nmodes%5:5;

        fgets(line,120,fp_g92_output);
        word = strtok(line," ");
        for (int j=0;j<nrows;j++)
        {
            count1++;
            if (strtol(word, NULL, 10) != count1)
            {
                fprintf(stderr,"Error parsing freq, count=%d\n",count1);
                fprintf(stderr,"Error parsing freq, look at %s\n",outfilename);
                return -1;
            }
            word=strtok(NULL," ");
        }
        fgets(line,120,fp_g92_output); fgets(line,120,fp_g92_output);
        word = strtok(line+22," ");
        for (j=0;j<nrows;j++)
        {
            frequencies.set_element(count2,strtod(word,NULL));
            word=strtok(NULL," ");
            count2++;
        }
        fgets(line,120,fp_g92_output);fgets(line,120,fp_g92_output);
        fgets(line,120,fp_g92_output);fgets(line,120,fp_g92_output);
        fgets(line,120,fp_g92_output);fgets(line,120,fp_g92_output);
        for (int atoms=0; atoms<natoms; atoms++)
        {
            for (int coords=0; coords<3; coords++)
            {
                fgets(line,120,fp_g92_output);

                if (strtol(strtok(line," "),NULL,10)!=coords+1)
                {
                    fprintf(stderr," Error parsing G92 Freq output\n");
                    return -1;
                }
                if (strtol(strtok(NULL," "),NULL,10)!=atoms+1)
                {
                    fprintf(stderr," Error parsing G92 Freq output\n");
                    return -1;
                }
                if (strtol(strtok(NULL," "),NULL,10)!=
                    mole->atom(atoms).element().number())
                {
                    fprintf(stderr," Error parsing G92 Freq output\n");
                    return -1;
                }
                
                for (int j=0;j<nrows;j++)
                {
                    int nmoffset=((i*5+j)*natoms*3+atoms*3+coords);
                    normalmodes.set_element(nmoffset,
                                            strtod(strtok(NULL," "),NULL));
                    total_count++;
                }
            }
        }
    }
    
    printf("Nmodes=%d\n",nmodes);
    // Now let's look for the archive entry
    while(fgets(line,120,fp_g92_output))
    {
        if (strstr(line,"1\\1\\"))
            goto found_archive;
    }
    fprintf(stderr,"Could not find archive entry in G92 scf output\n");
    fprintf(stderr,"Look at %s\n",outfilename);
    return -1;
    
    // Finally, trim the archive entry out 
    
  found_archive:
    int archive_size=10000;
    int arch_ptr=0;
    int done=0;
    char *archive=(char *) malloc(archive_size);
    
    while (!done)
    {
        if (strstr(line,"@")) done=1;
        
        word=strtok(line," \n");
        while (word)
        {
            strcpy(archive+arch_ptr,word);
            arch_ptr+=strlen(word);
            if (arch_ptr+120 > archive_size)
            {
                archive_size *= 2;
                archive=(char*)realloc(archive, archive_size);
            }
            word=strtok(NULL," \n");
        }
        if (!fgets(line,120,fp_g92_output))
        {
            fprintf(stderr,"Premature end of archive entry\n");
            fprintf(stderr,"Look at %s\n",outfilename);
            return -1;
        }
    }

    // Now parse out the appropriate energy term
    char *tmp = strstr(archive, parse_string);
    char *tmp2=strchr(tmp,'=');
    if (!tmp2)
    {
        fprintf(stderr,"Error parsing energy from archive entry\n");
        fprintf(stderr,"Look at %s\n",outfilename);
        return -1;
    }        
    energy=strtod(strtok(tmp2,"=\\"),NULL);

    // Free archive memory
    free(archive);
    
    // Close and Free filesnames 
    fclose(fp_g92_output);
    free(outfilename);

    return 0;
}

#ifdef TEST_MAIN
// Simple main to test run_g92
main()
{
    char *name="1mbi";
    int natoms=19;
    char *atoms[]={"C","C","N","C","N","C","C","H","C","C",
                   "N","H","H","C","H","H","H","H","H"};
    int charge=1;
    int multiplicity=1;
    double energy;
    LocalSCDimension dim(19);
    LocalSCDimension dimfreq(51);
    LocalSCDimension dim2(2907);
    RefSCVector gradient = new LocalSCVector(&dim);
    RefSCVector geom = new LocalSCVector(&dim);
    RefSCVector freq = new LocalSCVector(&dimfreq);
    RefSCVector normalmodes = new LocalSCVector(&dim2);
    
    AtomicCenter a1("C"  , -0.3605513998 ,  -0.0367641686  ,  0.0147926861 ,NULL);
    AtomicCenter a2("C"  , -0.3147667922 ,  -0.0218332569  ,  2.8476409608 ,NULL);
    AtomicCenter a3("N"  ,  1.8832840718 ,  -0.0339070724  , -0.9335183676 ,NULL);
    AtomicCenter a4("C"  , -2.7589829999 ,   0.0431473081  , -1.2585755617 ,NULL);
    AtomicCenter a5("N"  ,  2.0297956234 ,  -0.0001455266  ,  3.4807891018 ,NULL);
    AtomicCenter a6("C"  ,  3.4008141314 ,  -0.3069073243  ,  1.1156665213 ,NULL);
    AtomicCenter a7("C"  , -2.5959947997 ,   0.0911156949  ,  4.3063127063 ,NULL);
    AtomicCenter a8("H"  , -2.8188912524 ,   0.0238432165  , -3.2842659790 ,NULL);
    AtomicCenter a9("C"  , -4.8270216002 ,   0.1379091399  ,  0.1941044063 ,NULL);
    AtomicCenter a10("C" ,  3.5647013727 ,  -0.2221123911  ,  5.6805410998 ,NULL);
    AtomicCenter a11("N" ,  5.6135737940 ,  -1.0591150798  ,  1.4662313919 ,NULL);
    AtomicCenter a12("H" , -2.5308437439 ,   0.1269304909  ,  6.3315010948 ,NULL);
    AtomicCenter a13("H" , -6.6619683656 ,   0.2095229711  , -0.6724861034 ,NULL);
    AtomicCenter a14("C" , -4.7496587247 ,   0.1704403688  ,  2.9840462310 ,NULL);
    AtomicCenter a15("H" ,  2.5775316127 ,  -1.1076734509  ,  7.2210660123 ,NULL);
    AtomicCenter a16("H" ,  4.3963259769 ,   1.5874896722  ,  6.1462787592 ,NULL);
    AtomicCenter a17("H" ,  6.8507089561 ,  -1.3574062121  ,  0.0658293996 ,NULL);
    AtomicCenter a18("H" ,  5.7350153987 ,  -1.1717726768  ,  3.6404193236 ,NULL);
    AtomicCenter a19("H" , -6.5282305620 ,   0.2786975609  ,  3.9567978488 ,NULL);

    // Set up molecule
    RefMolecule mole = new Molecule();

    mole->add_atom(0,a1);
    mole->add_atom(1,a2);
    mole->add_atom(2,a3);
    mole->add_atom(3,a4);
    mole->add_atom(4,a5);
    mole->add_atom(5,a6);
    mole->add_atom(6,a7);
    mole->add_atom(7,a8);
    mole->add_atom(8,a9);
    mole->add_atom(9,a10);
    mole->add_atom(10,a11);
    mole->add_atom(11,a12);
    mole->add_atom(12,a13);
    mole->add_atom(13,a14);
    mole->add_atom(14,a15);
    mole->add_atom(15,a16);
    mole->add_atom(16,a17);
    mole->add_atom(17,a18);
    mole->add_atom(18,a19);
    
    // Set up keyval
    char *filename="run_g92.in";
    const RefKeyVal refg92_keyval(new ParsedKeyVal(filename));
    
    int nimag, nmodes;
    char *method=NULL;

    for (int j=0; j<1; j++)
    {
        int retval=g92_freq_driver(name, mole, refg92_keyval, method, energy,
                                   gradient, freq, normalmodes, nmodes, nimag);
        printf("Run #%d E=%lf\n",j,energy);
    }
    printf(" Nimag=%d, nmodes=%d\n",nimag, nmodes);
    for (int i=0; i<nmodes; i++)
    {
        printf("%d %15.10lf \n",i,freq.get_element(i));
        for (int j=0; j<natoms*3; j++)
            printf("Component  #%d = %lf\n",j,normalmodes.get_element(i*natoms*3+j));
    }
    fflush(stdout);
}
#endif // TEST_MAIN
