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
run_g92(char *name_in, const RefKeyVal& g92_keyval, RefMolecule& mole,
        double& energy, RefSCVector& gradient)
{
    // Read necessary pieces from KeyVal
    char *g92_method;
    if (g92_keyval->exists("g92_method"))
        g92_method=g92_keyval->pcharvalue("g92_method");
    else
    {
        fprintf(stderr,"Could not find g92_method in input\n");
        return -1;
    }
    int  charge=0;
    
    charge=g92_keyval->intvalue("charge");
    
    if (g92_keyval->error() != KeyVal::OK &&
        !g92_keyval->exists("multiplicity"))
    {
        // make sure that this molecule has a closed shell neutral
        // configuration
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
    run_g92_calc(name, runtype, basis, memory,
                 use_checkpoint_guess, scratch_dir, g92_dir, mole,
                 charge, multiplicity);

    if (parse_g92(name, g92_calc[runtype].parse_string, mole->natom(),
                  energy, gradient))
    {
        fprintf(stderr,"Error parsing G92 Force calculation\n");
        fprintf(stderr,"Check output in %s.out\n",name);
        return -1;
    }

    // Clean up output files
    char commandstr[100];
    if (keep_g92_log)
    {
        sprintf(commandstr,"cat %s.g92.out >> %s.log",name,name);
        system(commandstr);
    }
    else if (!use_checkpoint_guess)
    {
        sprintf(commandstr,"%s.chk",name);
        unlink(commandstr);
    }
    sprintf(commandstr,"%s.g92.out",name);
    unlink(commandstr);

    return 0;
}    

int
run_g92_calc(char *prefix, int runtype, char *basis, int memory,
             int chk_guess, char *scratch, char *g92_dir, RefMolecule &mole,
             int charge, int multiplicity)
{
    FILE *fp_g92_input;
    char *infilename=(char*) malloc(strlen(prefix)+5);
    char *outfilename=(char*) malloc(strlen(prefix)+9);
    char commandstr[1000];
    
    strcat(strcpy(infilename,prefix),".com");
    strcat(strcpy(outfilename,prefix),".g92.out");
    //printf("G92 input filename: %s\n",infilename);  
    //printf("G92 output filename: %s\n",outfilename);
    
    // First assemble the input file 
    fp_g92_input=fopen(infilename, "w");
    
    // Write out required headers 
    fprintf(fp_g92_input,"%%chk=%s\n",prefix);
    fprintf(fp_g92_input,"%%mem=%d\n",memory);
    fprintf(fp_g92_input,"#p units=au Force %s", g92_calc[runtype].command); 
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

    // assemble and execute scf command 
    sprintf(commandstr,"%s/g92 < %s > %s",g92_dir,infilename, outfilename);
    int ret = system(commandstr);
    
    // Free filesnames
    free(outfilename);
    free(infilename);

    return ret;
}

int
parse_g92(char *prefix, char *parse_string, int natoms, double & energy,
          RefSCVector& gradient)
{
    /* read and parse output file */
    FILE *fp_g92_output;
    char *outfilename=(char*) malloc(strlen(prefix)+9);
    char *word;
    char line[120];
    
    strcat(strcpy(outfilename,prefix),".g92.out");
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
                   "Center     Atomic                   Forces (Hartrees/Bohr)"
                   ))
            goto found_forces;
    }
    fprintf(stderr,"Could not find forces in G92 scf output\n");
    return -1;
    
  found_forces:
    // Skip next two lines
    fgets(line,120,fp_g92_output); fgets(line,120,fp_g92_output);
    
    // Now read in natoms lines of forces
    for (int i=0; i<natoms; i++)
    {
        fgets(line,120,fp_g92_output);
        word = strtok(line," ");
        if (strtol(word, NULL, 10) != i+1)
        {
            fprintf(stderr,"Error parsing forces, look at %s\n",outfilename);
            return -1;
        }
        word=(strtok(NULL," "),strtok(NULL," "));
        gradient.set_element(i*3,strtod(word,NULL));
        word=strtok(NULL," ");
        gradient.set_element(i*3+1,strtod(word,NULL));
        word=strtok(NULL," ");
        gradient.set_element(i*3+2,strtod(word,NULL));
    }
    
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
int
main()
{
    char *name="water";
    int natoms=3;
    char *atoms[]={"O","H","H"};
    int charge=0;
    int multiplicity=1;
    double energy;
    LocalSCDimension dim(9);
    RefSCVector gradient = new LocalSCVector(&dim);
    RefSCVector geom = new LocalSCVector(&dim);
    
    AtomicCenter O("O",0.,0.,0.,NULL);
    AtomicCenter H1("H",0.,0.,1.79018,NULL);
    AtomicCenter H2("H",1.72508,0.,-0.47838,NULL);

    // Set up molecule
    RefMolecule mole = new Molecule();

    mole->add_atom(0,O);
    mole->add_atom(1,H1);
    mole->add_atom(2,H2);
    
    // Set up keyval
    char *filename="run_g92.in";
    ParsedKeyVal g92_keyval(filename);

    int memory;
    char *basis;
    for (int j=0; j<10; j++)
    {
        int retval=run_g92(name, mole, g92_keyval, energy, gradient);
        printf("Run #%d E=%lf\n",j,energy);
    }
    printf(" Energy = %lf\n",energy);
    for (int i=0; i<natoms; i++)
        printf("%d   %15.10lf   %15.10lf   %15.10lf\n",i,gradient.get_element(i*3),
               gradient.get_element(i*3+1),gradient.get_element(i*3+2));

}
#endif // TEST_MAIN



