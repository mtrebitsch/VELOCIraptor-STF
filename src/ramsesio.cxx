/*! \file ramsesio.cxx
 *  \brief this file contains routines for ramses snapshot file io
 *
 * \todo need to check if amr file quantity ngrid_current is actually the number of cells in the file as
 * an example fortran code I have been given seems to also use the ngridlevel array, which stores the number of cells
 * at a given resolution level.
 * \todo need to add in ability for multiple read threads and sends between read threads
 *
 *
 * Edited by:    Rodrigo Ca\~nas
 *               rodrigo.canas@icrar.org
 *
 * Last edited:  7 - Jun - 2017
 *
 *
 */

/* TODO MT (not in order):
   - MPI with baryons [this has something to do with mpi_nsend_readthread_baryon not being initialized]
   - Read BH  -> compiles, seems to work
   - convert stellar birth_time in ages (if possible?)
   - zoom with refinement scalar (only count gas cells above a given threshold?)
   - check what ninputoffset and inputoffsetinfile are

   NOTES:
   - the ids are not unique: a star and a DM particle can share an ID... -> in progress for particles
   - Ref scalar: either read in get_nbodies (ugly) or resize after ReadData
*/

//-- RAMSES SPECIFIC IO

#include <sstream>

#include "stf.h"

#include "ramsesitems.h"
#include "endianutils.h"


int RAMSES_fortran_read(fstream &F, int &i){
    int dummy,byteoffset=0;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    F.read((char*)&i,sizeof(int)); byteoffset+=dummy;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    return byteoffset;
}
int RAMSES_fortran_read(fstream &F, RAMSESFLOAT &f){
    int dummy,byteoffset=0;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    F.read((char*)&f,sizeof(RAMSESFLOAT)); byteoffset+=dummy;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    return byteoffset;
}
int RAMSES_fortran_read(fstream &F, int *i){
    int dummy,byteoffset=0;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    F.read((char*)i,dummy); byteoffset+=dummy;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    return byteoffset;
}
int RAMSES_fortran_read(fstream &F, unsigned int *i){
    int dummy,byteoffset=0;
    F.read((char*)&dummy, sizeof(dummy)); byteoffset += sizeof(int);
    F.read((char*)i,dummy);               byteoffset += dummy;
    F.read((char*)&dummy, sizeof(dummy)); byteoffset += sizeof(int);
    return byteoffset;
}
int RAMSES_fortran_read(fstream &F, long long *i){
    int dummy,byteoffset=0;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    F.read((char*)i,dummy); byteoffset+=dummy;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    return byteoffset;
}
int RAMSES_fortran_read(fstream &F, RAMSESFLOAT *f){
    int dummy,byteoffset=0;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    F.read((char*)f,dummy); byteoffset+=dummy;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    return byteoffset;
}
int RAMSES_fortran_read(fstream &F, char *c){  // This is needed to read fortran int(1)
    int dummy,byteoffset=0;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    F.read((char*)c,dummy); byteoffset+=dummy;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    return byteoffset;
}

int RAMSES_fortran_skip(fstream &F, int nskips){
    int dummy,byteoffset=0;
    for (int i=0;i<nskips;i++) {
        F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
        F.seekg(dummy,ios::cur); byteoffset+=dummy;
        F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    }
    return byteoffset;
}

vector<string> RAMSES_read_descriptor(char *bufname){
    fstream Fdesc;
    string fieldname, fieldtype, tmpnum, stringbuf;
    int fieldnum;
    vector<string> vfields;
    stringstream fieldstream;
    Fdesc.open(bufname, ios::in);
    getline(Fdesc,stringbuf); // version comment
    getline(Fdesc,stringbuf); // header comment
    while (getline(Fdesc, stringbuf))
    {
	fieldstream.str(stringbuf);
	// Read the field ID, name and type
	getline(fieldstream, tmpnum, ',');
	getline(fieldstream, fieldname, ',');
	getline(fieldstream, fieldtype, ',');
	fieldstream.clear();
	// Clean a bit the strings, because C++
	fieldname.erase( remove_if( fieldname.begin(), fieldname.end(), ::isspace ), fieldname.end() );
	fieldtype.erase( remove_if( fieldtype.begin(), fieldtype.end(), ::isspace ), fieldtype.end() );
	fieldnum = stoi(tmpnum);
	// We need to store those in an array (fieldname and fieldtype are probably enough)
	vfields.push_back(fieldname);
    }
    Fdesc.close();
    return vfields;
}

Int_t RAMSES_get_nbodies(char *fname, int ptype, Options &opt)
{
    char buf[2000],buf1[2000],buf2[2000];
    double * dummy_age, * dummy_mass;
    char * dummy_type;
    double OmegaM, OmegaB;
    int totalghost = 0;
    int totalstars = 0;
    int totaldm    = 0;
    int alltotal   = 0;
    int ghoststars = 0;
    string stringbuf;
    int ninputoffset = 0;
    sprintf(buf1,"%s/amr_%s.out00001",fname,opt.ramsessnapname);
    sprintf(buf2,"%s/amr_%s.out",fname,opt.ramsessnapname);
    if (FileExists(buf1)) sprintf(buf,"%s",buf1);
    else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
    else {
        printf("Error. Can't find AMR data \nneither as `%s'\nnor as `%s'\n\n", buf1, buf2);
        exit(9);
    }
    //if gas searched in some fashion then load amr/hydro data
    if (opt.partsearchtype==PSTGAS||opt.partsearchtype==PSTALL||(opt.partsearchtype==PSTDARK&&opt.iBaryonSearch)) {
    sprintf(buf1,"%s/hydro_%s.out00001",fname,opt.ramsessnapname);
    sprintf(buf2,"%s/hydro_%s.out",fname,opt.ramsessnapname);
    if (FileExists(buf1)) sprintf(buf,"%s",buf1);
    else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
    else {
        printf("Error. Can't find Hydro data \nneither as `%s'\nnor as `%s'\n\n", buf1, buf2);
        exit(9);
    }
    }
    sprintf(buf1,"%s/part_%s.out00001",fname,opt.ramsessnapname);
    sprintf(buf2,"%s/part_%s.out",fname,opt.ramsessnapname);
    if (FileExists(buf1)) sprintf(buf,"%s",buf1);
    else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
    else {
        printf("Error. Can't find Particle data \nneither as `%s'\nnor as `%s'\n\n", buf1, buf2);
        exit(9);
    }

    fstream Framses;
    RAMSES_Header ramses_header_info;
    //buffers to load data
    int intbuff[NRAMSESTYPE];
    long long longbuff[NRAMSESTYPE];
    int i,j,k,ireaderror=0;
    Int_t nbodies=0;
    //IntType inttype;
    int dummy;

    int nusetypes,usetypes[NRAMSESTYPE];

    if (ptype==PSTALL) {nusetypes=4;usetypes[0]=RAMSESGASTYPE;usetypes[1]=RAMSESDMTYPE;usetypes[2]=RAMSESSTARTYPE;usetypes[3]=RAMSESBHTYPE;}
    else if (ptype==PSTDARK) {nusetypes=1;usetypes[0]=RAMSESDMTYPE;}
    else if (ptype==PSTGAS) {nusetypes=1;usetypes[0]=RAMSESGASTYPE;}
    else if (ptype==PSTSTAR) {nusetypes=1;usetypes[0]=RAMSESSTARTYPE;}
    else if (ptype==PSTBH) {nusetypes=1;usetypes[0]=RAMSESBHTYPE;}


    for (j = 0; j < NRAMSESTYPE; j++) ramses_header_info.npartTotal[j] = 0;

    //Open the specified file and the specified dataset in the file.
    //first open amr data
    sprintf(buf1,"%s/amr_%s.out00001",fname,opt.ramsessnapname);
    sprintf(buf2,"%s/amr_%s.out",fname,opt.ramsessnapname);
    if (FileExists(buf1)) sprintf(buf,"%s",buf1);
    else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
    Framses.open(buf, ios::binary|ios::in);
    //read header info
    //in ramses, the AMR files have the following header
    /*
    integer: ncpu                      -> number of files in total
    integer: ndim                      -> number of dimensions (3)
    integer: nx,ny,nz                  -> number of coarse cells in each dimension
    integer: nlevelmax                 -> highest level allowed
    integer: ngridmax                  -> size of the oct arrays
    integer: nboundary                 -> number of boundaries (> 0 if non-periodic)
    integer: ngrid_current             -> number of octs used
    float: boxlen                      -> box length in code units
    integer: noutput, ioutput, ifout   -> nb of outputs, increments for outputs
    float[noutput]: tout               -> output times
    float[noutput]: aout               -> output expansion factors
    float: t                           -> current time
    float[nlevelmax]: dtold            -> old timestep at each level
    float[nlevelmax]: dtnew            -> new timestep at each level
    integer: nstep, nstep_coarse       -> timestep (fine and coarse)
    float: einit, mass_tot_0, rho_tot  -> initial total energy and mass?, mean density
    float: omega_m, omega_l, omega_k, omega_b, h0, aexp_ini, boxlen_ini -> cosmology
    float: aexp, hexp, aexp_old, epot_tot_int, epot_tot_old
    float: mass_sph                    -> reference mass for baryons (mdm*Ob/Om in cosmo)
    */

    // Number of files
    RAMSES_fortran_read(Framses, ramses_header_info.num_files);
    opt.num_files=ramses_header_info.num_files;

    // Number of dimensions
    RAMSES_fortran_read(Framses, ramses_header_info.ndim);

    // Number of coarse cells in each dimension
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.nx, sizeof(int));
    Framses.read((char*)&ramses_header_info.ny, sizeof(int));
    Framses.read((char*)&ramses_header_info.nz, sizeof(int));
    Framses.read((char*)&dummy, sizeof(dummy));

    // Maximum refinement level
    RAMSES_fortran_read(Framses, ramses_header_info.nlevelmax);

    // Max number of grids in each array
    RAMSES_fortran_read(Framses, ramses_header_info.ngridmax);

    // Number of boundaries
    RAMSES_fortran_read(Framses, ramses_header_info.nboundary);

    //this would be number of active grids but ignore for now
    RAMSES_fortran_skip(Framses);

    // Boxsize
    RAMSES_fortran_read(Framses, ramses_header_info.BoxSize);

    //now skip 10 blocks (until mass_sph)
    RAMSES_fortran_skip(Framses, 10);
    
    // mass_sph (not really the gas mass, though...)
    RAMSES_fortran_read(Framses, ramses_header_info.mass[RAMSESGASTYPE]);

    Framses.close();

    //reopen to get number of amr cells might need to alter to read grid information and what cells have no so-called son cells
    if (opt.partsearchtype==PSTGAS||opt.partsearchtype==PSTALL||(opt.partsearchtype==PSTDARK&&opt.iBaryonSearch)) {
	for (i=0;i<ramses_header_info.num_files;i++) {
	    sprintf(buf1,"%s/amr_%s.out%05d",fname,opt.ramsessnapname,i+1);
	    sprintf(buf2,"%s/amr_%s.out",fname,opt.ramsessnapname);
	    if (FileExists(buf1)) sprintf(buf,"%s",buf1);
	    else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
	    Framses.open(buf, ios::binary|ios::in);
	    // Skip the first 6 entries
	    RAMSES_fortran_skip(Framses, 6);
	    // Read the number of oct used
	    RAMSES_fortran_read(Framses, ramses_header_info.npart[RAMSESGASTYPE]);
	    // There are 2**ndim cells per oct
	    ramses_header_info.npart[RAMSESGASTYPE] *= pow(2,ramses_header_info.ndim);
	    ramses_header_info.npartTotal[RAMSESGASTYPE] += ramses_header_info.npart[RAMSESGASTYPE];
	    Framses.close();
	}
	
	//now hydro header data
	sprintf(buf1,"%s/hydro_%s.out00001",fname,opt.ramsessnapname);
	sprintf(buf2,"%s/hydro_%s.out",fname,opt.ramsessnapname);
	if (FileExists(buf1)) sprintf(buf,"%s",buf1);
	else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
	Framses.open(buf, ios::binary|ios::in);

	// Skip ncpu
	RAMSES_fortran_skip(Framses);
	// Read nvar
	RAMSES_fortran_read(Framses, ramses_header_info.nvarh);
	// Skip ndim, nlevelmax, nboundary
	RAMSES_fortran_skip(Framses, 3);
	// Read gamma
	RAMSES_fortran_read(Framses, ramses_header_info.gamma_index);
	Framses.close();
    }

    // List existing fields in part_file_descriptor.txt
    vector<string> partfields;
    sprintf(buf1,"%s/part_file_descriptor.txt", fname);
    partfields = RAMSES_read_descriptor(buf1);

    // Read particle info
    for (i=0;i<ramses_header_info.num_files;i++)
    {
        sprintf(buf1,"%s/part_%s.out%05d",fname,opt.ramsessnapname,i+1);
        sprintf(buf2,"%s/part_%s.out",fname,opt.ramsessnapname);
        if (FileExists(buf1)) sprintf(buf,"%s",buf1);
        else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
        Framses.open(buf, ios::binary|ios::in); // Framses now reads from PART files

        ramses_header_info.npart[RAMSESDMTYPE]   = 0;
        ramses_header_info.npart[RAMSESSTARTYPE] = 0;
        ramses_header_info.npart[RAMSESSINKTYPE] = 0;

	// Skip ncpu and ndim
	RAMSES_fortran_skip(Framses, 2);
        // Total number of LOCAL particles (npart)
	RAMSES_fortran_read(Framses, ramses_header_info.npartlocal);
        // Skip Random seeds
	RAMSES_fortran_skip(Framses, 1);
        // Total number of Stars over all processors
	RAMSES_fortran_read(Framses, ramses_header_info.nstarTotal);
        // Skip total mass of stars and total lost mass of stars
	RAMSES_fortran_skip(Framses, 2);
        // Number of sink particles over the whole simulation (all are included in
        // all processors)
	RAMSES_fortran_read(Framses, ramses_header_info.npartTotal[RAMSESSINKTYPE]);

	//allocate memory to store types
	dummy_type = new char [ramses_header_info.npartlocal];
	
	for(j=0; j<partfields.size(); ++j)
	{
	    if (partfields[j] == "family") RAMSES_fortran_read(Framses, dummy_type);
	    else RAMSES_fortran_skip(Framses, 1);
	}

        ghoststars = 0;
        for (j = 0; j < ramses_header_info.npartlocal; j++)
        {
	    if (dummy_type[j] == 1)
	    {
		// This is a DM particle
                ramses_header_info.npart[RAMSESDMTYPE]++;
	    }
            else if (dummy_type[j] == 2)
	    {
		// Any form of star particles
		ramses_header_info.npart[RAMSESSTARTYPE]++;
	    }
	    else
	    {
		// Can be a cloud, or whatever
		// TODO: account for them?
		ghoststars++;
	    }
	}
        delete [] dummy_type;
        Framses.close();

        totalghost += ghoststars;
        totalstars += ramses_header_info.npart[RAMSESSTARTYPE];
        totaldm    += ramses_header_info.npart[RAMSESDMTYPE];
        alltotal   += ramses_header_info.npartlocal;

        //now with information loaded, set totals
        ramses_header_info.npartTotal[RAMSESDMTYPE]+=ramses_header_info.npart[RAMSESDMTYPE];
        ramses_header_info.npartTotal[RAMSESSTARTYPE]+=ramses_header_info.npart[RAMSESSTARTYPE];
        ramses_header_info.npartTotal[RAMSESSINKTYPE]+=ramses_header_info.npart[RAMSESSINKTYPE];
    }
    for(j=0, nbodies=0; j<nusetypes; j++) {
        k=usetypes[j];
        nbodies+=ramses_header_info.npartTotal[k];
        //nbodies+=((long long)(ramses_header_info.npartTotalHW[k]) << 32);
    }

    for (j=0;j<NPARTTYPES;j++) opt.numpart[j]=0;
    if (ptype==PSTALL || ptype==PSTDARK) opt.numpart[DARKTYPE]=ramses_header_info.npartTotal[RAMSESDMTYPE];
    if (ptype==PSTALL || ptype==PSTGAS) opt.numpart[GASTYPE]=ramses_header_info.npartTotal[RAMSESGASTYPE];
    if (ptype==PSTALL || ptype==PSTSTAR) opt.numpart[STARTYPE]=ramses_header_info.npartTotal[RAMSESSTARTYPE];
    if (ptype==PSTALL || ptype==PSTBH) opt.numpart[BHTYPE]=ramses_header_info.npartTotal[RAMSESSINKTYPE];

    return nbodies;

}

/// Reads a ramses file. If cosmological simulation uses cosmology (generally
/// assuming LCDM or small deviations from this) to estimate the mean interparticle
/// spacing and scales physical linking length passed by this distance. Also reads
/// header and overrides passed cosmological parameters with ones stored in header.
void ReadRamses(Options &opt, vector<Particle> &Part, const Int_t nbodies, Particle *&Pbaryons, Int_t nbaryons)
{
    char buf[2000],buf1[2000],buf2[2000];
    string stringbuf,orderingstring;
    fstream Finfo;
    fstream *Famr;
    fstream *Fhydro;
    fstream *Fpart;
#ifdef BHON
    fstream Fsink;  // There will be only one of these
#endif
    RAMSES_Header *header;
    int intbuff[NRAMSESTYPE];
    long long longbuff[NRAMSESTYPE];
    int i,j,k,n,idim,ivar,igrid,ireaderror=0;
    Int_t count2,bcount2,gascount;
    //IntType inttype;
    int dummy,byteoffset;
    Double_t MP_DM=MAXVALUE,LN=1.0,N_DM,MP_B=0;
    double z,aadjust,Hubble,Hubbleflow;
    Double_t mscale,lscale,lvscale,rhoscale;
    Double_t mtemp,utemp,rhotemp,Ztemp,Ttemp;
#ifdef HIGHRES
    Double_t reftemp = -1.0;
#endif
    Coordinate xpos,vpos;
    RAMSESFLOAT xtemp[3],vtemp[3];
    RAMSESIDTYPE idval, offsetID;
    int typeval;
    RAMSESFLOAT ageval,metval;
    int *ngridlevel,*ngridbound,*ngridfile;
    int *ncache, *numbl, *numbb;
    int lmin=1000000,lmax=0;

    int ninputoffset = 0;
    int ifirstfile=0,*ireadfile,ibuf=0;
    Int_t ibufindex;
    int *ireadtask,*readtaskID;
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
    ireadfile=new int[opt.num_files];
    for (i=0;i<opt.num_files;i++) ireadfile[i]=1;
#else
    MPI_Bcast (&(opt.num_files), sizeof(opt.num_files), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Barrier (MPI_COMM_WORLD);
#endif
    ///\todo because of the stupid fortran format, easier if chunksize is BIG so that
    ///number of particles local to a file are smaller
    Int_t chunksize=RAMSESCHUNKSIZE,nchunk;
    RAMSESFLOAT *xtempchunk, *vtempchunk, *mtempchunk, *sphtempchunk, *agetempchunk, *mettempchunk, *hydrotempchunk;
    RAMSESIDTYPE *idvalchunk, *levelchunk;
    char *typechunk;
    int *icellchunk;

    Famr       = new fstream[opt.num_files];
    Fhydro     = new fstream[opt.num_files];
    Fpart      = new fstream[opt.num_files];
    header     = new RAMSES_Header[opt.num_files];

    Particle *Pbuf;
#ifdef USEMPI
    MPI_Status status;
    MPI_Comm mpi_comm_read;

    vector<Particle> *Preadbuf;
    //for parallel io
    Int_t BufSize=opt.mpiparticlebufsize;
    Int_t Nlocalbuf,*Nbuf, *Nreadbuf,*nreadoffset;
    Int_t *Nlocalthreadbuf,Nlocaltotalbuf;
    int *irecv, sendTask,recvTask,irecvflag, *mpi_irecvflag;
    MPI_Request *mpi_request;
    Int_t inreadsend,totreadsend;
    Int_t *mpi_nsend_baryon;
    Int_t *mpi_nsend_readthread;
    Int_t *mpi_nsend_readthread_baryon;
    if (opt.iBaryonSearch) mpi_nsend_baryon=new Int_t[NProcs*NProcs];
    if (opt.nsnapread>1) {
        mpi_nsend_readthread=new Int_t[opt.nsnapread*opt.nsnapread];
        if (opt.iBaryonSearch) mpi_nsend_readthread_baryon=new Int_t[opt.nsnapread*opt.nsnapread];
    }

    Nbuf=new Int_t[NProcs];
    nreadoffset=new Int_t[opt.nsnapread];
    ireadtask=new int[NProcs];
    readtaskID=new int[opt.nsnapread];
    MPIDistributeReadTasks(opt,ireadtask,readtaskID);
    MPI_Comm_split(MPI_COMM_WORLD, (ireadtask[ThisTask]>=0), ThisTask, &mpi_comm_read);
    if (ireadtask[ThisTask]>=0)
    {
        //to temporarily store data
        Pbuf=new Particle[BufSize*NProcs];
        Nreadbuf=new Int_t[opt.nsnapread];
        for (int j=0;j<NProcs;j++) Nbuf[j]=0;
        for (int j=0;j<opt.nsnapread;j++) Nreadbuf[j]=0;
        if (opt.nsnapread>1){
            Preadbuf=new vector<Particle>[opt.nsnapread];
            for (int j=0;j<opt.nsnapread;j++) Preadbuf[j].reserve(BufSize);
        }
        //to determine which files the thread should read
        ireadfile=new int[opt.num_files];
        ifirstfile=MPISetFilesRead(opt,ireadfile,ireadtask);
        inreadsend=0;
        for (int j=0;j<opt.num_files;j++) inreadsend+=ireadfile[j];
        MPI_Allreduce(&inreadsend,&totreadsend,1,MPI_Int_t,MPI_MIN,mpi_comm_read);
    }
    else {
        Nlocalthreadbuf=new Int_t[opt.nsnapread];
        irecv=new int[opt.nsnapread];
        mpi_irecvflag=new int[opt.nsnapread];
        for (i=0;i<opt.nsnapread;i++) irecv[i]=1;
        mpi_request=new MPI_Request[opt.nsnapread];
    }
    Nlocal=0;
    if (opt.iBaryonSearch) Nlocalbaryon[0]=0;

    if (ireadtask[ThisTask]>=0) {
#endif

    //first read cosmological information
    sprintf(buf1,"%s/info_%s.txt",opt.fname,opt.ramsessnapname);
    Finfo.open(buf1, ios::in);
    getline(Finfo,stringbuf);//nfiles
    getline(Finfo,stringbuf);//ndim
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].levelmin;
    Finfo.ignore();
    getline(Finfo,stringbuf);//lmax
    getline(Finfo,stringbuf);//ngridmax
    getline(Finfo,stringbuf);//nstep
    getline(Finfo,stringbuf);//blank

    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].BoxSize;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].time;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].aexp;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].HubbleParam;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].Omegam;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].OmegaLambda;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].Omegak;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].Omegab;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].scale_l;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].scale_d;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].scale_t;
    Finfo.ignore();

    //convert boxsize to comoving kpc/h
    header[ifirstfile].BoxSize*=header[ifirstfile].scale_l/3.086e21/header[ifirstfile].aexp*header[ifirstfile].HubbleParam/100.0;
    getline(Finfo,stringbuf);
    Finfo>>stringbuf>>orderingstring;
    Finfo.ignore();
    getline(Finfo,stringbuf);
    Finfo.close();

    opt.p            = header[ifirstfile].BoxSize;
    opt.a            = header[ifirstfile].aexp;
    opt.Omega_m      = header[ifirstfile].Omegam;
    opt.Omega_Lambda = header[ifirstfile].OmegaLambda;
    opt.Omega_b      = header[ifirstfile].Omegab;
    opt.h            = header[ifirstfile].HubbleParam/100.0;
    opt.Omega_cdm    = opt.Omega_m-opt.Omega_b;
    //set hubble unit to km/s/kpc
    opt.H = 0.1;
    //set Gravity to value for kpc (km/s)^2 / solar mass
    opt.G = 4.30211349e-6;
    //and for now fix the units
    opt.lengthtokpc=opt.velocitytokms=opt.masstosolarmass=1.0;

    //Hubble flow
    if (opt.comove) aadjust=1.0;
    else aadjust=opt.a;
    CalcOmegak(opt);
    Hubble=GetHubble(opt, aadjust);
    CalcCriticalDensity(opt, aadjust);
    CalcBackgroundDensity(opt, aadjust);
    CalcVirBN98(opt,aadjust);
    //if opt.virlevel<0, then use virial overdensity based on Bryan and Norman 1997 virialization level is given by
    if (opt.virlevel<0) opt.virlevel=opt.virBN98;
    PrintCosmology(opt);

    //adjust length scale so that convert from 0 to 1 (box units) to kpc comoving
    //to scale mpi domains correctly need to store in opt.lengthinputconversion the box size in comoving little h value
    //opt.lengthinputconversion= opt.p*opt.h/opt.a;
    opt.lengthinputconversion = header[ifirstfile].BoxSize;
    //adjust velocity scale to that ramses is converted to km/s from which you can convert again;
    opt.velocityinputconversion = header[ifirstfile].scale_l/header[ifirstfile].scale_t*1e-5;

    //convert mass from code units to Solar Masses
    mscale   = header[ifirstfile].scale_d * (header[ifirstfile].scale_l * header[ifirstfile].scale_l * header[ifirstfile].scale_l) / 1.988e+33;

    //convert length from code units to kpc (this lscale includes the physical box size)
    lscale   = header[ifirstfile].scale_l/3.086e21;

    //convert velocity from code units to km/s
    lvscale  = header[ifirstfile].scale_l/header[ifirstfile].scale_t*1e-5;

    //convert density to Msun/kpc^3
    rhoscale = mscale/(lscale*lscale*lscale);

    //ignore hubble flow
    Hubbleflow=0.;

    //for (int j=0;j<NPARTTYPES;j++) nbodies+=opt.numpart[j];
    cout<<"Particle system contains "<<nbodies<<" particles (of interest) at is at time "<<opt.a<<" in a box of size "<<opt.p<<endl;

    //number of DM particles
    //NOTE: this assumes a uniform box resolution. However this is not used in the rest of this function
    N_DM = opt.Neff*opt.Neff*opt.Neff;

    // Define offset for the particle IDs
    // for 1000 particles in total:
    // - gas would be indexed 1, 2, 3
    // - DM particles would be indexed 1001, 1002, 1003, etc
    // - star would be indexed 2001, 2002, etc
    // - BH (?) would be indexed 3001, 3002, 3003, etc
    offsetID = powf(10.0f, ceilf(log10f(nbodies)+1));

    //grab from the first particle file the dimensions of the arrays and also the number of cpus (should be number of files)
    sprintf(buf1,"%s/part_%s.out00001",opt.fname,opt.ramsessnapname);
    Fpart[ifirstfile].open(buf1, ios::binary|ios::in);
    RAMSES_fortran_read(Fpart[ifirstfile],header[ifirstfile].nfiles);
    RAMSES_fortran_read(Fpart[ifirstfile],header[ifirstfile].ndim);
    //adjust the number of files
    opt.num_files=header[ifirstfile].nfiles;
    Fpart[ifirstfile].close();
#ifdef USEMPI
    //now read tasks prepped and can read files to send information
    }
#endif

    //if not only gas being searched open particle data
    count2=bcount2=0;
    if (opt.partsearchtype!=PSTGAS && opt.partsearchtype!=PSTBH) {
#ifdef USEMPI
	if (ireadtask[ThisTask]>=0) {
	    inreadsend=0;
#endif

	// List existing particle fields from part_file_descriptor.txt
	vector<string> partfields;
	sprintf(buf1,"%s/part_file_descriptor.txt", opt.fname);
	partfields = RAMSES_read_descriptor(buf1);

	//read particle files consists of positions,velocities, mass, id, and level (along with ages and met if some flags set)
	for (i=0;i<opt.num_files;i++) {
	    if (ireadfile[i]) {
		sprintf(buf1,"%s/part_%s.out%05d",opt.fname,opt.ramsessnapname,i+1);
		sprintf(buf2,"%s/part_%s.out",opt.fname,opt.ramsessnapname);
		if (FileExists(buf1)) sprintf(buf,"%s",buf1);
		else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
		Fpart[i].open(buf, ios::binary|ios::in);
		
		//skip header information in each file save for number in the file
		// //@{
		byteoffset=0;
		// ncpus
		byteoffset+=RAMSES_fortran_skip(Fpart[i]);
		// ndims
		byteoffset+=RAMSES_fortran_skip(Fpart[i]);
		// store number of particles locally in file
		byteoffset+=RAMSES_fortran_read(Fpart[i],header[i].npartlocal);
		// skip local seeds, nstartot, mstartot, mstarlost, nsink
		byteoffset+=RAMSES_fortran_skip(Fpart[i],5);
		//byteoffset now stores size of header offset for particles

		//data loaded into memory in chunks
		chunksize    = nchunk = header[i].npartlocal;
		ninputoffset = 0;
		xtempchunk   = new RAMSESFLOAT  [3*chunksize];
		vtempchunk   = new RAMSESFLOAT  [3*chunksize];
		mtempchunk   = new RAMSESFLOAT  [chunksize];
		idvalchunk   = new RAMSESIDTYPE [chunksize];
		levelchunk   = new RAMSESIDTYPE [chunksize];
		typechunk    = new char         [chunksize];
		agetempchunk = new RAMSESFLOAT  [chunksize];
		mettempchunk = new RAMSESFLOAT  [chunksize];


		for (j=0;j<partfields.size();++j)
		{
		    // The order should not be important, since we skip the fields we do not read
		    if      (partfields[j] == string("position_x"))  {RAMSES_fortran_read(Fpart[i],&xtempchunk[0*nchunk]);}
		    else if (partfields[j] == string("position_y"))  {RAMSES_fortran_read(Fpart[i],&xtempchunk[1*nchunk]);}
		    else if (partfields[j] == string("position_z"))  {RAMSES_fortran_read(Fpart[i],&xtempchunk[2*nchunk]);}
		    else if (partfields[j] == string("velocity_x"))  {RAMSES_fortran_read(Fpart[i],&vtempchunk[0*nchunk]);}
		    else if (partfields[j] == string("velocity_y"))  {RAMSES_fortran_read(Fpart[i],&vtempchunk[1*nchunk]);}
		    else if (partfields[j] == string("velocity_z"))  {RAMSES_fortran_read(Fpart[i],&vtempchunk[2*nchunk]);}
		    else if (partfields[j] == string("mass"))        {RAMSES_fortran_read(Fpart[i],mtempchunk)           ;}
		    else if (partfields[j] == string("identity"))    {RAMSES_fortran_read(Fpart[i],idvalchunk)           ;}
		    else if (partfields[j] == string("levelp"))      {RAMSES_fortran_read(Fpart[i],levelchunk)           ;}
		    else if (partfields[j] == string("family"))      {RAMSES_fortran_read(Fpart[i],typechunk)            ;}
		    else if (partfields[j] == string("birth_time"))  {RAMSES_fortran_read(Fpart[i],agetempchunk)         ;}
		    else if (partfields[j] == string("metallicity")) {RAMSES_fortran_read(Fpart[i],mettempchunk)         ;}
		    else {RAMSES_fortran_skip(Fpart[i]);}
		}

		// Finally create particle
		for (int nn=0;nn<nchunk;nn++)
		{
		    if ((typechunk[nn] != 1) && (typechunk[nn] != 2))
		    {
			//  GHOST PARTICLE!!!
			// Basically, neither star nor DM: can be cloud, dust, whatever
			// used to be: ((typechunk[nn] != 1) && (agetempchunk[nn] == 0.0))
			// MT: this could be a cloud particle (which should? have type == 3)
			// MT: FIXME BH
		    }
		    else
		    {
			xtemp[0] = xtempchunk[nn];
			xtemp[1] = xtempchunk[nn+nchunk];
			xtemp[2] = xtempchunk[nn+2*nchunk];

			vtemp[0] = vtempchunk[nn];
			vtemp[1] = vtempchunk[nn+nchunk];
			vtemp[2] = vtempchunk[nn+2*nchunk];

			idval = idvalchunk[nn];

			///Need to check this for correct 'endianness'
//             for (int kk=0;kk<3;kk++) {xtemp[kk]=LittleRAMSESFLOAT(xtemp[kk]);vtemp[kk]=LittleRAMSESFLOAT(vtemp[kk]);}
#ifndef NOMASS
			mtemp=mtempchunk[nn];
#else
			mtemp=1.0;
#endif
			ageval = agetempchunk[nn];
			metval = mettempchunk[nn];
			if (typechunk[nn] == 1) typeval = DARKTYPE;
			else if (typechunk[nn] == 2) typeval = STARTYPE;
			else typeval=-1; // FIXME: is -1 ok as a undefined value?
#ifdef HIGHRES
                        // For high-res simulation, find mass of highres DM particles
                        if (mtemp>0 && mtemp < MP_DM && typeval==DARKTYPE) MP_DM=mtemp;
#endif

		
#ifdef USEMPI
			//determine processor this particle belongs on based on its spatial position
			ibuf=MPIGetParticlesProcessor(xtemp[0],xtemp[1],xtemp[2]);
			ibufindex=ibuf*BufSize+Nbuf[ibuf];
#endif
			//reset hydro quantities of buffer
#ifdef USEMPI
#ifdef GASON
			Pbuf[ibufindex].SetU(0);
#ifdef STARON
			Pbuf[ibufindex].SetSFR(0);
			Pbuf[ibufindex].SetZmet(metval);  // should reset to 0?
#endif
#endif
#ifdef STARON
			Pbuf[ibufindex].SetZmet(metval);  // should reset to 0?
			Pbuf[ibufindex].SetTage(ageval); // careful, this is in weird units
#endif
#ifdef BHON
#endif
#endif
		
			if (opt.partsearchtype==PSTALL) {
#ifdef USEMPI
			    Pbuf[ibufindex]=Particle(mtemp*mscale,xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
						     vtemp[0]*opt.V+Hubbleflow*xtemp[0],
						     vtemp[1]*opt.V+Hubbleflow*xtemp[1],
						     vtemp[2]*opt.V+Hubbleflow*xtemp[2],
						     count2,typeval);
			    Pbuf[ibufindex].SetPID(idval + typeval*offsetID);
#ifdef EXTRAINPUTINFO
			    if (opt.iextendedoutput)
			    {
				Part[ibufindex].SetInputFileID(i);
				Part[ibufindex].SetInputIndexInFile(nn+ninputoffset);
			    }
#endif
			    Nbuf[ibuf]++;
			    MPIAddParticletoAppropriateBuffer(ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, Nlocal, Part.data(), Nreadbuf, Preadbuf);
#else
			    Part[count2]=Particle(mtemp*mscale,
						  xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
						  vtemp[0]*opt.V+Hubbleflow*xtemp[0],
						  vtemp[1]*opt.V+Hubbleflow*xtemp[1],
						  vtemp[2]*opt.V+Hubbleflow*xtemp[2],
						  count2,typeval);
			    Part[count2].SetPID(idval + typeval*offsetID);
#ifdef EXTRAINPUTINFO
			    if (opt.iextendedoutput)
			    {
				Part[count2].SetInputFileID(i);
				Part[count2].SetInputIndexInFile(nn+ninputoffset);
			    }
#endif
#endif
			    count2++;
			}
			else if (opt.partsearchtype==PSTDARK) {
			    if (!(typeval==STARTYPE||typeval==BHTYPE)){
				// Dealing with a DM particle
#ifdef USEMPI
			        Pbuf[ibufindex]=Particle(mtemp*mscale,
				                         xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
				                         vtemp[0]*opt.V+Hubbleflow*xtemp[0],
				                         vtemp[1]*opt.V+Hubbleflow*xtemp[1],
				                         vtemp[2]*opt.V+Hubbleflow*xtemp[2],
				                         count2,DARKTYPE);
				Pbuf[ibufindex].SetPID(idval + typeval*offsetID);
#ifdef EXTRAINPUTINFO
				if (opt.iextendedoutput)
				{
				    Pbuf[ibufindex].SetInputFileID(i);
				    Pbuf[ibufindex].SetInputIndexInFile(nn+ninputoffset);
				}
#endif
				//ensure that store number of particles to be sent to other reading threads
				Nbuf[ibuf]++;
				MPIAddParticletoAppropriateBuffer(ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, Nlocal, Part.data(), Nreadbuf, Preadbuf);
#else
				Part[count2]=Particle(mtemp*mscale,
						      xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
						      vtemp[0]*opt.V+Hubbleflow*xtemp[0],
						      vtemp[1]*opt.V+Hubbleflow*xtemp[1],
						      vtemp[2]*opt.V+Hubbleflow*xtemp[2],
						      count2,typeval);
				Part[count2].SetPID(idval + typeval*offsetID);
#ifdef EXTRAINPUTINFO
				if (opt.iextendedoutput)
				{
				    Part[count2].SetInputFileID(i);
				    Part[count2].SetInputIndexInFile(nn+ninputoffset);
				}
#endif
#endif
				count2++;
			    }
			    else if (opt.iBaryonSearch) {
				// Baryon particle AND we search for them
#ifdef USEMPI
				Pbuf[ibufindex]=Particle(mtemp*mscale,
							 xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
							 vtemp[0]*opt.V+Hubbleflow*xtemp[0],
							 vtemp[1]*opt.V+Hubbleflow*xtemp[1],
							 vtemp[2]*opt.V+Hubbleflow*xtemp[2],
							 bcount2);  //count2 -> bcount2?
				Pbuf[ibufindex].SetPID(idval + typeval*offsetID);
#ifdef EXTRAINPUTINFO
				if (opt.iextendedoutput)
				{
				    Pbuf[ibufindex].SetInputFileID(i);
				    Pbuf[ibufindex].SetInputIndexInFile(nn+ninputoffset);
				}
#endif
				if (typeval==STARTYPE) Pbuf[ibufindex].SetType(STARTYPE);
				else if (typeval==BHTYPE) Pbuf[ibufindex].SetType(BHTYPE);
				//ensure that store number of particles to be sent to the reading threads
				Nbuf[ibuf]++;
				if (ibuf==ThisTask) {
				    if (typeval==STARTYPE) Nlocalbaryon[2]++;
				    else if (typeval==BHTYPE) Nlocalbaryon[3]++;
				}
				MPIAddParticletoAppropriateBuffer(ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, Nlocalbaryon[0], Pbaryons, Nreadbuf, Preadbuf);
#else
				Pbaryons[bcount2]=Particle(mtemp*mscale,
							   xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
							   vtemp[0]*opt.V+Hubbleflow*xtemp[0],
							   vtemp[1]*opt.V+Hubbleflow*xtemp[1],
							   vtemp[2]*opt.V+Hubbleflow*xtemp[2],
							   bcount2,typeval); //count2 -> bcount2?
				Pbaryons[bcount2].SetPID(idval + typeval*offsetID);
#ifdef EXTRAINPUTINFO
				if (opt.iextendedoutput)
				{
				    Pbaryons[bcount2].SetInputFileID(i);
				    Pbaryons[bcount2].SetInputIndexInFile(nn+ninputoffset);
				}
#endif
#endif
				bcount2++;
			    }
			}
			else if (opt.partsearchtype==PSTSTAR) {
			    if (typeval==STARTYPE) {
#ifdef USEMPI
				//if using MPI, determine proccessor and place in ibuf, store particle in particle buffer and if buffer full, broadcast data
				//unless ibuf is 0, then just store locally
				Pbuf[ibufindex]=Particle(mtemp*mscale,
							 xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
							 vtemp[0]*opt.V+Hubbleflow*xtemp[0],
							 vtemp[1]*opt.V+Hubbleflow*xtemp[1],
							 vtemp[2]*opt.V+Hubbleflow*xtemp[2],
							 count2,STARTYPE);
				//ensure that store number of particles to be sent to the reading threads
				Pbuf[ibufindex].SetPID(idval + typeval*offsetID);
#ifdef EXTRAINPUTINFO
				if (opt.iextendedoutput)
				{
				    Pbuf[ibufindex].SetInputFileID(i);
				    Pbuf[ibufindex].SetInputIndexInFile(nn+ninputoffset);
				}
#endif
				Nbuf[ibuf]++;
				MPIAddParticletoAppropriateBuffer(ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, Nlocal, Part.data(), Nreadbuf, Preadbuf);
#else
				Part[count2]=Particle(mtemp*mscale,
						      xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
						      vtemp[0]*opt.V+Hubbleflow*xtemp[0],
						      vtemp[1]*opt.V+Hubbleflow*xtemp[1],
						      vtemp[2]*opt.V+Hubbleflow*xtemp[2],
						      count2,typeval);
				Part[count2].SetPID(idval + typeval*offsetID);
#ifdef EXTRAINPUTINFO
				if (opt.iextendedoutput)
				{
				    Part[count2].SetInputFileID(i);
				    Part[count2].SetInputIndexInFile(nn+ninputoffset);
				}
#endif
#endif
				count2++;
			    }
			}
		    }//end of ghost particle check
		}//end of loop over chunk
	
		delete[] xtempchunk;
		delete[] vtempchunk;
		delete[] mtempchunk;
		delete[] idvalchunk;
		delete[] typechunk;
		delete[] agetempchunk;
		delete[] levelchunk;
		delete[] mettempchunk;
		Fpart[i].close();
#ifdef USEMPI

		//send information between read threads
		if (opt.nsnapread>1&&inreadsend<totreadsend){
		    MPI_Allgather(Nreadbuf, opt.nsnapread, MPI_Int_t, mpi_nsend_readthread, opt.nsnapread, MPI_Int_t, mpi_comm_read);
		    // There should be something similar here for mpi_nsend_readthread_baryon, which is NOT initialized so far: next has an issue because of that
		    MPISendParticlesBetweenReadThreads(opt, Preadbuf, Part.data(), ireadtask, readtaskID, Pbaryons, mpi_comm_read, mpi_nsend_readthread, mpi_nsend_readthread_baryon);
		    inreadsend++;
		    for(ibuf = 0; ibuf < opt.nsnapread; ibuf++) Nreadbuf[ibuf]=0;
		}
#endif
	    }//end of whether reading a file
	}//end of loop over file
#ifdef USEMPI
    //once finished reading the file if there are any particles left in the buffer broadcast them
    for(ibuf = 0; ibuf < NProcs; ibuf++) if (ireadtask[ibuf]<0)
    {
        MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
        if (Nbuf[ibuf]>0) {
            MPI_Ssend(&Pbuf[ibuf*BufSize], sizeof(Particle)*Nbuf[ibuf], MPI_BYTE, ibuf, ibuf, MPI_COMM_WORLD);
            Nbuf[ibuf]=0;
            //last broadcast with Nbuf[ibuf]=0 so that receiver knows no more particles are to be broadcast
            MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t,ibuf,ibuf+NProcs,MPI_COMM_WORLD);
        }
    }
    if (opt.nsnapread>1){
        MPI_Allgather(Nreadbuf, opt.nsnapread, MPI_Int_t, mpi_nsend_readthread, opt.nsnapread, MPI_Int_t, mpi_comm_read);
        MPISendParticlesBetweenReadThreads(opt, Preadbuf, Part.data(), ireadtask, readtaskID, Pbaryons, mpi_comm_read, mpi_nsend_readthread, mpi_nsend_readthread_baryon);
    }
    }//end of ireadtask[ibuf]>0
#endif
#ifdef USEMPI
    else {
        MPIReceiveParticlesFromReadThreads(opt,Pbuf,Part.data(),readtaskID, irecv, mpi_irecvflag, Nlocalthreadbuf, mpi_request,Pbaryons);
    }
#endif
    }

#ifdef BHON
    if (opt.partsearchtype==PSTALL || opt.partsearchtype==PSTBH || (opt.partsearchtype==PSTDARK && opt.iBaryonSearch)) {
	//Read the first sink file, as all of the are identical
	sprintf(buf1,"%s/sink_%s.out00001",opt.fname,opt.ramsessnapname);
	Fsink.open(buf1, ios::binary|ios::in);

	// List existing sink fields from sink_file_descriptor.txt
	vector<string> sinkfields;
	sprintf(buf1,"%s/sink_file_descriptor.txt", opt.fname);
	sinkfields = RAMSES_read_descriptor(buf1);

	// Skip header: nsink and nindsink
	RAMSES_fortran_skip(Fsink, 2);

	/* Read relevant variables:
	   - "identity"   -> sink ID
	   - "mass"       -> sink mass
	   - "position_x" -> position (same with y and z)
	   - "velocity_x" -> velocity (same with y and z)
	   - "birth_time" -> time of formation
	*/
	nchunk = opt.numpart[BHTYPE];
	xtempchunk   = new RAMSESFLOAT  [3*nchunk];
	vtempchunk   = new RAMSESFLOAT  [3*nchunk];
	mtempchunk   = new RAMSESFLOAT  [nchunk];
	idvalchunk   = new RAMSESIDTYPE [nchunk];
	agetempchunk = new RAMSESFLOAT  [nchunk];

	if (opt.iverbose) cout<<ThisTask<<" Reading "<<nchunk<<" sinks  ... "<<endl;
	for (j=0;j<sinkfields.size();++j)
	{
	    // The order should not be important, since we skip the fields we do not read
	    if      (sinkfields[j] == string("identity"))    {RAMSES_fortran_read(Fsink,idvalchunk)           ;}
	    else if (sinkfields[j] == string("mass"))        {RAMSES_fortran_read(Fsink,mtempchunk)           ;}
	    else if (sinkfields[j] == string("position_x"))  {RAMSES_fortran_read(Fsink,&xtempchunk[0*nchunk]);}
	    else if (sinkfields[j] == string("position_y"))  {RAMSES_fortran_read(Fsink,&xtempchunk[1*nchunk]);}
	    else if (sinkfields[j] == string("position_z"))  {RAMSES_fortran_read(Fsink,&xtempchunk[2*nchunk]);}
	    else if (sinkfields[j] == string("velocity_x"))  {RAMSES_fortran_read(Fsink,&vtempchunk[0*nchunk]);}
	    else if (sinkfields[j] == string("velocity_y"))  {RAMSES_fortran_read(Fsink,&vtempchunk[1*nchunk]);}
	    else if (sinkfields[j] == string("velocity_z"))  {RAMSES_fortran_read(Fsink,&vtempchunk[2*nchunk]);}
	    else if (sinkfields[j] == string("birth_time"))  {RAMSES_fortran_read(Fsink,agetempchunk)         ;}
	    else {RAMSES_fortran_skip(Fsink);}
	}
	Fsink.close();

	// Finally create particles
	for (int nn=0;nn<nchunk;nn++) {
	    xtemp[0] = xtempchunk[nn];
	    xtemp[1] = xtempchunk[nn+nchunk];
	    xtemp[2] = xtempchunk[nn+2*nchunk];
	
	    vtemp[0] = vtempchunk[nn];
	    vtemp[1] = vtempchunk[nn+nchunk];
	    vtemp[2] = vtempchunk[nn+2*nchunk];
	
	    idval = idvalchunk[nn];
	
#ifndef NOMASS
	    mtemp=mtempchunk[nn];
#else
	    mtemp=1.0;
#endif
	    ageval = agetempchunk[nn];
		
#ifdef USEMPI
	    // //determine processor this particle belongs on based on its spatial position
	    // ibuf=MPIGetParticlesProcessor(xtemp[0],xtemp[1],xtemp[2]);
	    // ibufindex=ibuf*BufSize+Nbuf[ibuf];
#endif
	    //reset hydro quantities of buffer
#ifdef USEMPI
// #ifdef GASON
// 			Pbuf[ibufindex].SetU(0);
// #ifdef STARON
// 			Pbuf[ibufindex].SetSFR(0);
// 			Pbuf[ibufindex].SetZmet(0);
// #endif
// #endif
// #ifdef STARON
// 			Pbuf[ibufindex].SetZmet(0);
// 			Pbuf[ibufindex].SetTage(ageval); // careful, this is in weird units
// #endif
#endif

	    if (opt.partsearchtype==PSTALL) {
#ifdef USEMPI
// 	    Pbuf[ibufindex]=Particle(mtemp*mscale,xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
// 				     vtemp[0]*opt.V+Hubbleflow*xtemp[0],
// 				     vtemp[1]*opt.V+Hubbleflow*xtemp[1],
// 				     vtemp[2]*opt.V+Hubbleflow*xtemp[2],
// 				     count2,BHTYPE);
// 	    Pbuf[ibufindex].SetPID(idval + BHTYPE*offsetID);
// #ifdef EXTRAINPUTINFO
// 	    if (opt.iextendedoutput)
// 	    {
// 		Part[ibufindex].SetInputFileID(i);
// 		Part[ibufindex].SetInputIndexInFile(nn+ninputoffset);
// 	    }
// #endif
// 	    Nbuf[ibuf]++;
// 	    MPIAddParticletoAppropriateBuffer(ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, Nlocal, Part.data(), Nreadbuf, Preadbuf);
#else
		Part[count2]=Particle(mtemp*mscale,
				      xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
				      vtemp[0]*opt.V+Hubbleflow*xtemp[0],
				      vtemp[1]*opt.V+Hubbleflow*xtemp[1],
				      vtemp[2]*opt.V+Hubbleflow*xtemp[2],
				      count2,BHTYPE);
		Part[count2].SetPID(idval + BHTYPE*offsetID);
#ifdef EXTRAINPUTINFO
		if (opt.iextendedoutput)
		{
		    Part[count2].SetInputFileID(1);
		    Part[count2].SetInputIndexInFile(nn+ninputoffset);
		}
#endif
#endif
		count2++;
	    }
	    else if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
#ifdef USEMPI
// 	    Pbuf[ibufindex]=Particle(mtemp*mscale,
// 				     xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
// 				     vtemp[0]*opt.V+Hubbleflow*xtemp[0],
// 				     vtemp[1]*opt.V+Hubbleflow*xtemp[1],
// 				     vtemp[2]*opt.V+Hubbleflow*xtemp[2],
// 				     bcount2);  //count2 -> bcount2?
// 	    Pbuf[ibufindex].SetPID(idval + BHTYPE*offsetID);
// #ifdef EXTRAINPUTINFO
// 	    if (opt.iextendedoutput)
// 	    {
// 		Pbuf[ibufindex].SetInputFileID(1);
// 		Pbuf[ibufindex].SetInputIndexInFile(nn+ninputoffset);
// 	    }
// #endif
// 	    Pbuf[ibufindex].SetType(BHTYPE);
// 	    //ensure that store number of particles to be sent to the reading threads
// 	    Nbuf[ibuf]++;
// 	    if (ibuf==ThisTask) {
// 		if (k==RAMSESSTARTYPE) Nlocalbaryon[2]++;
// 		else if (k==RAMSESSINKTYPE) Nlocalbaryon[3]++;
// 	    }
// 	    MPIAddParticletoAppropriateBuffer(ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, Nlocalbaryon[0], Pbaryons, Nreadbuf, Preadbuf);
#else
		Pbaryons[bcount2]=Particle(mtemp*mscale,
					   xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
					   vtemp[0]*opt.V+Hubbleflow*xtemp[0],
					   vtemp[1]*opt.V+Hubbleflow*xtemp[1],
					   vtemp[2]*opt.V+Hubbleflow*xtemp[2],
					   bcount2,BHTYPE);
		Pbaryons[bcount2].SetPID(idval + BHTYPE*offsetID);
#ifdef EXTRAINPUTINFO
		if (opt.iextendedoutput)
		{
		    Pbaryons[bcount2].SetInputFileID(1);
		    Pbaryons[bcount2].SetInputIndexInFile(nn+ninputoffset);
		}
#endif
#endif
		bcount2++;
	    }
	    else if (opt.partsearchtype==PSTBH) {
#ifdef USEMPI
// 	    //if using MPI, determine proccessor and place in ibuf, store particle in particle buffer and if buffer full, broadcast data
// 	    //unless ibuf is 0, then just store locally
// 	    Pbuf[ibufindex]=Particle(mtemp*mscale,
// 				     xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
// 				     vtemp[0]*opt.V+Hubbleflow*xtemp[0],
// 				     vtemp[1]*opt.V+Hubbleflow*xtemp[1],
// 				     vtemp[2]*opt.V+Hubbleflow*xtemp[2],
// 				     count2,STARTYPE);
// 	    //ensure that store number of particles to be sent to the reading threads
// 	    Pbuf[ibufindex].SetPID(idval);
// #ifdef EXTRAINPUTINFO
// 	    if (opt.iextendedoutput)
// 	    {
// 		Pbuf[ibufindex].SetInputFileID(1);
// 		Pbuf[ibufindex].SetInputIndexInFile(nn+ninputoffset);
// 	    }
// #endif
// 	    Nbuf[ibuf]++;
// 	    MPIAddParticletoAppropriateBuffer(ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, Nlocal, Part.data(), Nreadbuf, Preadbuf);
#else
		Part[count2]=Particle(mtemp*mscale,
				      xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
				      vtemp[0]*opt.V+Hubbleflow*xtemp[0],
				      vtemp[1]*opt.V+Hubbleflow*xtemp[1],
				      vtemp[2]*opt.V+Hubbleflow*xtemp[2],
				      count2,BHTYPE);
		Part[count2].SetPID(idval);
#ifdef EXTRAINPUTINFO
		if (opt.iextendedoutput)
		{
		    Part[count2].SetInputFileID(1);
		    Part[count2].SetInputIndexInFile(nn+ninputoffset);
		}
#endif
#endif
		count2++;
	    } // end of various case searches
	} // end loop over particles
    
	delete[] xtempchunk;
	delete[] vtempchunk;
	delete[] mtempchunk;
	delete[] idvalchunk;
	delete[] agetempchunk;

	// Need to deal with MPI exchanges as for the regular particles

    } // End of reading BH
#endif

    

    //if gas searched in some fashion then load amr/hydro data
    gascount=0;
    if (opt.partsearchtype==PSTGAS||opt.partsearchtype==PSTALL||(opt.partsearchtype==PSTDARK&&opt.iBaryonSearch)) {
#ifdef USEMPI
    if (ireadtask[ThisTask]>=0) {
    inreadsend=0;
#endif
    // List existing hydro fields from hydro_file_descriptor.txt
    vector<string> hfields;
    sprintf(buf1,"%s/hydro_file_descriptor.txt", opt.fname);
    hfields = RAMSES_read_descriptor(buf1);
    // define the indices for the hydro variables (default = standard ramses)
    int IDrho=0, IDvx=1, IDvy=2, IDvz=3, IDU=4, IDZ=5;
#ifdef HIGHRES
    int IDref=-1;
#endif
    for (ivar=0;ivar<hfields.size();++ivar) {
	if (hfields[ivar] == string("density")) IDrho = ivar;
	if (hfields[ivar] == string("velocity_x")) IDvx = ivar;
	if (hfields[ivar] == string("velocity_y")) IDvy = ivar;
	if (hfields[ivar] == string("velocity_z")) IDvz = ivar;
	if (hfields[ivar] == string("pressure")) IDU = ivar;
	if (hfields[ivar] == string("metallicity")) IDZ = ivar;
#ifdef HIGHRES
        if (hfields[ivar] == string("refinement_scalar")) IDref = ivar;
#endif
    }

    for (i=0;i<opt.num_files;i++) if (ireadfile[i]) {
        sprintf(buf1,"%s/amr_%s.out%05d",opt.fname,opt.ramsessnapname,i+1);
        sprintf(buf2,"%s/amr_%s.out",opt.fname,opt.ramsessnapname);
        if (FileExists(buf1)) sprintf(buf,"%s",buf1);
        else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
        Famr[i].open(buf, ios::binary|ios::in);
        sprintf(buf1,"%s/hydro_%s.out%05d",opt.fname,opt.ramsessnapname,i+1);
        sprintf(buf2,"%s/hydro_%s.out",opt.fname,opt.ramsessnapname);
        if (FileExists(buf1)) sprintf(buf,"%s",buf1);
        else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
        Fhydro[i].open(buf, ios::binary|ios::in);
        header[i].BoxSize = header[ifirstfile].BoxSize;

        //read some of the amr header till get to number of cells in current file
        //@{
        byteoffset=0;
	byteoffset+=RAMSES_fortran_read(Famr[i],header[i].nfiles);  // read ncpu
        byteoffset+=RAMSES_fortran_read(Famr[i],header[i].ndim);
        header[i].twotondim=pow(2,header[i].ndim);
        Famr[i].read((char*)&dummy, sizeof(dummy));
        Famr[i].read((char*)&header[i].nx, sizeof(int));
        Famr[i].read((char*)&header[i].ny, sizeof(int));
        Famr[i].read((char*)&header[i].nz, sizeof(int));
        Famr[i].read((char*)&dummy, sizeof(dummy));
        byteoffset+=RAMSES_fortran_read(Famr[i],header[i].nlevelmax);
        byteoffset+=RAMSES_fortran_read(Famr[i],header[i].ngridmax);
        byteoffset+=RAMSES_fortran_read(Famr[i],header[i].nboundary); // nboundary == 0 i periodic
        byteoffset+=RAMSES_fortran_read(Famr[i],header[i].npart[RAMSESGASTYPE]);

        //Skip until taill (included)
        RAMSES_fortran_skip(Famr[i], 14);
        if (lmin>header[i].nlevelmax) lmin=header[i].nlevelmax;
        if (lmax<header[i].nlevelmax) lmax=header[i].nlevelmax;
        //@}
        //read header info from hydro files
        //@{
        RAMSES_fortran_skip(Fhydro[i]);  // skip ncpu
        RAMSES_fortran_read(Fhydro[i],header[i].nvarh);
        RAMSES_fortran_skip(Fhydro[i]);  // skip ndim
        RAMSES_fortran_skip(Fhydro[i]);  // skip nlevelmax
        RAMSES_fortran_skip(Fhydro[i]);  // skip nboundary
        RAMSES_fortran_read(Fhydro[i],header[i].gamma_index);
        //@}
    }


    // Define some convenience variables
    int ncpu, nboundary;    
    for (i=0;i<opt.num_files;i++) if (ireadfile[i]) {
	ncpu = header[i].nfiles;
	nboundary = header[i].nboundary;
	// First, deal with AMR files
	ncache=new int[(ncpu+nboundary)*header[i].nlevelmax];
	// ncache should contains numbl and numbb

	
	// Read number of grids in level: numbl(1:ncpu, 1:nlevelmax)
	numbl=new int[ncpu*header[i].nlevelmax];
	RAMSES_fortran_read(Famr[i],numbl);

	// Fill the first part of ncache (ibound <= ncpu)
	for (j=0;j<header[i].nlevelmax;j++) {
	    for (k=0;k<ncpu;k++) {
		ncache[(ncpu+nboundary)*j + k] = numbl[ncpu*j+k];
	    }
	}
	
	// skip total number of grids per level: numbtot(1:10, 1:nlevelmax)
	RAMSES_fortran_skip(Famr[i]);
	if (header[i].nboundary>0) {
	    // simple_boundary case, non periodic
	    // skip headb, tailb
	    RAMSES_fortran_skip(Famr[i], 2);
	    // read numbb
	    numbb=new int[header[i].nboundary*header[i].nlevelmax];
	    RAMSES_fortran_read(Famr[i],numbb);
	    // If needed, fill the rest of ncache
	    for (j=0;j<header[i].nlevelmax;j++) {
		for (k=0;k<nboundary;k++) {
		    ncache[(ncpu+nboundary)*j + k+ncpu] = numbl[nboundary*j+k];
		}
	    }
	}
	// skip free memory (headf,tailf,numbf,used_mem,used_mem_tot) and ordering
	RAMSES_fortran_skip(Famr[i], 2);
	// skip keys
	if (orderingstring==string("bisection")) RAMSES_fortran_skip(Famr[i],5);
	else RAMSES_fortran_skip(Famr[i]);
	// skip coarse levels (son, flag1, cpu_map)
	RAMSES_fortran_skip(Famr[i], 3);

        ninputoffset=0;

	// start loop on levels for hydro and AMR
	for (j=0;j<header[i].nlevelmax;j++) {
	    double dx = pow(0.5, j+1); // local cell size
	    // start loop on boundaries for both
	    for (k=0;k<(ncpu+nboundary);k++) {
		/* AMR structure:
		   --------------
		   (only if ncache>0)
		   integer: ind_grid(1..ncache)
		   integer: next(1..ncache)
		   integer: prev(1..ncache)
		   for idim in [1,ndim]:
		     float: xg(1..ncache, idim)
		   integer: father(1..ncache)
		   for ind in [1,2*ndim]:
		     integer: nbor(1..ncache, ind)
		   for ind in [1,2**ndim]:
		     integer: son(1..ncache, ind)  // with iskip
		   for ind in [1,2**ndim]:
		     integer: cpu_map(1..ncache, ind)  // with iskip
		   for ind in [1,2**ndim]:
		     integer: flag1(1..ncache, ind)  // with iskip

		   Hydro structure:
		   ----------------
		   integer: ilevel
		   integer: ncache
		   if ncache>0:
		     for ind in [1,2**ndim]:
		       iskip <- ncoarse + ind*ngridmax
		       write density[ind_grid(i)+iskip]
		       write vx[ind_grid(i)+iskip]
		       ...
		   
		*/

		// get ncache
		chunksize = ncache[(ncpu+nboundary)*j + k];
		if (chunksize>0) {
		    // We want the cell positions from the AMR files
		    // Skip ind_grid, next and prev
		    RAMSES_fortran_skip(Famr[i], 3);
		    // Store grid centres
		    xtempchunk=new RAMSESFLOAT[3*chunksize];
		    for (idim=0;idim<header[i].ndim;idim++) {
                        RAMSES_fortran_read(Famr[i],&xtempchunk[idim*chunksize]);
                    }
		    // Skip father (1) and nbor (2*ndim)
		    RAMSES_fortran_skip(Famr[i], 1+2*header[i].ndim);
		    // Read son index (2**ndim), needed to identify leaf cells
		    icellchunk=new int[header[i].twotondim*chunksize];
		    for (idim=0;idim<header[i].twotondim;idim++) {
                        RAMSES_fortran_read(Famr[i],&icellchunk[idim*chunksize]);
                    }
                    //skip cpu map and refinement map (2**ndim * 2)
                    RAMSES_fortran_skip(Famr[i],2*header[i].twotondim);
		} // chunksize > 0

		// We can now start working with the hydro files
		// Skip ilevel and ncache (already known)
		RAMSES_fortran_skip(Fhydro[i], 2);
		
		if (chunksize>0) {
		    // Define an array to receive hydro data: size is (ncache)*(2**ndim)*(nvar)
		    hydrotempchunk=new RAMSESFLOAT[chunksize*header[i].twotondim*header[i].nvarh];
		    // loop over cells in octs (per dimension)
                    for (idim=0;idim<header[i].twotondim;idim++) {
			// loop over variables
                        for (ivar=0;ivar<header[i].nvarh;ivar++) {
			    // Read a block of size ncache, for each variable, for each sub-oct
			    RAMSES_fortran_read(Fhydro[i], &hydrotempchunk[chunksize*(idim*header[i].nvarh+ivar)]);			    
			} // loop over variables
		    } // loop over cells in octs (per dimension)

		    // At this stage, hydrotempchunk should contain the hydro data
		    // for a given level and a given cpu

		    // We can now start to store things in particle structures
                    for (idim=0;idim<header[i].twotondim;idim++) {  // loop over cells in octs (per dimension)
			int ix=0, iy=0, iz=0;
			for (igrid=0;igrid<chunksize;igrid++) {  // loop over cells in the chunk
			    // Select only leaf cells (internal cells or at maximum level)
			    if (icellchunk[idim*chunksize+igrid]==0 || j==header[i].nlevelmax-1) {
				// Deal with positions
				iz = idim/4;                  // z-position of the cell in the oct
				iy = (idim - (4*iz))/2;       // y-position of the cell in the oct
				ix = idim - (2*iy) - (4*iz);  // x-position of the cell in the oct
				// cell centres are xc=(dble(ix)-0.5D0)*dx, same for yc and zc
				// then, x=(xg(ind_grid_new(i),1)+xc(ind,1)-skip_loc(1))*scale
				// here: (xtempchunk[igrid] + (double(ix)-0.5)*dx)*scale
				// Added jitter should be (U(0,1)-0.5)*dx = U(0,1)*dx - 0.5*dx
				// the boxsize is (?) already included in the "lscale" used later
				// pos  = [ --------------- jitter -------------- ] + [ ------ grid centre ------ ] + [ position in oct ]
				xpos[0] = ( (float)rand()/(float)RAND_MAX - 0.5)*dx + xtempchunk[igrid+0*chunksize] + (double(ix)-0.5)*dx ;
				xpos[1] = ( (float)rand()/(float)RAND_MAX - 0.5)*dx + xtempchunk[igrid+1*chunksize] + (double(iy)-0.5)*dx ;
				xpos[2] = ( (float)rand()/(float)RAND_MAX - 0.5)*dx + xtempchunk[igrid+2*chunksize] + (double(iz)-0.5)*dx ;
				
#ifdef USEMPI

				//CHECK: determine processor this particle belongs on based on its spatial position
				ibuf=MPIGetParticlesProcessor(xpos[0],xpos[1],xpos[2]);
				ibufindex=ibuf*BufSize+Nbuf[ibuf];
#endif
				// Deal with velocities
				vpos[0] = hydrotempchunk[chunksize*(idim*header[i].nvarh + IDvx) + igrid];
				vpos[1] = hydrotempchunk[chunksize*(idim*header[i].nvarh + IDvy) + igrid];
				vpos[2] = hydrotempchunk[chunksize*(idim*header[i].nvarh + IDvz) + igrid];

				// Density
				rhotemp = hydrotempchunk[chunksize*(idim*header[i].nvarh + IDrho) + igrid];
				// Mass
				mtemp = rhotemp * dx*dx*dx;

				// Internal energy (P/rho)
				//TODO: check if should be divided by rho
				utemp=hydrotempchunk[chunksize*(idim*header[i].nvarh + IDU) + igrid] / rhotemp / (header[i].gamma_index-1.0);

				// Metallicity
				//TODO: check if should be divided by rho
				Ztemp=hydrotempchunk[chunksize*(idim*header[i].nvarh + IDZ) + igrid];

#ifdef HIGHRES
                                // Refinement scalar
                                if (IDref > 0) reftemp=hydrotempchunk[chunksize*(idim*header[i].nvarh + IDref) + igrid];
#endif

				// Set density scale
				rhotemp = rhotemp*rhoscale;

				// Create particles
// #ifdef HIGHRES
//                                 if (reftemp>=0.5) {  // ARBITRARY THRESHOLD, FOR TEST, SHOULD BE A PARAMETER (EVEN 0?)
// #endif
				if (opt.partsearchtype==PSTALL || opt.partsearchtype==PSTGAS) {
				    // TODO: PSTGAS
#ifdef USEMPI
				    Pbuf[ibufindex]=Particle(mtemp*mscale,
							     xpos[0]*lscale,xpos[1]*lscale,xpos[2]*lscale,
							     vpos[0]*opt.V+Hubbleflow*xpos[0],
							     vpos[1]*opt.V+Hubbleflow*xpos[1],
							     vpos[2]*opt.V+Hubbleflow*xpos[2],
							     count2,GASTYPE);
				    Pbuf[ibufindex].SetPID(gascount + GASTYPE*offsetID);  // GASCOUNT NEEDS TO BE SHARED WITH MPI?
#ifdef GASON
				    Pbuf[ibufindex].SetU(utemp);
				    Pbuf[ibufindex].SetSPHDen(rhotemp);
#ifdef STARON
				    Pbuf[ibufindex].SetZmet(Ztemp);
#endif
#endif
#ifdef EXTRAINPUTINFO
				    if (opt.iextendedoutput)
				    {
					Pbuf[ibufindex].SetInputFileID(i);
					Pbuf[ibufindex].SetInputIndexInFile(idim*chunksize*header[i].nvarh+0*chunksize+igrid+ninputoffset);
				    }
#endif
				    //ensure that store number of particles to be sent to the threads involved with reading snapshot files
				    Nbuf[ibuf]++;
				    MPIAddParticletoAppropriateBuffer(ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, Nlocal, Part.data(), Nreadbuf, Preadbuf);
#else
				    Part[count2]=Particle(mtemp*mscale,
							  xpos[0]*lscale,xpos[1]*lscale,xpos[2]*lscale,
							  vpos[0]*opt.V+Hubbleflow*xpos[0],
							  vpos[1]*opt.V+Hubbleflow*xpos[1],
							  vpos[2]*opt.V+Hubbleflow*xpos[2],
							  count2,GASTYPE);
				    // Stupid counter for the gas particles (idval -> gascount)
				    Part[count2].SetPID(gascount + GASTYPE*offsetID);
#ifdef GASON
				    Part[count2].SetU(utemp);
				    Part[count2].SetSPHDen(rhotemp);
#ifdef STARON
				    Part[count2].SetZmet(Ztemp);
#endif
#endif
#ifdef EXTRAINPUTINFO
				    if (opt.iextendedoutput)
				    {
					Part[count2].SetInputFileID(i);
					Part[count2].SetInputIndexInFile(idim*chunksize*header[i].nvarh+0*chunksize+igrid+ninputoffset);
				    }
#endif
				    
#endif
				    count2++;
				    gascount++;
				}
				else if (opt.partsearchtype==PSTDARK&&opt.iBaryonSearch) {
#ifdef USEMPI
				    Pbuf[ibufindex]=Particle(mtemp*mscale,
							     xpos[0]*lscale,xpos[1]*lscale,xpos[2]*lscale,
							     vpos[0]*opt.V+Hubbleflow*xpos[0],
							     vpos[1]*opt.V+Hubbleflow*xpos[1],
							     vpos[2]*opt.V+Hubbleflow*xpos[2],
							     bcount2,GASTYPE); //count2 -> bcount2
				    Pbuf[ibufindex].SetPID(gascount + GASTYPE*offsetID);  // GASCOUNT NEEDS TO BE SHARED WITH MPI?
#ifdef GASON
				    Pbuf[ibufindex].SetU(utemp);
				    Pbuf[ibufindex].SetSPHDen(rhotemp);
#ifdef STARON
				    Pbuf[ibufindex].SetZmet(Ztemp);
#endif
#endif
#ifdef EXTRAINPUTINFO
				    if (opt.iextendedoutput)
				    {
					Pbuf[ibufindex].SetInputFileID(i);
					Pbuf[ibufindex].SetInputIndexInFile(chunksize*(idim*header[i].nvarh + IDrho)+ igrid +ninputoffset);
				    }
#endif
				    //ensure that store number of particles to be sent to the reading threads
				    if (ibuf==ThisTask) {
					Nlocalbaryon[1]++;
				    }
				    MPIAddParticletoAppropriateBuffer(ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, Nlocalbaryon[0], Pbaryons, Nreadbuf, Preadbuf);
#else
				    Pbaryons[bcount2]=Particle(mtemp*mscale,
							       xpos[0]*lscale,xpos[1]*lscale,xpos[2]*lscale,
							       vpos[0]*opt.V+Hubbleflow*xpos[0],
							       vpos[1]*opt.V+Hubbleflow*xpos[1],
							       vpos[2]*opt.V+Hubbleflow*xpos[2],
							       bcount2,GASTYPE);  // count2 -> bcount2
				    // Stupid counter for the gas particles (idval -> gascount + GASTYPE*offsetID)
				    Pbaryons[bcount2].SetPID(gascount + GASTYPE*offsetID);
#ifdef GASON
				    Pbaryons[bcount2].SetU(utemp);
				    Pbaryons[bcount2].SetSPHDen(rhotemp);
#ifdef STARON
				    Pbaryons[bcount2].SetZmet(Ztemp);
#endif
#endif
#ifdef EXTRAINPUTINFO
				    if (opt.iextendedoutput)
				    {
					Pbaryons[bcount2].SetInputFileID(i);
					Pbaryons[bcount2].SetInputIndexInFile(chunksize*(idim*header[i].nvarh + IDrho) + igrid +ninputoffset);
				    }
#endif
				    
#endif
				    bcount2++;
				    gascount++;
				}
// #ifdef HIGHRES
//                                 }
// #endif

			    } // leaf cells
			} // loop over cells in the chunk
		    } // loop over cells in octs (per dimension)

		} // chunksize > 0
		
		if (chunksize>0) {
                    delete[] xtempchunk;
                    delete[] hydrotempchunk;
		    delete[] icellchunk;
		}
		
	    } // loop on ncpu+nboudary
	} // loop on levels
	Famr[i].close();
	Fhydro[i].close();

#ifdef USEMPI
        //send information between read threads
        if (opt.nsnapread>1&&inreadsend<totreadsend){
            MPI_Allgather(Nreadbuf, opt.nsnapread, MPI_Int_t, mpi_nsend_readthread, opt.nsnapread, MPI_Int_t, mpi_comm_read);
            MPISendParticlesBetweenReadThreads(opt, Preadbuf, Part.data(), ireadtask, readtaskID, Pbaryons, mpi_comm_read, mpi_nsend_readthread, mpi_nsend_readthread_baryon);
            inreadsend++;
            for(ibuf = 0; ibuf < opt.nsnapread; ibuf++) Nreadbuf[ibuf]=0;
        }
#endif
	} // loop on files, if read
    
#ifdef USEMPI
    //once finished reading the file if there are any particles left in the buffer broadcast them
    for(ibuf = 0; ibuf < NProcs; ibuf++) if (ireadtask[ibuf]<0)
    {
        MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
        if (Nbuf[ibuf]>0) {
            MPI_Ssend(&Pbuf[ibuf*BufSize], sizeof(Particle)*Nbuf[ibuf], MPI_BYTE, ibuf, ibuf, MPI_COMM_WORLD);
            Nbuf[ibuf]=0;
            //last broadcast with Nbuf[ibuf]=0 so that receiver knows no more particles are to be broadcast
            MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t,ibuf,ibuf+NProcs,MPI_COMM_WORLD);
        }
    }
    if (opt.nsnapread>1){
        MPI_Allgather(Nreadbuf, opt.nsnapread, MPI_Int_t, mpi_nsend_readthread, opt.nsnapread, MPI_Int_t, mpi_comm_read);
        MPISendParticlesBetweenReadThreads(opt, Preadbuf, Part.data(), ireadtask, readtaskID, Pbaryons, mpi_comm_read, mpi_nsend_readthread, mpi_nsend_readthread_baryon);
    }
    }//end of reading task
#endif
#ifdef USEMPI
    //if not reading information than waiting to receive information
    else {
        MPIReceiveParticlesFromReadThreads(opt,Pbuf,Part.data(),readtaskID, irecv, mpi_irecvflag, Nlocalthreadbuf, mpi_request,Pbaryons);
    }
#endif
    }//end of check if gas loaded

    //update info
    //in gadgetio.cxx: opt.a/opt.h if not comove, otherwise /opt.h
    // by default, opt.p is boxsize, expressed in kpc/h
    // FIXME: check what should be used.
    opt.p*=opt.a/opt.h;

    // By default, LN should be based on levelmin (?)
    LN   = (lscale/pow(2.0, header[ifirstfile].levelmin));
#ifdef HIGHRES
    opt.zoomlowmassdm=MP_DM*mscale * (1.0001);  // Add small offset ("a la HaloMaker")
    cout<<"Lowest DM particle mass: "<<opt.zoomlowmassdm<<" Msun"<<endl;

    if (opt.Neff==-1) {
	if  (opt.partsearchtype==PSTDARK||opt.partsearchtype==PSTALL) {
	    // Ideally, we would want this to be the "min level of the high res region"
	    // This can be inferred from (total_DM_mass/smallest_DM_mass)**(1./3.)
	    // luckily, this is (1./MP_DM)**(1./3.) in code units, with a factor Om/Ocdm
	    LN = lscale * pow(MP_DM*opt.Omega_m/opt.Omega_cdm, 1./3.);
	}
	else {
	    // not looking for DM particles: we use levelmax information
	    LN   = (lscale/pow(2.0, header[ifirstfile].nlevelmax));
	}
    }
    else {
	LN   = (lscale/(double)opt.Neff);
    }
#endif
    opt.ellxscale = LN;
    opt.uinfo.eps*=LN;

#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
    //update cosmological data and boundary in code units
    MPI_Bcast(&(opt.p),sizeof(opt.p),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.a),sizeof(opt.a),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.Omega_cdm),sizeof(opt.Omega_cdm),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.Omega_b),sizeof(opt.Omega_b),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.Omega_m),sizeof(opt.Omega_m),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.Omega_Lambda),sizeof(opt.Omega_Lambda),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.h),sizeof(opt.h),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.rhocrit),sizeof(opt.rhocrit),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.rhobg),sizeof(opt.rhobg),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.virlevel),sizeof(opt.virlevel),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.virBN98),sizeof(opt.virBN98),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.ellxscale),sizeof(opt.ellxscale),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.lengthinputconversion),sizeof(opt.lengthinputconversion),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.velocityinputconversion),sizeof(opt.velocityinputconversion),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.massinputconversion),sizeof(opt.massinputconversion),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.G),sizeof(opt.G),MPI_BYTE,0,MPI_COMM_WORLD);
#ifdef NOMASS
    MPI_Bcast(&(opt.MassValue),sizeof(opt.MassValue),MPI_BYTE,0,MPI_COMM_WORLD);
#endif
    MPI_Bcast(&(Ntotal),sizeof(Ntotal),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&opt.zoomlowmassdm,sizeof(opt.zoomlowmassdm),MPI_BYTE,0,MPI_COMM_WORLD);
#endif

    //a bit of clean up
#ifdef USEMPI
    MPI_Comm_free(&mpi_comm_read);
    if (opt.iBaryonSearch) delete[] mpi_nsend_baryon;
    if (opt.nsnapread>1) {
        delete[] mpi_nsend_readthread;
        if (opt.iBaryonSearch) delete[] mpi_nsend_readthread_baryon;
        if (ireadtask[ThisTask]>=0) delete[] Preadbuf;
    }
    delete[] Nbuf;
    if (ireadtask[ThisTask]>=0) {
        delete[] Nreadbuf;
        delete[] Pbuf;
        delete[] ireadfile;
    }
    delete[] ireadtask;
    delete[] readtaskID;
#endif
}
