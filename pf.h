/*************************************************************************
 > File Name: pf.h
 > Author: Weidong, ZHANG; Chang Cui
 > Mail: zhangwd@pku.edu.cn, professorcui@pku.edu.cn
 > Created Time: Fri Apr 14 19:44:03 2017
 ************************************************************************/

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <memory>
#include <cstdio>
#include <iomanip>
#include <cmath>
#include <string>
#include <cstdarg>
#include <mpi.h>
#include <valarray>

using namespace std;

namespace pf {
// Const definition, which don't appear in configure file
static const int MAXN = 430;
static const int MAXUNBO = MAXN * (MAXN - 1) / 2;
// Statistics: sample number of histogram
static const int nbinmax = 105;
static const int nEbinmax = 105;
static const int nRbinmax = 105;
// Consts used to strength of interactions (调整相互作用的强度参数)
static const double eps = 1.0e-3, eps1 = 1.0 - 1.0e-6, eps2 = 1.0e-4,
      epstht = 1.0-1.0e-12;

class Logger {
private:
    std::map<std::string, std::ofstream> outputfstream;
public:
    Logger() {};
    ~Logger();
    int info(string filename, string msg);
    int debug(string filename, string msg);
    int warn(string filename, string msg);
    int err(string filename, string msg);

    string format(const char* fmt, ...);

private:
    static std::string outputdir;
};

class statistics {
public:
    double PFBbin[nbinmax][nbinmax], PFbin[nbinmax], PBbin[nbinmax];
    double PFBEbin[nbinmax][nbinmax][nEbinmax];
    double PQbEbbin[nbinmax][nEbinmax];
    double PFBRbin[nbinmax][nbinmax][nRbinmax];
    double PFFBbin[nbinmax][nbinmax][nbinmax];
    double PQwbin[nbinmax], PQwQbbin[nbinmax][nbinmax],
           PQwQfbin[nbinmax][nbinmax], PQwRbin[nbinmax][nRbinmax][2],
           PQwEbbin[nbinmax][nEbinmax][2];
} statis;
/* Class Parameter for maintaining parameters
 * main functions: Parameter(), Parse()
 * main members: parameters enumerating in input.***.test.dat
 */
class Parameter {
public:
    Parameter();
    Parameter(string filename);
    int parse(string filename);
    void display();
public:
    // Parameters in the head of configure file
    // conf_filename = input.****.test.dat pdb_filename = pdb
    // initialconform_filename = InitialConf.dat
    // nativeconform_filename = NativeConf.dat
    // filename5 = appNCS.dat
    string  title, date;
    string conf_filename, pdb_filename;
    string nativeconform_filename, initialconform_filename, appNCS_filename;
    string term, errmess;
    // common/Bnum/npart1, iFlagMov, nResGap, npartM, nparttol, nunbond
    int npart1, iFlagMov, nResGap, npartM, nparttol, nunbond;


    // BUG: maybe some int is float, some float is int
    // Parameters appear in configure file
    string rand_num;
    int ioctsd[4];
    double epsil, epsil1, epsil2, enscale;
    double ck_r, ck_tht, ck_phi1, ck_phi3;
    double sigma_ij, amass, gamma;
    double avsn0, boltz;
    double gr0, gdr;
    double s_dt, s_gm; int  iFixKr;
    int k_sol, n_sol, m_sol;
    double epsilon_p, epsilon_pp;
    double temp;
    double ga1_f1, ga2_f1, gQ0_f1;
    double ga1_f2, ga2_f2, gQ0_f2;
    double ga1_b, ga2_b, gQ0_b;
    double ga1_w, ga2_w, gQ0_w, alpha_Qw;
    int Is_Solvation;
    int nConform, nCon0, nRunConf;
    double gQbnativ, gQbdenatural, gQf1nativ, gQf1denatural, gQf2nativ, 
           gQf2denatural;
    int nsnap, nstep;
    int nConformOutput, nOutput0, ndOutput;
    double outQf1_i, outQf1_f, outQf2_i, outQf2_f, outQb_i, outQb_f;
    int nbinsnap, nbinsnap0;
    int nbin_f, nbin_b;
    double dbin_f, dbin_b, vbin0;
    int IsEbin, nEbin; double dEbin, vEbin0;
    int IsEbbin, nEbbin; double dEbbin, vEbbin0;
    int IsRbin; double nRbin, dRbin, vRbin0;
    int IsWbin, nWbin; double dWbin, vWbin0, cri_Qb;
    double pL, dl;
    double Alpha1, Alpha2, Beta;
    double Delta, CritR_non;

    // The following parameters don't appear in configure file. I'm not sure 
    // whether parts of these arguments are initialized in "block data"
    double gyr;
    double pi;
    int nadim;
    double dt, gm;
    int iun[MAXUNBO], jun[MAXUNBO], kunbond[MAXUNBO], nQnative_f1, nQnative_f2,
        nQnative_b;
    double c_0, c_1, c_2, randconst;
    double gQ_f, gQ_f1, gQ_f2, gQ_b, eGr;
    double gQ_w;
    double iQb_i, iQb_f;
    double pseudoQ_f, pseudoQ_b, dr_sol, ddr_sol;
    double R;
    double gQ_non_f, gQ_non_b;

    // bond length, non-bond length, bond angle, dihedral angle
    // (键长/键角/二面角和非键距离)
    // Arrays keep unchanged after execute nativeinformation()
    double rbond_nat[MAXN], runbond_nat[MAXUNBO], theta_nat[MAXN], \
        dihedral_nat[MAXN];

}; // end of class Parameter

/* Class Particle for solving force fields 
 * main functions: Force(), e_bond_tot()
 */
class Particle {
public:
    Particle(double x, double y, double z);

public:
    // Arguments which don't appear in configure file
    // double xinit[MAXN], yinit[MAXN], zinit[MAXN];
    // double xsave[MAXN], ysave[MAXN], zsave[MAXN];
    double xinit, yinit, zinit;
    double xsave, ysave, zsave;

    // Particle's attributes which are prefixed with common
    // common/coordinates/x(MAXN),y(MAXN),z(MAXN)
    // common/velocity/vx(MAXN),vy(MAXN),vz(MAXN)
    // common/forcetotalnew/fx(MAXN),fy(MAXN),fz(MAXN)
    // common/forcetotalold/fxo(MAXN),fyo(MAXN),fzo(MAXN)
    // common/forcestretching/fxr(MAXN),fyr(MAXN),fzr(MAXN)
    // common/forcebending/fxth(MAXN),fyth(MAXN),fzth(MAXN)
    // common/forcetorsion/fxph(MAXN),fyph(MAXN),fzph(MAXN)
    // common/forceunbonded/fxun(MAXN),fyun(MAXN),fzun(MAXN)
    // common/forcerandnew/frandx(MAXN),frandy(MAXN),frandz(MAXN)
    // common/forcerandold/frandxo(MAXN),frandyo(MAXN),frandzo(MAXN)
    // common/nativeinfo/rbond_nat(MAXN),runbond_nat(MAXUNBO),theta_nat(MAXN),
    // dihedral_nat(MAXN)
    double x, y, z;
    double vx, vy, vz;
    double fx, fy, fz;
    // *o: *old, for the f_{k + 1} = f_{k} + f_{k - 1}
    double fxo, fyo, fzo;
    // *r, *th, *ph, *un: rbond, theta, phi, unbond
    double fxr, fyr, fzr;
    double fxth, fyth, fzth;
    double fxph, fyph, fzph;
    double fxun, fyun, fzun;
    // *rand*: force by Brownian Motion (布朗运动产生的随机力)
    double frandx, frandy, frandz;
    double frandxo, frandyo, frandzo;

}; // end of class Particle

/* Class Force for solving force fields 
 * main functions: Force(), e_bond_tot()
 */
class Force {

public:
    // composite force
    Force(Parameter &param);

    int force(vector<Particle> &particle_list, double &e_pot,
            double &e_unbond_tot, double &e_bind_tot, double &e_tors_tot,
            double &e_bend_tot, double &e_bond_tot, int start_idx, int end_idx);
    // component forces
    double fbond(vector<Particle> &particle_list, double &e_bond_tot,
            int start_idx, int end_idx);
    double fbend(vector<Particle> &particle_list, double &e_bend_tot,
            int start_idx, int end_idx);
    double ftorsion(vector<Particle> &particle_list, double &e_tors_tot,
            int start_idx, int end_idx);
    double funbond(vector<Particle> &particle_list, double &e_unbond_tot,
            double &e_bind_tot, int start_idx, int end_idx);

    double funbond_with(vector<Particle> &particle_list, double &e_unbond_tot,
            double &e_bind_tot, int start_idx, int end_idx);
    double funbond_without(vector<Particle> &particle_list, double
            &e_unbond_tot, double &e_bind_tot, int start_idx, int end_idx);

public:
    Parameter &param;

}; // end of class Force

/* Class Simulation for simulating protein folding process 
 * main functions: Simulation(), start()
 */
class Simulation {

public:
    Simulation(string conf_filename);

    // Initialize parameter reference from fortran
    int intpar(double &enerkin);
    int start_simulation(int rank, int size);

    // Init particle list
    int init_particle(string filename);

    // Generate output filenames
    // Output the intermediate status
    int formate_output_filenames();
    int open_ir_file();

    // Output the Simulation information
    int output_info();

    // Read particles from file
    // Read appNCS_****.dat
    int read_nativeconform(string filename);
    int read_appNCS(string filename);
    int read_initalconform(string filename);

    // Out of outermost loop
    int nativeinformation();
    // In outermost loop
    // Algin conformation to coordinate system
    int origin_adjust();
    int InitVel(double &enerkin);

    // Other functions
    // Related to Random generator
    int RANTERM();
    // verlet: one sort of dynamics model
    int verlet(double &enerkin, double &e_pot, int start_idx, int end_idx);
    int pbc_shift();
    int output_conformation(int &nOutputGap, int &nOutputCount);
    int write_mol(int &nOutputCount);
    int write_target(int &nOutputCount);
    int write_histogram();

    // mathematical related function
    // TL: When to execute these functions: rnyubd, savern is unused
    double gauss(double xsi);
    // Initialize the m(4) and l(4)
    int rnyubd();
    double rannyu();
    int setrn(int iseed[4]);
    int savern(int iseed[4]);

    // if mpi is enabled
    int exchange_particle(int rank, int size);

    // Unused function
    double calculate_gyrationradius();

public:
    vector<Particle> particle_list;
    Parameter param;
    Force force;
    Logger log;
    // Output filenames
    vector<string> output_filenames;
    // kinetic energy
    double enerkin;
    // TL: what's the mean of xsi, a1-9
    double xsi = 0.0;
    const double pi = 3.141592654, a1 = 3.949846138, a3 = 0.252408784,
          a5 = 0.076542912, a7 = 0.008355968, a9 = 0.029899776;
    // TL: when to initialize the m(4), while l(4) is init in setrn()
    int m[4], l[4];

    // Statistic and profiling
    double start_time;

    static string inputdir;

}; // end of class Simulation

/*
 * Class MPI_Core for supporting MPI Utility
 */
class MPI_Util {
public:
    MPI_Util(int argc, char** argv);
    ~MPI_Util();
    bool is_available();
    int size();
    int rank();

private:
    bool mpi_available;
    int world_size;
    int world_rank;
}; // end of class MPI_Core

/*
 * Class CUDA_Core for supporting CUDA Utility
 */
class CUDA_Util {
public:
    bool init();
    int rank();
    bool finish();

}; // end of class CUDA_Core
} // end of namespace::pf
