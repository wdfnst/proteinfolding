/*************************************************************************
 > File Name: pf.cpp
 > Author: Weidong, ZHANG; Chang Cui
 > Mail: zhangwd@pku.edu.cn; professorcui@pku.edu.cn
 > Created Time: Fri Apr 14 19:43:04 2017
 > For Protein Folding
 ************************************************************************/

#include "pf.h"

using namespace std;

string pf::Logger::format(const char* fmt, ...){
    int size = 512;
    char* buffer = 0;
    buffer = new char[size];
    va_list vl;
    va_start(vl, fmt);
    int nsize = vsnprintf(buffer, size, fmt, vl);
    if(size<=nsize){ //fail delete buffer and try again
        delete[] buffer;
        buffer = 0;
        buffer = new char[nsize+1]; //+1 for /0
        nsize = vsnprintf(buffer, size, fmt, vl);
    }
    std::string ret(buffer);
    va_end(vl);
    delete[] buffer;
    return ret;
}

int pf::Logger::info(string filename, string msg){
    return 0;
}

// Constructors of class Parameter
pf::Parameter::Parameter() {}
pf::Parameter::Parameter(string filename) {
    parse(filename);
}

/////////////////////////////////////////////////////////////////////////
/*******************Definition of class Parameter members***************/
/////////////////////////////////////////////////////////////////////////
// Parse the configure file
int pf::Parameter::parse(string filename) {
    // Open file for reading
    fstream fin(filename, std::ifstream::in);
    if (!fin.is_open()) {
        cerr << "failed to open " << filename << '\n';
        return 1;
    }

    // Read file head
    fin >> title >> date;
    fin >> pdb_filename;
    fin >> npart1 >> nparttol >> nResGap >> iFlagMov;
    fin >> initialconform_filename;
    fin >> nativeconform_filename;
    fin >> appNCS_filename;
    if (0 == iFlagMov || 1 == iFlagMov) {
        npartM = npart1;
    }
    else {
        npartM = nparttol;
    }
    if (npart1 > nparttol) cerr << "Err: npart1 > nparttol\n";
    if (npart1 < 1 || npart1 > MAXN) cerr << "Err: npart1 < 1 or npart1 > MAXN\n";
    if (nparttol > MAXN) cerr << "Err: nparttol > MAXN\n";

    // Read separator between file head and file body
    fin >> term;
    // Read file body of configure file
    string comments;
    fin >> rand_num;
    // Split the rand_num to 4 integer
    sscanf(rand_num.c_str(), "%4d%4d%4d%4d", &ioctsd[0], &ioctsd[1], 
            &ioctsd[2], &ioctsd[3]);
    // Skip the end of the comments
    getline(fin, comments);

    fin >> epsil >> epsil1 >> epsil2 >> enscale;
    getline(fin, comments);
    fin >> ck_r >> ck_tht >> ck_phi1 >> ck_phi3;
    getline(fin, comments);
    fin >> sigma_ij >> amass >> gamma;
    getline(fin, comments);
    fin >> avsn0 >> boltz;
    getline(fin, comments);
    fin >> gr0 >> gdr;
    getline(fin, comments);
    fin >> s_dt >> s_gm >> iFixKr;
    getline(fin, comments);
    fin >> k_sol >> n_sol >> m_sol;
    getline(fin, comments);
    fin >> epsilon_p >> epsilon_pp;
    getline(fin, comments);
    fin >> temp;
    getline(fin, comments);
    fin >> ga1_f1 >> ga2_f1 >> gQ0_f1;
    getline(fin, comments);
    fin >> ga1_f2 >> ga2_f2 >> gQ0_f2;
    getline(fin, comments);
    fin >> ga1_b >> ga2_b >> gQ0_b;
    getline(fin, comments);
    fin >> ga1_w >> ga2_w >> gQ0_w >> alpha_Qw;
    getline(fin, comments);
    fin >> Is_Solvation;
    getline(fin, comments);
    fin >> nConform >> nCon0 >> nRunConf;
    getline(fin, comments);
    fin >> gQbnativ >> gQbdenatural >> gQf1nativ >> gQf1denatural >> gQf2nativ 
        >> gQf2denatural;
    getline(fin, comments);
    fin >> nsnap >> nstep;
    getline(fin, comments);
    fin >> nConformOutput >> nOutput0 >> ndOutput;
    getline(fin, comments);
    fin >> outQf1_i >> outQf1_f >> outQf2_i >> outQf2_f >> outQb_i >> outQb_f;
    getline(fin, comments);
    fin >> nbinsnap >> nbinsnap0;
    getline(fin, comments);
    fin >> nbin_f >> nbin_b;
    getline(fin, comments);
    fin >> dbin_f >> dbin_b >> vbin0;
    getline(fin, comments);
    fin >> IsEbin >> nEbin >> dEbin >> vEbin0;
    getline(fin, comments);
    fin >> IsEbbin >> nEbbin >> dEbbin >> vEbbin0;
    getline(fin, comments);
    fin >> IsRbin >> nRbin >> dRbin >> vRbin0;
    getline(fin, comments);
    fin >> IsWbin >> nWbin >> dWbin >> vWbin0 >> cri_Qb;
    getline(fin, comments);
    fin >> pL >> dl;
    getline(fin, comments);
    fin >> Alpha1 >> Alpha2 >> Beta;
    getline(fin, comments);
    fin >> Delta >> CritR_non;
    getline(fin, comments);

    // Reasonability check
    if(nbin_f > nbinmax) cerr << "Err: nbin_f > nbinmax\n";
    if(nbin_b > nbinmax) cerr << "Err: nbin_b > nbinmax\n";
    if(nEbin > nEbinmax) cerr << "Err: nEbin > nEbinmax\n";
    if(nRbin > nRbinmax) cerr << "Err: nRbin > nRbinmax\n";

    fin.close();
    return 0;
}

void pf::Parameter::display() {
    // file header 
    cout << title << '\t' << date << endl;
    cout << pdb_filename << endl;
    cout << npart1 << '\t' << nparttol << '\t' << nResGap << '\t' << iFlagMov 
        << endl;
    cout << initialconform_filename << endl;
    cout << nativeconform_filename << endl;
    cout << appNCS_filename << endl;

    // separator
    cout << term << endl;

    // file body 
    cout << '\t' << rand_num << '\t' << endl;
    cout << '\t' << ioctsd[0] << ' ' << ioctsd[1] << ' ' << ioctsd[2] << ' ' <<
        ioctsd[3] << '\t' << endl;
    cout << '\t' << epsil << '\t' << epsil1 << '\t' << epsil2 << '\t' << 
        enscale << endl;
    cout << '\t' << ck_r << '\t' << ck_tht << '\t' << ck_phi1 << '\t' <<
        ck_phi3 << endl;
    cout << '\t' << sigma_ij << '\t' << amass << '\t' << gamma << endl;
    cout << '\t' << avsn0 << '\t' << boltz << endl;
    cout << '\t' << gr0 << '\t' << gdr << endl;
    cout << '\t' << s_dt << '\t' << s_gm << '\t' << iFixKr << endl;
    cout << '\t' << k_sol << '\t' << n_sol << '\t' << m_sol << endl;
    cout << '\t' << epsilon_p << '\t' << epsilon_pp << endl;
    cout << '\t' << temp << endl;
    cout << '\t' << ga1_f1 << '\t' << ga2_f1 << '\t' << gQ0_f1 << endl;
    cout << '\t' << ga1_f2 << '\t' << ga2_f2 << '\t' << gQ0_f2 << endl;
    cout << '\t' << ga1_b << '\t' << ga2_b << '\t' << gQ0_b << endl;
    cout << '\t' << ga1_w << '\t' << ga2_w << '\t' << gQ0_w << '\t' << alpha_Qw
        << endl;
    cout << '\t' << Is_Solvation << endl;
    cout << '\t' << nConform << '\t' << nCon0 << '\t' << nRunConf << endl;
    cout << '\t' << gQbnativ << '\t' << gQbdenatural << '\t' << gQf1nativ << 
        '\t' << gQf1denatural << '\t' << gQf2nativ << '\t' << gQf2denatural <<
        endl;
    cout << '\t' << nsnap << '\t' << nstep << endl;
    cout << '\t' << nConformOutput << '\t' << nOutput0 << '\t' << ndOutput <<
        endl;
    cout << '\t' << outQf1_i << '\t' << outQf1_f << '\t' << outQf2_i << '\t' <<
        outQf2_f << '\t' << outQb_i << '\t' << outQb_f << endl;
    cout << '\t' << nbinsnap << '\t' << nbinsnap0 << endl;
    cout << '\t' << nbin_f << '\t' << nbin_b << endl;
    cout << '\t' << dbin_f << '\t' << dbin_b << '\t' << vbin0 << endl;
    cout << '\t' << IsEbin << '\t' << nEbin << '\t' << dEbin << '\t' << vEbin0
        << endl;
    cout << '\t' << IsEbbin << '\t' << nEbbin << '\t' << dEbbin << '\t' <<
        vEbbin0 << endl;
    cout << '\t' << IsRbin << '\t' << nRbin << '\t' << dRbin << '\t' << vRbin0
        << endl;
    cout << '\t' << IsWbin << '\t' << nWbin << '\t' << dWbin << '\t' << vWbin0
        << '\t' << cri_Qb << endl;
    cout << '\t' << pL << '\t' << dl << endl;
    cout << '\t' << Alpha1 << '\t' << Alpha2 << '\t' << Beta << endl;
    cout << '\t' << Delta << '\t' << CritR_non << endl;
}
/////////////////////////////////////////////////////////////////////////
/******************End of definition of class Parameter*****************/
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
/**********************Definition of class Particle*********************/
/////////////////////////////////////////////////////////////////////////
pf::Particle::Particle(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
}
/////////////////////////////////////////////////////////////////////////
/******************End of definition of class Particle******************/
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
/*************************Definition of class Force*********************/
/////////////////////////////////////////////////////////////////////////

pf::Force::Force(Parameter &param_) :param(param_){ }

int pf::Force::force(vector<Particle> &particle_list, double &e_pot, 
        double &e_unbond_tot, double &e_bind_tot, double &e_tors_tot, 
        double &e_bend_tot, double &e_bond_tot) {
    e_pot=0.0;

    for (int i = 0; i < param.npartM; i++) {
        particle_list[i].fx = 0;
        particle_list[i].fy = 0;
        particle_list[i].fz = 0;
        particle_list[i].fxr = 0;
        particle_list[i].fyr = 0;
        particle_list[i].fzr = 0;
        particle_list[i].fxth = 0;
        particle_list[i].fyth = 0;
        particle_list[i].fzth = 0;
        particle_list[i].fxph = 0;
        particle_list[i].fyph = 0;
        particle_list[i].fzph = 0;
        particle_list[i].fxun = 0;
        particle_list[i].fyun = 0;
        particle_list[i].fzun = 0;
    }

    e_bond_tot = 0.0;
    e_bend_tot = 0.0;
    e_tors_tot = 0.0;
    e_unbond_tot = 0.0;
    e_bind_tot = 0.0;

    fbond(particle_list, e_bond_tot);
    fbend(particle_list, e_bend_tot);
    ftorsion(particle_list, e_tors_tot);
    funbond(particle_list, e_unbond_tot, e_bind_tot);

    for (int i = 0; i < param.npartM; i++) {
        particle_list[i].fx += (particle_list[i].fxr + particle_list[i].fxth +
            particle_list[i].fxph + particle_list[i].fxun);
        particle_list[i].fy += (particle_list[i].fyr + particle_list[i].fyth +
                particle_list[i].fyph + particle_list[i].fyun);
        particle_list[i].fz += (particle_list[i].fzr + particle_list[i].fzth +
                particle_list[i].fzph + particle_list[i].fzun);
    }
    e_pot = e_bond_tot + e_bend_tot + e_tors_tot + e_unbond_tot;
    return 0;
}

double pf::Force::fbond(vector<Particle> &particle_list, double &e_bond_tot) {
    for (int i = 0; i < param.npartM - 1; i++) {
        if (i == param.npart1) break;
        int j = i + 1;
        double xij = particle_list[i].x - particle_list[j].x;
        double yij = particle_list[i].y - particle_list[j].y;
        double zij = particle_list[i].z - particle_list[j].z;
        double rij = sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2));

        if (rij < eps1) break;

        double e_bond = param.ck_r * pow(rij - param.rbond_nat[i], 2);
        e_bond_tot += e_bond;
        double e_bond_drv =
            -2 * param.ck_r * (rij - param.rbond_nat[i]);
        double drix = xij / rij;
        double driy = yij / rij;
        double driz = zij / rij;
        double drjx = -drix;
        double drjy = -driy;
        double drjz = -driz;

        if (i != param.npart1) {
            particle_list[i].fxr += e_bond_drv * drix;
            particle_list[i].fyr += e_bond_drv * driy;
            particle_list[i].fzr += e_bond_drv * driz;

            particle_list[j].fxr += e_bond_drv * drjx;
            particle_list[j].fyr += e_bond_drv * drjy;
            particle_list[j].fzr += e_bond_drv * drjz;
        }
    }

    return 0;
}

double pf::Force::fbend(vector<Particle> &particle_list, double &e_bend_tot) {
    for (int i = 0; i < param.npartM - 2; i++) {
        int j = i + 1;
        int k = i + 2;
        double xij = particle_list[i].x - particle_list[j].x;
        double yij = particle_list[i].y - particle_list[j].y;
        double zij = particle_list[i].z - particle_list[j].z;

        double xkj = particle_list[k].x - particle_list[j].x;
        double ykj = particle_list[k].y - particle_list[j].y;
        double zkj = particle_list[k].z - particle_list[j].z;

        double rij = sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2));
        double rkj = sqrt(pow(xkj, 2) + pow(ykj, 2) + pow(zkj, 2));

        if (rij < eps || rkj < eps) continue;
        double costheta = (xij * xkj + yij * ykj + zij * zkj) / (rij * rkj);
        // fortran: if(abs(costheta).gt.epstht) costheta=sign(epstht,costheta)
        if(fabs(costheta) > epstht) {
            costheta = costheta > 0 ? epstht : (-1 * epstht);
        }

        double theta = acos(costheta);
        double delta = theta - param.theta_nat[i];


        double tem = 0.0;
        if (i <= param.npart1 - 2) {
            tem = param.Alpha1;
        } else if (i > param.npart1) {
            tem = param.Alpha2;
        }

        double e_bend = tem * param.ck_tht * pow(delta, 2);
        double e_bend_drv = -1 * tem * 2 * param.ck_tht * delta;

        e_bend_tot += e_bend;
        double sintinv = 1.0 / sin(theta);

        double drix = sintinv * (costheta * xij / rij - xkj / rkj) / rij;
        double driy = sintinv * (costheta * yij / rij - ykj / rkj) / rij;
        double driz = sintinv * (costheta * zij / rij - zkj / rkj) / rij;

        double drkx = sintinv * (costheta * xkj / rkj - xij / rij) / rkj;
        double drky = sintinv * (costheta * ykj / rkj - yij / rij) / rkj;
        double drkz = sintinv * (costheta * zkj / rkj - zij / rij) / rkj;

        particle_list[i].fxth += e_bend_drv * drix;
        particle_list[i].fyth += e_bend_drv * driy;
        particle_list[i].fzth += e_bend_drv * driz;

        // fortran code: fxth(j)=fxth(j)+e_bend_drv*(-drix-drkx)
        // TL: -drix
        particle_list[j].fxth += e_bend_drv *(-drix - drkx);
        particle_list[j].fyth += e_bend_drv *(-driy - drky);
        particle_list[j].fzth += e_bend_drv *(-driz - drkz);

        particle_list[k].fxth += e_bend_drv * drkx;
        particle_list[k].fyth += e_bend_drv * drky;
        particle_list[k].fzth += e_bend_drv * drkz;
    }

    return 0;
}

double pf::Force::ftorsion(vector<Particle> &particle_list,
        double &e_tors_tot) {
    for (int i = 0; i < param.npartM - 3; i++) {
        int j = i + 1;
        int k = i + 2;
        int l = i + 3;

        double xij = particle_list[i].x - particle_list[j].x;
        double yij = particle_list[i].y - particle_list[j].y;
        double zij = particle_list[i].z - particle_list[j].z;

        double xkj = particle_list[k].x - particle_list[j].x;
        double ykj = particle_list[k].y - particle_list[j].y;
        double zkj = particle_list[k].z - particle_list[j].z;

        double xkl = particle_list[k].x - particle_list[l].x;
        double ykl = particle_list[k].y - particle_list[l].y;
        double zkl = particle_list[k].z - particle_list[l].z;

        // rij is never used
        //double rij = sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2));
        double rkj = sqrt(pow(xkj, 2) + pow(ykj, 2) + pow(zkj, 2));
        double rkl = sqrt(pow(xkl, 2) + pow(ykl, 2) + pow(zkl, 2));

        double xmj = yij * zkj - ykj * zij;
        double ymj = zij * xkj - zkj * xij;
        double zmj = xij * ykj - xkj * yij;

        double xnk = ykj * zkl - ykl * zkj;
        double ynk = zkj * xkl - zkl * xkj;
        double znk = xkj * ykl - xkl * ykj;

        double xil = ymj * znk - ynk * zmj;
        double yil = zmj * xnk - znk * xmj;
        double zil = xmj * ynk - xnk * ymj;

        double rnk = sqrt(pow(xnk, 2) + pow(ynk, 2) + pow(znk, 2));
        double rmj = sqrt(pow(xmj, 2) + pow(ymj, 2) + pow(zmj, 2));

        if (pow(rnk, 2) < eps || pow(rmj, 2) < eps) break;

        double phi = (xnk * xmj + ynk * ymj + znk * zmj) / (rnk * rmj);
        if (fabs(phi) > eps1) {
            particle_list[l].x += eps2;
            particle_list[l].y += eps2;
            particle_list[l].z += eps2;

            xkl = particle_list[k].x - particle_list[l].x;
            ykl = particle_list[k].y - particle_list[l].y;
            zkl = particle_list[k].z - particle_list[l].z;

            rkl = sqrt(pow(xkl, 2) + pow(ykl, 2) + pow(zkl, 2));

            xnk = ykj * zkl - ykl * zkj;
            ynk = zkj * xkl - zkl * xkj;
            znk = xkj * ykl - xkl * ykj;

            xil = ymj * znk - ynk * zmj;
            yil = zmj * xnk - znk * xmj;
            zil = xmj * ynk - xnk * ymj;

            rnk = sqrt(pow(xnk, 2) + pow(ynk, 2) + pow(znk, 2));
            rmj = sqrt(pow(xmj, 2) + pow(ymj, 2) + pow(zmj, 2));

            if (pow(rnk, 2) < eps || pow(rmj, 2) < eps) continue;
            phi = (xnk * xmj + ynk * ymj + znk * zmj) / (rnk * rmj);
        }

        if(phi > eps1) phi = eps1;
        if(phi < -1.0 * eps1) phi = eps1;

        phi = acos(phi);
        double tmpvalue = xkj * xil + ykj * yil + zkj * zil; 
        phi = tmpvalue > 0 ? phi : (-1 * phi);
        double phi_0 = param.dihedral_nat[i];


        double e_tors = param.ck_phi1 * (1.0 - cos(phi - phi_0)) +
            param.ck_phi3 * (1.0 - cos(3.0 * (phi - phi_0)));

        double drv1 = param.ck_phi1 * sin(phi - phi_0);
        double drv2 = 3.0 * param.ck_phi3 * (4.0 * pow(cos(phi), 2.0) *
                sin(phi - 3.0 * phi_0) + 3.0 * cos(phi) * sin(3.0 * phi_0) -
                sin(phi) * cos(3.0 * phi_0));

        double e_tors_drv = -(drv1 + drv2);

        double tem = 0.0;
        if(i <= param.npart1 - 3) {
            tem = param.Alpha1;
        }
        else if(i > param.npart1) {
            tem = param.Alpha2;
        }
        e_tors_tot += tem * e_tors;

        double rijrkj = (xij * xkj + yij * ykj + zij * zkj) / pow(rkj, 2);
        double rklrkj = (xkl * xkj + ykl * ykj + zkl * zkj) / pow(rkj, 2);

        double drix = xmj * rkj / pow(rmj, 2);
        double driy = ymj * rkj / pow(rmj, 2);
        double driz = zmj * rkj / pow(rmj, 2);

        double drlx = -xnk * rkj / pow(rnk, 2);
        double drly = -ynk * rkj / pow(rnk, 2);
        double drlz = -znk * rkj / pow(rnk, 2);

        double drjx = (rijrkj - 1) * drix - rklrkj * drlx;
        double drjy = (rijrkj - 1) * driy - rklrkj * drly;
        double drjz = (rijrkj - 1) * driz - rklrkj * drlz;

        double drkx = (rklrkj - 1) * drlx - rijrkj * drix;
        double drky = (rklrkj - 1) * drly - rijrkj * driy;
        double drkz = (rklrkj - 1) * drlz - rijrkj * driz;

        particle_list[i].fxph += tem * e_tors_drv * drix;
        particle_list[i].fyph += tem * e_tors_drv * driy;
        particle_list[i].fzph += tem * e_tors_drv * driz;

        particle_list[j].fxph  += tem * e_tors_drv * drjx;
        particle_list[j].fyph  += tem * e_tors_drv * drjy;
        particle_list[j].fzph  += tem * e_tors_drv * drjz;

        particle_list[k].fxph += tem * e_tors_drv * drkx;
        particle_list[k].fyph += tem * e_tors_drv * drky;
        particle_list[k].fzph += tem * e_tors_drv * drkz;

        particle_list[l].fxph += tem * e_tors_drv * drlx;
        particle_list[l].fyph += tem * e_tors_drv * drly;
        particle_list[l].fzph += tem * e_tors_drv * drlz;
    }

    return 0;
}

double pf::Force::funbond(vector<Particle> &particle_list, double &e_unbond_tot,
        double &e_bind_tot) {

    if(param.Is_Solvation == 0) {
        funbond_without(particle_list, e_unbond_tot, e_bind_tot);
    }
    else if(param.Is_Solvation == 1) {
        funbond_with(particle_list, e_unbond_tot, e_bind_tot);
    }
    return 0;
}

/**
 * Solve unbond force with solvent
 */
double pf::Force::funbond_with(vector<Particle> &particle_list,
        double &e_unbond_tot, double &e_bind_tot) {
    param.gQ_f     = 0.0;
    param.gQ_f1    = 0.0;
    param.gQ_f2    = 0.0;
    param.gQ_b     = 0.0;
    param.gQ_w     = 0.0;
    double sum_rij = 0.0;
    param.gQ_non_f = 0.0;
    param.gQ_non_b = 0.0;

    // !     r''-r'
    param.dr_sol = 3.0;
    // !     r*-r'
    param.ddr_sol = 1.5;

    double xij = 0.0;
    double yij = 0.0;
    double zij = 0.0;
    double rij = 0.0;
    double tem = 0.0;
    // !  unbonded force for A:
    for (int k = 0; k < param.nunbond; k++) {
        int i = param.iun[k] - 1;
        int j = param.jun[k] - 1;

        double e_unbond = 0;
        double e_unbond_drv = 0;

        xij = (particle_list[i].x - particle_list[j].x);
        yij = (particle_list[i].y - particle_list[j].y);
        zij = (particle_list[i].z - particle_list[j].z);

        rij = sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2));

        if ((param.kunbond[k] == 2 && rij >= (2 * param.CritR_non + 3.0))
                || (param.kunbond[k] == 0 && rij >= (1.8 * param.sigma_ij))
                || (i > param.npartM && j > param.npartM)) continue;

        if (param.kunbond[k] == 1) {
            if (rij / param.kunbond[k] - param.gr0 < 8 * param.gdr) {
                tem = 1 / (1 + exp((rij / (param.runbond_nat[k] +
                                    param.ddr_sol) - param.gr0) / param.gdr));
                if (i <= param.npart1 && j <= param.npart1) {
                    param.gQ_f1 += tem;
                } else if (i <= param.npart1 && j > param.npart1) {
                    param.gQ_b += tem;
                } else if (param.iFlagMov == 2){
                    param.gQ_f2 += tem;
                }
            }
            if (i <= param.npart1 && j > param.npart1) {
                sum_rij += (rij - param.runbond_nat[k]);
            }

            if (rij < param.runbond_nat[k]) {
                double Zr = pow((param.runbond_nat[k] / rij), param.k_sol);
                double dZr = -param.k_sol * Zr / rij;
                e_unbond = param.epsil1 * Zr * (Zr - 2.0);
                e_unbond_drv = param.epsil1 * 2 * (Zr - 1) * dZr;
            } else if (rij < param.runbond_nat[k] + param.ddr_sol) {
                double Yr = pow((rij - (param.runbond_nat[k] + param.ddr_sol)),
                        2);
                double Yrn = pow(Yr, param.n_sol);
                double dYr = 2 * (rij - (param.runbond_nat[k] + param.ddr_sol));
                double dr2n = pow(param.ddr_sol, (2 * param.n_sol));
                // !  Note that C should be .../(r*-r')^4n
                double CC = 4 * param.n_sol * (param.epsil1
                        + param.epsilon_pp) / pow(dr2n, 2);
                e_unbond = CC * Yrn * (0.5 * Yrn - dr2n) / (2 * param.n_sol) +
                    param.epsilon_pp;
                e_unbond_drv = CC / (2 * param.n_sol) * (Yrn - dr2n) *
                    (param.n_sol * Yrn / Yr) * dYr;
            } else {
                double Yr = pow((rij - (param.runbond_nat[k] + param.ddr_sol)),
                        2);
                double dYr = 2 * (rij - (param.runbond_nat[k] + param.ddr_sol));
                double ddr = param.dr_sol - param.ddr_sol;
                double BB = param.epsilon_p * param.m_sol * pow(ddr,
                        (2 * (param.m_sol - 1)));
                double h1 = (1 - 1.0 / param.m_sol) * pow(ddr, 2) /
                    (param.epsilon_p / param.epsilon_pp + 1.0);
                double h2 = (param.m_sol - 1) * pow(ddr, (2 * param.m_sol)) /
                    (1.0 + param.epsilon_pp / param.epsilon_p);
                e_unbond = -BB * (Yr - h1) / (pow(Yr, param.m_sol) + h2);
                e_unbond_drv = -BB * dYr / (pow(Yr, param.m_sol) + h2) +
                    BB * (Yr - h1) / pow((pow(Yr, param.m_sol) + h2), 2) *
                    param.m_sol * pow(Yr, (param.m_sol - 1)) * dYr;
            }
            e_unbond_drv = -e_unbond_drv;

            if(i <= param.npart1 && j <= param.npart1) {
                e_unbond = param.Alpha1 * e_unbond;
                e_unbond_drv = param.Alpha1 * e_unbond_drv;
            }
            else if(i <= param.npart1 && j > param.npart1) {
                e_unbond = param.Beta * e_unbond;
                e_unbond_drv = param.Beta * e_unbond_drv;
            }
            else if(param.iFlagMov == 2) {
                e_unbond = param.Alpha2 * e_unbond;
                e_unbond_drv = param.Alpha2 * e_unbond_drv;
            }
        } else if (param.kunbond[k] == 2) {
            if((rij / param.CritR_non - param.gr0) < 8 * param.gdr) {
                tem = 1 / (1 +
                        exp((rij / param.CritR_non - param.gr0) / param.gdr));
                if(i <= param.npart1 && j > param.npart1) {
                    param.gQ_non_b = param.gQ_non_b + tem;
                } else if (param.iFlagMov == 2 || j <= param.npart1) {
                    param.gQ_non_f = param.gQ_non_f + tem;
                }
            }
            double rij12 = pow((param.sigma_ij / rij), 12.0);
            double rHP = exp(-1 * pow((rij - param.CritR_non), 2) / 2);
            double Atem = 0.0;
            if(i <= param.npart1 && j <= param.npart1) {
                Atem = param.Alpha1;
            }
            else if(i <= param.npart1 && j > param.npart1) {
                // ! Atem = sqrt(Alpha1*Alpha2)
                Atem = param.Beta;
            }
            else {
                Atem = param.Alpha2;
            }
            e_unbond = Atem * param.epsil2 * rij12 - param.Delta *
                param.epsil1 * rHP;
            e_unbond_drv = Atem * 12 * param.epsil2 * rij12 / rij -
                param.Delta * param.epsil1 * rHP * (rij - param.CritR_non);
        } else {
            double rij12 = pow((param.sigma_ij / rij), 12.0);
            double Atem = 0.0;
            if(i <= param.npart1 && j <= param.npart1) {
                Atem = param.Alpha1;
            }
            else if (i <= param.npart1 && j > param.npart1) {
                // ! Atem = sqrt(Alpha1*Alpha2)
                Atem = param.Beta;
            }
            else {
                Atem = param.Alpha2;
            }
            e_unbond = Atem * param.epsil2 * rij12;
            e_unbond_drv = Atem * 12 * param.epsil2 * rij12 / rij;
        }

        e_unbond_tot += e_unbond;
        if(i <= param.npart1 && j > param.npart1) e_bind_tot += e_unbond;

        if(i <= param.npartM) {
            particle_list[i].fxun += e_unbond_drv * xij / rij;
            particle_list[i].fyun += e_unbond_drv * yij / rij;
            particle_list[i].fzun += e_unbond_drv * zij / rij;
        }
        if(j <= param.npartM) {
          particle_list[j].fxun -= e_unbond_drv * xij / rij;
          particle_list[j].fyun -= e_unbond_drv * yij / rij;
          particle_list[j].fzun -= e_unbond_drv * zij / rij;
        }
    } // fortran code: 40 continue

    if (param.iFlagMov != 2) param.gQ_f2 = 0;
    param.gQ_f = param.gQ_f1 + param.gQ_f2;

    // !  Lagrange constrant potential:
    // !  Fixed B
    double eGr_f1 = param.ga1_f1 * (param.gQ_f1 - param.gQ0_f1) +
        param.ga2_f1 * pow(param.gQ_f1 - param.gQ0_f1, 2);
    double eGr_f2 = param.ga1_f2 * (param.gQ_f2 - param.gQ0_f2) +
        param.ga2_f2 * pow(param.gQ_f2 - param.gQ0_f2, 2);
    double eGr_b = param.ga1_b * (param.gQ_b - param.gQ0_b ) +
        param.ga2_b * pow(param.gQ_b - param.gQ0_b, 2);
    double eGr_f = eGr_f1 + eGr_f2;
    double eGr_w = param.ga1_w * (param.gQ_w - param.gQ0_w ) +
        param.ga2_w * pow(param.gQ_w - param.gQ0_w, 2);
    double eGr = eGr_f + eGr_b + eGr_w;

    if (eGr < 1e-5 && eGr > -1e-5) return 0;

    for (int k = 0; k < param.nunbond; k++) {
        int i = param.iun[k] - 1;
        int j = param.jun[k] - 1;

        double e_unbond_drv = 0;

        if (param.kunbond[k] == 1) {
        xij = particle_list[i].x - particle_list[j].x;
        yij = particle_list[i].y - particle_list[j].y;
        zij = particle_list[i].z - particle_list[j].z;
        rij = sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2));

        if( i <= param.npart1 && j <= param.npart1) {
            e_unbond_drv = param.ga1_f1 + 2 * param.ga2_f1 * (param.gQ_f1 -
                    param.gQ0_f1);
        }
        else if(i <= param.npart1 && j > param.npart1) {
          e_unbond_drv = param.ga1_b + 2 * param.ga2_b *
              (param.gQ_b - param.gQ0_b ) + param.ga1_w +
              2 * param.ga2_w * (param.gQ_w - param.gQ0_w);
        }
        else {
            e_unbond_drv = param.ga1_f2 + 2 * param.ga2_f2 * (param.gQ_f2 -
                    param.gQ0_f2);
        }

        if ((rij / (param.runbond_nat[k] + param.ddr_sol) - param.gr0) <
                -8 * param.gdr || (rij / (param.runbond_nat[k] + param.ddr_sol)
                    - param.gr0) > 8 * param.gdr) {
            e_unbond_drv = 0.0;
        } else {
            tem = exp((rij / (param.runbond_nat[k] + param.ddr_sol) - param.gr0)
                    / param.gdr);
            e_unbond_drv = e_unbond_drv * tem / pow(1 + tem, 2) /
                ((param.runbond_nat[k] + param.ddr_sol) * param.gdr);
        }

        if (i <= param.npartM) {
            particle_list[i].fxun += e_unbond_drv * xij / rij;
            particle_list[i].fyun += e_unbond_drv * yij / rij;
            particle_list[i].fzun += e_unbond_drv * zij / rij;
        }
        if (j <= param.npartM) {
            particle_list[j].fxun -= e_unbond_drv * xij / rij;
            particle_list[j].fyun -= e_unbond_drv * yij / rij;
            particle_list[j].fzun -= e_unbond_drv * zij / rij;
        }
    }
}
    return 0;
}

/**
 * Solve unbond force without solvent
 */
double pf::Force::funbond_without(vector<Particle> &particle_list,
        double &e_unbond_tot, double &e_bind_tot) {
    param.gQ_f     = 0.0;
    param.gQ_f1    = 0.0;
    param.gQ_f2    = 0.0;
    param.gQ_b     = 0.0;
    param.gQ_w     = 0.0;
    double sum_rij = 0.0;
    param.gQ_non_f = 0.0;
    param.gQ_non_b = 0.0;

    // !     r''-r'
    param.dr_sol = 3.0;
    // !     r*-r'
    param.ddr_sol = 1.5;

    // Some temporary variables
    double xij = 0.0;
    double yij = 0.0;
    double zij = 0.0;
    double rij = 0.0;
    double tem = 0.0;
    double e_unbond = 0;
    double e_unbond_drv = 0;
    // !  unbonded force for A:
    for (int k = 0; k < param.nunbond; k++) {
        int i = param.iun[k] - 1;
        int j = param.jun[k] - 1;

        xij = (particle_list[i].x - particle_list[j].x);
        yij = (particle_list[i].y - particle_list[j].y);
        zij = (particle_list[i].z - particle_list[j].z);

        // ! cancel for Q_w
        // !      if ( kunbond(k) == 1 .and. rij >= (2*runbond_nat(k)+3.0) .or. &
        if ((param.kunbond[k] == 2 && rij >= (2 * param.CritR_non + 3.0))
                || (param.kunbond[k] == 0 && rij >= (1.8 * param.sigma_ij))
                || (i > param.npartM && j > param.npartM)) continue;

        if (param.kunbond[k] == 1) {
            if((rij / param.runbond_nat[k] - param.gr0) < 8 * param.gdr) {
                tem = 1 / (1 + exp((rij / param.runbond_nat[k] - param.gr0) /
                            param.gdr));
                if (i <= param.npart1 && j <= param.npart1) {
                    param.gQ_f1 += tem;
                } else if (i <= param.npart1 && j > param.npart1) {
                    param.gQ_b += tem;
                } else if (param.iFlagMov == 2) {
                    param.gQ_f2 += tem;
                }
            } 
            if(i <= param.npart1 && j > param.npart1) {
                sum_rij += (rij - param.runbond_nat[k]);
            }
            double rij10 = pow((param.runbond_nat[k] / rij), 10.0);
            double rij12 = pow((param.runbond_nat[k] / rij), 12.00);
            if (i <= param.npart1 && j <= param.npart1) {
                e_unbond = param.Alpha1 * param.epsil1 * (5 * rij12 -
                        6 * rij10);
                e_unbond_drv = param.Alpha1 * 60 * param.epsil1 * (rij12 -
                        rij10) / rij;
            }
            else if (i <= param.npart1 && j > param.npart1) {
                e_unbond = param.Beta * param.epsil1 * (5 * rij12 - 6 * rij10);
                e_unbond_drv = param.Beta * 60 * param.epsil1 * (rij12 -
                        rij10) / rij;
            }
            else if (param.iFlagMov == 2) {
                e_unbond = param.Alpha2 * param.epsil1 * ( 5 * rij12 -
                        6 * rij10 );
                e_unbond_drv = param.Alpha2 * 60 * param.epsil1 * ( rij12 -
                        rij10 ) / rij;
            }
        } else if (param.kunbond[k] == 2) {
            if ((rij / param.CritR_non - param.gr0) < 8 * param.gdr) {
                tem = 1 / (1 + exp((rij / param.CritR_non - param.gr0) /
                            param.gdr));
                if (i <= param.npart1 && j > param.npart1) {
                    param.gQ_non_b += tem;
                }
                else if (param.iFlagMov == 2 || j <= param.npart1) {
                    param.gQ_non_f += tem;
                }
            }

            double rij12 = pow(param.sigma_ij / rij, 12.0);
            double rHP  = exp(-1 * pow(rij - param.CritR_non, 2) / 2);
            double Atem = 0.0;
            if (i <= param.npart1 && j <= param.npart1) {
                Atem = param.Alpha1;
            } else if (i <= param.npart1 && j > param.npart1) {
                // ! Atem = sqrt(Alpha1*Alpha2)
                Atem = param.Beta;
            } else {
               Atem = param.Alpha2;
            }
            e_unbond = Atem * param.epsil2 * rij12 - param.Delta *
                param.epsil1 * rHP;
            e_unbond_drv = Atem * 12 * param.epsil2 * rij12 / rij -
                param.Delta * param.epsil1 * rHP * (rij - param.CritR_non);
        } else {
            double rij12 = pow(param.sigma_ij / rij, 12.0);
            double Atem = 0.0;
            if (i <= param.npart1 && j <= param.npart1) {
                Atem = param.Alpha1;
            } else if (i <= param.npart1 && j > param.npart1) {
                // ! Atem = sqrt(Alpha1*Alpha2)
                Atem = param.Beta;
            } else {
                Atem = param.Alpha2;
            }
            e_unbond = Atem * param.epsil2 * rij12;
            e_unbond_drv = Atem * 12 * param.epsil2 * rij12 / rij;
        }

        e_unbond_tot += e_unbond;
        if (i <= param.npart1 && j > param.npart1) e_bind_tot += e_unbond;

        if (i <= param.npartM) {
            particle_list[i].fxun += e_unbond_drv * xij / rij;
            particle_list[i].fyun += e_unbond_drv * yij / rij;
            particle_list[i].fzun += e_unbond_drv * zij / rij;
        }
        if (j <= param.npartM) {
            particle_list[j].fxun -= e_unbond_drv * xij / rij;
            particle_list[j].fyun -= e_unbond_drv * yij / rij;
            particle_list[j].fzun -= e_unbond_drv * zij / rij;
        }
    } // fortran code: 40 continue

    if (param.iFlagMov != 2) param.gQ_f2 = 0;
    param.gQ_f = param.gQ_f1 + param.gQ_f2;

    // !  Lagrange constrant potential:
    // !  Fixed B
    param.gQ_w = param.gQ_b - param.alpha_Qw * sum_rij;
    double eGr_f1 = param.ga1_f1 * (param.gQ_f1 - param.gQ0_f1) +
        param.ga2_f1 * pow(param.gQ_f1 - param.gQ0_f1, 2);
    double eGr_f2 = param.ga1_f2 * (param.gQ_f2 - param.gQ0_f2 ) +
        param.ga2_f2 * pow(param.gQ_f2 - param.gQ0_f2, 2);
    double eGr_b = param.ga1_b * (param.gQ_b - param.gQ0_b ) +
        param.ga2_b * pow(param.gQ_b - param.gQ0_b, 2);
    double eGr_f = eGr_f1 + eGr_f2;
    double eGr_w = param.ga1_w * (param.gQ_w - param.gQ0_w ) +
        param.ga2_w * pow(param.gQ_w - param.gQ0_w, 2);
    param.eGr = eGr_f + eGr_b + eGr_w;

    if(param.eGr < 1e-5 && param.eGr > -1e-5) return 0;

    for (int k = 0; k < param.nunbond; k++) {
        int i = param.iun[k] - 1;
        int j = param.jun[k] - 1;

        e_unbond_drv = 0;
        if (param.kunbond[k] == 1) {
            xij = particle_list[i].x - particle_list[j].x;
            yij = particle_list[i].y - particle_list[j].y;
            zij = particle_list[i].z - particle_list[j].z;
            rij = sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2));

            if (i <= param.npart1 && j <= param.npart1) {
                e_unbond_drv = param.ga1_f1 + 2 * param.ga2_f1 * (param.gQ_f1 -
                        param.gQ0_f1);
            }
            else if(i <= param.npart1 && j > param.npart1) {
                e_unbond_drv = param.ga1_b + 2 * param.ga2_b *
                    (param.gQ_b - param.gQ0_b ) + param.ga1_w +
                    2 * param.ga2_w * (param.gQ_w - param.gQ0_w);
            } else {
                e_unbond_drv = param.ga1_f2 + 2 * param.ga2_f2 * (param.gQ_f2 -
                        param.gQ0_f2);
            }

            if ((rij / (param.runbond_nat[k] + param.ddr_sol) - param.gr0) <
                    -8 * param.gdr || (rij / (param.runbond_nat[k] +
                            param.ddr_sol) - param.gr0) > 8 * param.gdr) {
                e_unbond_drv = 0.0;
            } else {
                tem = exp((rij / (param.runbond_nat[k] + param.ddr_sol) -
                            param.gr0) / param.gdr);
                e_unbond_drv = e_unbond_drv * tem / pow(1 + tem, 2) /
                    ((param.runbond_nat[k] + param.ddr_sol) * param.gdr);
            }

            if (i <= param.npart1 && j > param.npart1) {
                e_unbond_drv = e_unbond_drv + param.alpha_Qw * (param.ga1_w +
                        2 * param.ga2_w * (param.gQ_w - param.gQ0_w));
            }

            if (i <= param.npartM) {
                particle_list[i].fxun += e_unbond_drv * xij / rij;
                particle_list[i].fyun += e_unbond_drv * yij / rij;
                particle_list[i].fzun += e_unbond_drv * zij / rij;
            }
            if (j <= param.npartM) {
                particle_list[j].fxun -= e_unbond_drv * xij / rij;
                particle_list[j].fyun -= e_unbond_drv * yij / rij;
                particle_list[j].fzun -= e_unbond_drv * zij / rij;
            }
        }
    }

    return 0;
}
/////////////////////////////////////////////////////////////////////////
/**********************End of definition of class Force*****************/
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
/*******************Definition of class Simulation**********************/
/////////////////////////////////////////////////////////////////////////
pf::Simulation::Simulation(string conf_filename) : param(conf_filename), force(param){
}

int pf::Simulation::intpar(double &enerkin) {
    param.temp = param.temp * param.epsil / param.boltz;
    double timeunit = sqrt(param.amass / param.epsil) * param.sigma_ij;
    // !       dt=0.005*timeunit
    // !       gm=0.05/timeunit
    param.dt = param.s_dt * timeunit;
    param.gm = param.s_gm / timeunit;

    // c_*: intensity arguments
    param.c_0 = param.dt * param.gm / (2.0 * param.amass);
    param.c_1 = (1.0 - param.c_0) * (1.0 - param.c_0 + pow(param.c_0, 2));
    param.c_2 = (param.dt / (2 * param.amass)) *
      (1.0 - param.c_0 + pow(param.c_0, 2));
    param.randconst =
      sqrt(2.0 * param.boltz * param.temp * param.gm / param.dt);

    // !  no change ck_r    -ZRLiu
    // !       ck_r   = enscale*ck_r
    if (param.iFixKr == 0) param.ck_r = param.enscale * param.ck_r;
    param.ck_tht = param.enscale * param.ck_tht;
    param.ck_phi1= param.enscale * param.ck_phi1;
    param.ck_phi3= param.enscale * param.ck_phi3;
    param.epsil1 = param.enscale * param.epsil1;
    param.epsil2 = param.enscale * param.epsil2;
    param.epsilon_p= param.enscale * param.epsilon_p;
    param.epsilon_pp = param.enscale * param.epsilon_pp;

    // ! ccccccccccccccccccccccccccccccccc
    // !     initiliaze velocities
    // ! ccccccccccccccccccccccccccccccccc
    double sumvx = 0.0;
    double sumvy = 0.0;
    double sumvz = 0.0;

    double vm = sqrt(param.temp * param.boltz / param.amass);
    for (int i = 0; i < param.npartM; i++) {
      particle_list[i].vx = vm * gauss(xsi);
      particle_list[i].vy = vm * gauss(xsi);
      particle_list[i].vz = vm * gauss(xsi);
      sumvx = sumvx + particle_list[i].vx;
      sumvy = sumvy + particle_list[i].vy;
      sumvz = sumvz + particle_list[i].vz;
    }

    sumvx /= param.npartM;
    sumvy /= param.npartM;
    sumvz /= param.npartM;
    enerkin = 0.0;

    for (int i = 0; i < param.npartM; i++) {
        particle_list[i].vx = particle_list[i].vx - sumvx;
        particle_list[i].vy = particle_list[i].vy - sumvy;
        particle_list[i].vz = particle_list[i].vz - sumvz;
        enerkin += (pow(particle_list[i].vx, 2) + pow(particle_list[i].vy, 2) +
                pow(particle_list[i].vz, 2));
    }

    enerkin =  0.5 * enerkin * param.amass;
    // TL: It seems that tempins is unused.
    // double tempins = 2.0 * enerkin / (3.0 * param.npartM * param.boltz);

    // !  initialize histogram data
    for (int i = 0; i < param.nbin_f; i++) {
        statis.PFbin[i] = 0.0;
        for (int j = 0; j < param.nbin_b; j++) {
            statis.PBbin[j] = 0.0;
            statis.PFBbin[i][j] = 0.0;
            for (int k = 0; k < param.nEbin; k++) {
                statis.PFBEbin[i][j][k] = 0.0;
            } // fortran code: 50 continue
            for (int k = 0; k < param.nEbbin; k++) {
                statis.PQbEbbin[j][k] = 0.0;
            } // fortran code: 55 continue
            for (int l = 0; l < param.nRbin; l++) {
                statis.PFBRbin[i][j][l] = 0.0;
            } // fortran code: 60 continue
            for (int k = 0; k < param.nbin_f; k++) {
                statis.PFFBbin[i][k][j] = 0.0;
            } // fortran code: 70 continue
        } // fortran code: 40 continue
    } // fortran code: 30 continue

    // if ( iFlagMov.le.0 .or. IsWbin.eq.0 ) goto 180
    if (param.iFlagMov <= 0 || param.IsWbin == 0 ) return 0;
    for (int i = 0; i < param.nWbin; i++) {
        statis.PQwbin[i] = 0;
        for (int j = 0; j < param.nbin_b; j++) {
            statis.PQwQbbin[i][j] = 0;
        } // fortran code: 110 continue
        for (int l = 0; l < param.nRbin; l++) {
            // PQwRbin(i,l,1)=0 PQwRbin(i,l,2)=0
            statis.PQwRbin[i][l][0] = 0;
            statis.PQwRbin[i][l][1] = 0;
        } // fortran code: 120 continue
        for (int k = 0; k < param.nEbbin; k++) {
            // PQwEbbin(i,k,1)=0 PQwEbbin(i,k,2)=0
            statis.PQwEbbin[i][k][0] = 0;
            statis.PQwEbbin[i][k][1] = 0;
        } // fortran code: 130 continue
        for (int j = 0; j < param.nbin_f; j++) {
            statis.PQwQfbin[i][j] = 0;
        } // fortran code: 140 continue
    } // fortran code: 150 continue
    // fortran code: 180 continue

    return 0;
}

/**
 * Generate the filename of intermediate results
 */
int pf::Simulation::formate_output_filenames() {
    output_filenames.resize(35);
    // !  file to save output information: (Output*)
    string Output_filename = "Output" + param.pdb_filename;
    output_filenames[9] = Output_filename;
    // open ( 9, file = filename2(1:LF2), status = 'replace' )
    // !  file to save Histogram P(Q_b): (QbHist*)
    string QbHist_filename = "QbHist" + param.pdb_filename;
    output_filenames[14] = QbHist_filename;
    // if(iFlagMov.gt.0) open (14,file=filename2(1:LF2), status = 'replace' ) 
    // !  file to save Histogram P(Q_f): (QfHist*)
    string QfHist_filename = "QfHist" + param.pdb_filename;
    output_filenames[15] = QfHist_filename;
    // open (15,file=filename2(1:LF2), status = 'replace' )
    // !  file to save Histogram P(Q_f,Q_b): (QfbHist*)
    string QfbHist_filename = "QfbHist" + param.pdb_filename;
    output_filenames[19] = QfbHist_filename;
    // if(iFlagMov.gt.0) open (19,file=filename2(1:LF2), status = 'replace' )
    // !  file to save the (Q,E)-histogram: (QEHist*.21)
    string QEHist_filename = "QEHist" + param.pdb_filename;
    output_filenames[23] = QEHist_filename;
    // if(IsEbin.eq.1)  open (23,file=filename2(1:LF2), status = 'replace' )
    // !  file to save the (Q,R)-histogram: (QRHist*.21)
    string QRHist_filename = "QRHist" + param.pdb_filename;
    output_filenames[24] = QRHist_filename;
    // if(IsRbin.eq.1)  open (24,file=filename2(1:LF2), status = 'replace' )
    // !  file to save the (Qf1,Qf2,Qb)-histogram: (QffbHist*.21)
    string QffbHist_filename = "QffbHist" + param.pdb_filename;
    output_filenames[25] = QffbHist_filename;
    // if(iFlagMov.eq.2)  open (25,file=filename2(1:LF2), status = 'replace' )
    // !  file to save the output conformation, d*.mol
    string d_mol_filename = "d" + param.pdb_filename + ".mol";
    output_filenames[22] = d_mol_filename;
    // open (22,file=filename2(1:LF2), status = 'replace' )
    // !  file to save the output conformation coordinate for further
    // simulation, d*.X 
    string d_mol_dat_filename = "d" + param.pdb_filename + ".mol.dat";
    output_filenames[27] = d_mol_dat_filename;
    // open(27,file=filename2(1:LF2), status = 'replace' )
    // !  file to save the output conformation coordinate for further
    // simulation, d*.X 
    string d_X_filename = "d" + param.pdb_filename + ".X";
    output_filenames[26] = d_X_filename;
    // open (26,file=filename2(1:LF2), status = 'replace' )
    // !  file to save the binding/unbinding time and transmission coefficient,
    // d*.tc
    string d_tc_filename = "d" + param.pdb_filename + ".tc";
    output_filenames[28] = d_tc_filename;
    // open (28,file=filename2(1:LF2), status = 'replace' )
    // !  file to save Histogram P(Q_b,E_b): (QbEbHist*)
    string QbEbHist_filename = "QbEbHist" + param.pdb_filename;
    output_filenames[29] = QbEbHist_filename;
    // if(iFlagMov.gt.0 .and. IsEbbin.eq.1) open (29,file=filename2(1:LF2),
    // status = 'replace' )
    // !  file to save Histogram P(Q_w): (QwHist*)
    string QwHist_filename = "QwHist" + param.pdb_filename;
    output_filenames[30] = QwHist_filename;
    // if(iFlagMov.gt.0 .and. IsWbin.eq.1) open (30,file=filename2(1:LF2),
    // status = 'replace' )
    // !  file to save Histogram P(Q_w,Q_b): (QwQbHist*)
    string QwQbHist_filename = "QwQbHist" + param.pdb_filename;
    output_filenames[31] = QwQbHist_filename;
    // if(iFlagMov.gt.0 .and. IsWbin.eq.1) open (31,file=filename2(1:LF2),
    // status = 'replace' )
    // !  file to save Histogram P(Q_w,R): (QwRHist*)
    string QwRHist_filename = "QwRHist" + param.pdb_filename;
    output_filenames[32] = QwRHist_filename;
    // if(iFlagMov.gt.0 .and. IsWbin.eq.1) open (32,file=filename2(1:LF2),
    // status = 'replace' )
    // !  file to save Histogram P(Q_w,E_b): (QwEbHist*)
    string QwEbHist_filename = "QwEbHist" + param.pdb_filename;
    output_filenames[33] = QwEbHist_filename;
    // if(iFlagMov.gt.0 .and. IsWbin.eq.1) open (33,file=filename2(1:LF2),
    // status = 'replace' )
    // !  file to save Histogram P(Q_w,Q_f1): (QwEbHist*)
    string QwQfHist_filename = "QwQfHist" + param.pdb_filename;
    output_filenames[34] = QwQfHist_filename;
    // if(iFlagMov.gt.0 .and. IsWbin.eq.1) open (34,file=filename2(1:LF2),
    // status = 'replace' )

    return 0;
}

/**
 * Output the simulation information
 */
int pf::Simulation::output_info() {
    cout << setw(50) << "Length of chain 1 (npart1):" << param.npart1 << endl;
    cout << setw(50) << "Number of moving chains (iFlagMov):" << param.iFlagMov
        << endl;
    cout << setw(50) << "Total length (nparttol):" << param.nparttol << endl;
    cout << setw(50) << "Chain-Residue number gap (nResGap):" << param.nResGap
        << endl;
    cout << setw(50) << "interaction strength epsilon:" << param.epsil << endl;
    cout << setw(50) << "interaction strength epsilon_1:" << param.epsil1 <<
        endl;
    cout << setw(50) << "interaction strength epsilon_2:" << param.epsil2 <<
        endl;
    cout << setw(50) << "interaction strength ck_r:" << param.ck_r << endl;
    cout << setw(50) << "interaction strength ck_tht:" << param.ck_tht << endl;
    cout << setw(50) << "interaction strength ck_phi1:" << param.ck_phi1 <<
        endl;
    cout << setw(50) << "interaction strength ck_phi3:" << param.ck_phi3 <<
        endl;
    cout << setw(50) << "interaction strength scaling param:" << param.enscale
        << endl;
    cout << setw(50) << "Solvation k:" << param.k_sol << endl;
    cout << setw(50) << "Solvation n:" << param.n_sol << endl;
    cout << setw(50) << "Solvation m:" << param.m_sol << endl;
    cout << setw(50) << "Solvation epsilon pri:" << param.epsilon_p << endl;
    cout << setw(50) << "Solvation epsilon pripri:" << param.epsilon_pp << endl;
    cout << setw(50) << "dt:" << param.s_dt << endl;
    cout << setw(50) << "gamma:" << param.s_gm << endl;
    cout << setw(50) << "Is Fix Kr:" << param.iFixKr << endl;

    cout << setw(50) << "temperature:" << param.temp << endl;
    cout << setw(50) << "cut-off:" << param.gamma << endl;
    cout << setw(50) << "hard core distance:" << param.sigma_ij << endl;
    cout << setw(50) << "total simulation step:" << param.nstep << endl;
    cout << setw(50) << "frequency of snapshots:" << param.nsnap << endl;
    cout << setw(50) << "Avogadro number:" << param.avsn0 << endl;
    cout << setw(50) << "Boltzmann constant:" << param.boltz << endl;
    cout << setw(50) << "random seed:" << param.rand_num << endl;

    cout << setw(50) << "Lagrange constraint a1_f1:" << param.ga1_f1 << endl;
    cout << setw(50) << "Lagrange constraint a2_f1:" << param .ga2_f1 << endl;
    cout << setw(50) << "Lagrange constraint Q0_f1:" << param.gQ0_f1 << endl;

    cout << setw(50) << "Lagrange constraint a1_f2:" << param.ga1_f2 << endl;
    cout << setw(50) << "Lagrange constraint a2_f2:" << param.ga2_f2 << endl;
    cout << setw(50) << "Lagrange constraint Q0_f2:" << param.gQ0_f2 << endl;

    cout << setw(50) << "Lagrange constraint a1_b:" << param.ga1_b << endl;
    cout << setw(50) << "Lagrange constraint a2_b:" << param.ga2_b << endl;
    cout << setw(50) << "Lagrange constraint Q0_b:" << param.gQ0_b << endl;

    // The following three lines is unused
    // cout << setw(50) << "Lagrange constraint a1_b:" << param.ga1_2 << endl;
    // cout << setw(50) << "Lagrange constraint a2_b:" << param.ga2_2 << endl;
    // cout << setw(50) << "Lagrange constraint Q0_b:" << param.gQ0_2 << endl;
    cout << setw(50) << "Parameter for Q_w as alpha_Qw:" << param.alpha_Qw << 
        endl;

    cout << setw(50) << "Lagrange constraint r0:" << param.gr0 << endl;
    cout << setw(50) << "Lagrange constraint dr:" << param.gdr << endl;

    cout << setw(50) << "nConformOutput:" << param.nConformOutput << endl;
    cout << setw(50) << "nOutput0 :" << param.nOutput0 << endl;
    cout << setw(50) << "ndOutput:" << param.ndOutput << endl;
    cout << setw(50) << "outQf1_i:" << param.outQf1_i << endl;
    cout << setw(50) << "outQf1_f:" << param.outQf1_f << endl;
    cout << setw(50) << "outQf2_i:" << param.outQf2_i << endl;
    cout << setw(50) << "outQf2_f:" << param.outQf2_f << endl;
    cout << setw(50) << "outQb_i:" << param.outQb_i << endl;
    cout << setw(50) << "outQb_f:" << param.outQb_f << endl;

    cout << setw(50) << "Histogram nbin_f:" << param.nbin_f << endl;
    cout << setw(50) << "Histogram dbin_f:" << param.dbin_f << endl;
    cout << setw(50) << "Histogram nbin_b:" << param.nbin_b << endl;
    cout << setw(50) << "Histogram dbin_b:" << param.dbin_b << endl;

    cout << setw(50) << "num. of conformation:" << param.nConform << endl;
    cout << setw(50) << "num. of conformation to skip:" << param.nCon0  << endl;
    cout << setw(50) << "nRunConf:" << param.nRunConf << endl;
    cout << setw(50) << "gQbnativ:" << param.gQbnativ << endl;
    cout << setw(50) << "gQbdenatural:" << param.gQbdenatural << endl;
    cout << setw(50) << "gQf1nativ:" << param.gQf1nativ << endl;
    cout << setw(50) << "gQf1denatural:" << param.gQf1denatural << endl;
    cout << setw(50) << "gQf2nativ:" << param.gQf2nativ << endl;
    cout << setw(50) << "gQf2denatural:" << param.gQf2denatural << endl;


   cout << setw(50) << "Force scheme:" << param.Is_Solvation << endl;
   cout << setw(50) << "If save E-Histogram:" << param.IsEbin << endl;
   cout << setw(50) << "E-Histogram nbin:" << param.nEbin << endl;
   cout << setw(50) << "E-Histogram dbin:" << param.dEbin << endl;
   cout << setw(50) << "E-Histogram vbin0:" << param.vEbin0 << endl;
   cout << setw(50) << "If save R-Histogram:" << param.IsRbin << endl;
   cout << setw(50) << "R-Histogram nbin:" << param.nRbin << endl;
   cout << setw(50) << "R-Histogram dbin:" << param.dRbin << endl;
   cout << setw(50) << "R-Histogram vbin0:" << param.vRbin0 << endl;
   cout << setw(50) << "If save Qw-Histogram:" << param.IsWbin << endl;
   cout << setw(50) << "Qw-Histogram nbin:" << param.nWbin << endl;
   cout << setw(50) << "Qw-Histogram dbin:" << param.dWbin << endl;
   cout << setw(50) << "Qw-Histogram vbin0:" << param.vWbin0 << endl;
   cout << setw(50) << "cri_Qb for Qw-Histogram:" << param.cri_Qb << endl;
   cout << setw(50) << "periodic boundary condiction size:" << param.pL << endl;
   cout << setw(50) << "periodic boundary condiction buffer:" << param.dl <<
       endl;
   cout << setw(50) << "Ratio of intrachain interaction 1:" << param.Alpha1 <<
       endl;
   cout << setw(50) << "Ratio of intrachain interaction 2:" << param.Alpha2 <<
       endl;
   cout << setw(50) << "Ratio of interchain interaction:" << param.Beta << endl;
   cout << setw(50) << "Delta:" << param.Delta << endl;
   cout << setw(50) << "CritR_non:" << param.CritR_non << endl;
    return 0;
}

/**
 * BUG: initial conformation file may store multiple inital conformation
 * BUG: 只读了第一个构象
 * Read particles from file
 */
int pf::Simulation::read_nativeconform(string filename) {
    fstream fin(filename, std::ifstream::in);
    if (!fin.is_open()) {
        cerr << "failed to open " << filename << '\n';
        return 1;
    }

    // Repeatedly read the initial coordinates of particle
    for (int i = 0; i < param.nparttol; i++) {
        double x, y, z;
        fin >> x >> y >> z;
        Particle particle(x, y, z);
        particle_list.push_back(particle);
    }
    return particle_list.size();
}
/** Read appNCS_****.dat
 */
int pf::Simulation::read_appNCS(string filename) {
    fstream fin(filename, std::ifstream::in);
    if (!fin.is_open()) {
        cerr << "failed to open " << filename << '\n';
        return 1;
    }

    // Repeatedly read the initial coordinates of particle
    string skip_others;
    // fortran code: nunbond=0. (TL: I guess: unbond is similar to the indexi.)
    param.nunbond = -1;
    for (int i = 0; i < MAXUNBO; i++) {
        int item1 = -1, item2 = -1, item3 = -1;
        fin >> item1 >> item2 >> item3;
        if( item1 == -1)
            break;
        if (item1 > param.npart1) {
            if (item1 < param.npart1 + param.nResGap) {
                cerr << "contact map error 4.\n";
            }
            item1 -= param.nResGap;
        }
        if (item2 > param.npart1) {
            if (item2 < param.npart1 + param.nResGap) {
                cerr << "contact map error 5.\n";
            }
            item2 -= param.nResGap;
        }
        if (item1 <= 0 || item1 > param.nparttol) {
            cerr << "contact map error 1.\n";
        }
        if (item2 <= 0 || item2 > param.nparttol) {
            cerr << "contact map error 2.\n";
        }
        if (item1 == item2) cerr << "contact map error 3.\n";
        // fortran code: nunbond=nunbond+1
        param.nunbond += 1;
        if (item1 < item2) {
            param.iun[param.nunbond] = item1;
            param.jun[param.nunbond] = item2;
        } else {
            param.iun[param.nunbond] = item2;
            param.jun[param.nunbond] = item1;
        }
        param.kunbond[param.nunbond] = item3;
        // BUG: I can't understand the logical relationship
        if ((item1 > param.npartM && item2 > param.npartM) || 
                (param.iFlagMov == 0 && 
                 (item1 > param.npart1 || item2 > param.npart1))) {
            param.nunbond -= 1;
        }

        // Skip the last string
        getline(fin, skip_others);
    }
    cout << setw(50) << "Effective Number of unbond:" << param.nunbond << endl;
    return param.nunbond;
}

/**
 * Read initial conformation of protein
 * 修改：只读一个构象
 */
int pf::Simulation::read_initalconform(string filename) {
    fstream fin(filename, std::ifstream::in);
    if (!fin.is_open()) {
        cerr << "failed to open " << filename << '\n';
        return 1;
    }

    // TL:Crucial code for native conformation equals initial conformation
    // Fortran code:
    //        do 316 j=1,nCon0
    //          do 315 i=1,nparttol
    //            read(20,*) xinit(i),yinit(i),zinit(i)
    //  315      continue
    //  316    continue
    for (int j = 0; j < param.nparttol; j++) {
        fin >> particle_list[j].xinit >> particle_list[j].yinit >>
            particle_list[j].zinit;
    }

    return 0;
}

/**
 * Start the simulation
 */
int pf::Simulation::start_simulation() {
    string msg = "";

    //=======================Initial work start=====================/
    // !-------------read parameters-------------------
    // Parameter Reading and parsing are in the constructor of class Simulation
    
    // !-------------open files to save-------------------
    // The work of intermediate result is complished by class Logger
    // Formate_output_filenames
    formate_output_filenames();
    output_info();

    // Read native conformation: unit_20<-->initialConf, unit_2<-->NativeConf
    // fortran code: read (2,*)(x(j),y(j),z(j),j=1,nparttol)
    read_nativeconform(param.nativeconform_filename);
    param.nunbond = read_appNCS(param.appNCS_filename);

    // Invoke the function (TODO: supply functional describes)
    setrn(param.ioctsd);
    nativeinformation();
    intpar(enerkin);

    // !  total number of folding/unfolding events
    // Meaningful when dynamics simulation, otherwise they are useless
    // (动力学模拟时，才有意义；热力学模拟无意义)
    int nFoldTol   = 0;
    int nunFoldTol = 0;
    // Average steps of folding or unbolding (发生折叠或解折叠的平均步数)
    int aveNadim   = 0;

    read_initalconform(param.initialconform_filename);
    //=======================Initial Work End=======================/

    // TL: nOutputCount is used till end of main loop, what's it used for?
    int nOutputCount = 0;

    // !===================== simulation start =====================
    // outermost loop
    // fortran code: do 2000 nConfcount=nCon0+1, nConform
    for (int nConfcount = param.nCon0; nConfcount < param.nConform;
            nConfcount++) {
        // read the intial conformation
        // (in word, reverse initial and native conformation filenames) 
        // fortran code: read(20,*) xinit(i),yinit(i),zinit(i)
        read_initalconform(param.initialconform_filename);

        // Copy from fortran code, TL: means?
        // Utility counter for all trajectories
        // (对得到的所有轨迹进行统计：发生/未发生的百分比)
        int nFoldSub = 0;
        int nunFoldSub = 0;

        // Fortran code: do 1500 nRuncount=1, nRunConf
        for (int nRuncount = 0; nRuncount < param.nRunConf; nRuncount++) {
            // copy intial conformation to native conformation
            for (int j = 0; j < param.nparttol; j++) {
                particle_list[j].x = particle_list[j].xinit;
                particle_list[j].y = particle_list[j].yinit;
                particle_list[j].z = particle_list[j].zinit;
            }
            origin_adjust();
            InitVel(enerkin);
            // The following variable is used before definition, fortran code:
            // call force(e_pot,e_unbond_tot,e_bind_tot,e_tors_tot,e_bend_tot,
            // e_bond_tot)
            double e_pot        = 0.0;
            double e_bond_tot   = 0.0;
            double e_bend_tot   = 0.0;
            double e_tors_tot   = 0.0;
            double e_unbond_tot = 0.0;
            double e_bind_tot   = 0.0;
            force.force(particle_list, e_pot, e_unbond_tot, e_bind_tot,
                    e_tors_tot, e_bend_tot, e_bond_tot);
            // Generate a random coordinates
            RANTERM();

            for (int kl = 0; kl < param.npartM; kl++) {
                particle_list[kl].fxo     = particle_list[kl].fx;
                particle_list[kl].fyo     = particle_list[kl].fy;
                particle_list[kl].fzo     = particle_list[kl].fz;
                particle_list[kl].frandxo = particle_list[kl].frandx;
                particle_list[kl].frandyo = particle_list[kl].frandy;
                particle_list[kl].frandzo = particle_list[kl].frandz;
            }

            // Copy from fortran code, TL: means?
            nOutputCount = 0;
            // Control the output steps (控制在哪一步进行输出)
            int nOutputGap = -1;

            // Counter for the following loop (对下面循环进行计数)
            param.nadim = -1; // fortran code: nadim=-1
            // ! main dynamics cycle
            // fortran code: do 100 nadim_new=0,nstep
            for (int j = 0; j < param.nstep; j++) {
                param.nadim += 1;
                verlet(enerkin, e_pot);

                // !-------------save data for purpose of restore-------------
                // Because nadim_old and nOutGap_old is needed in the next
                // "if statement", so I bring them out of the "if statement"
                int nadim_old = param.nadim;
                int nOutGap_old = nOutputGap;
                // if(mod(nadim,nsnap).eq.0.and.vx(1).ge.-1e6.and.vx(1).le.1e6)
                // because particle_list's index startswith 0, so v(1)-->[0]
                if (param.nadim % param.nsnap == 0
                        && particle_list[0].vx >= -1e6
                        && particle_list[0].vx <= 1e6) {
                    for (int kl = 0; kl < param.nparttol; kl++) {
                        particle_list[kl].xsave = particle_list[kl].x;
                        particle_list[kl].ysave = particle_list[kl].y;
                        particle_list[kl].zsave = particle_list[kl].z;
                    }
                    // TL: It seems that nadim_old and nOutGap_old are unused.
                    nadim_old = param.nadim;
                    nOutGap_old = nOutputGap;
                }
                // ! restore, TL: why restore?
                if (!(particle_list[0].vx >= -1e6
                            && particle_list[0].vx <= 1e6)) {
                    for (int i = 0; i < param.nparttol; i++) {
                        particle_list[i].x = particle_list[i].xsave;
                        particle_list[i].y = particle_list[i].ysave;
                        particle_list[i].z = particle_list[i].zsave;
                    }

                    InitVel(enerkin);
                    force.force(particle_list, e_pot, e_unbond_tot, e_bind_tot,
                            e_tors_tot, e_bend_tot, e_bond_tot);
                    // Related to Randon number generator
                    RANTERM();

                    for (int kl = 0; kl < param.npartM; kl++) {
                        particle_list[kl].fxo = particle_list[kl].fx;
                        particle_list[kl].fyo = particle_list[kl].fy;
                        particle_list[kl].fzo = particle_list[kl].fz;
                        particle_list[kl].frandxo = particle_list[kl].frandx;
                        particle_list[kl].frandyo = particle_list[kl].frandy;
                        particle_list[kl].frandzo = particle_list[kl].frandz;
                    }
                    param.nadim = nadim_old;
                    nOutputGap = nOutGap_old;
                }
                // !-------------save data for purpose of restore end---------- 

                // !--------------output conformation-------------
                nOutputGap += 1;
                output_conformation(nOutputGap, nOutputCount);
                // !--------------output conformation end---------

                // To check whether end the loop
                if ((param.gQ_f1 <= param.gQf1denatural
                            && param.gQ_f2 <= param.gQf2denatural
                            && param.gQ_b <= param.gQbdenatural)
                        || (param.gQ_f1 >= param.gQf1nativ
                            && param.gQ_f2 >= param.gQf2nativ
                            && param.gQ_b >= param.gQbnativ)) {
                    break;
                }
            } // fortran code: 100 continue

            // Accumulate the number of folding (动力学模拟情况下: 折叠计数)
            aveNadim += param.nadim;
            int k = 0;
            if(!(particle_list[0].vx >= -1e8 && particle_list[0].vx <= 1e8)) {
                k = 0;
            } else if (param.gQ_f1 <= param.gQf1denatural
                    && param.gQ_f2 <= param.gQf2denatural
                    && param.gQ_b <= param.gQbdenatural) {
                k = -1;
                nunFoldSub += 1;
                nunFoldTol += 1;
            } else if (param.gQ_f1 >= param.gQf1nativ
                    && param.gQ_f2 >= param.gQf2nativ
                    && param.gQ_b >= param.gQbnativ) {
                k = 1;
                nFoldSub += 1;
                nFoldTol += 1;
            }
            msg = log.format("%5d %5d %5d %13d %10.3f %10.3f %10.3f\n",
                    nConfcount, nRuncount, k, param.nadim, param.gQ_b,
                    param.gQ_f1, param.gQ_f2);
            log.info(output_filenames[28], msg);
        } // end of fortran code: 1500     continue

        msg = log.format("# Conform No.:%5d, Folding times: %5d, Transmission \
                coefficient: %10.3f, unFolding times: %5d\n",
                nConfcount, nFoldSub, 1.0 * nFoldSub / param.nRunConf,
                nunFoldSub);
        log.info(output_filenames[28], msg);
    } // end of the outermost loop, fortran code: 2000 continue 

    aveNadim = aveNadim / ((param.nConform - param.nCon0) * param.nRunConf);
    msg = log.format("# Total folding times: %5d, Total unfolding times: %4d\n",
            nFoldTol, nunFoldTol);
    log.info(output_filenames[28], msg);
    msg = log.format("# Average Running Time: %14.4f\n", aveNadim);
    log.info(output_filenames[28], msg);
    // !===================== simulation end =====================

    // ! output target
    write_target(nOutputCount);
    write_histogram();
    return 0;
}

int pf::Simulation::RANTERM() {
    for (int i = 0; i < param.npartM; i++) {
        particle_list[i].frandx = gauss(xsi) * param.randconst;
        particle_list[i].frandy = gauss(xsi) * param.randconst;
        particle_list[i].frandz = gauss(xsi) * param.randconst;
    }

    return 0;
}

/**
 * verlet: one sort of dynamics model
 * 分子动力学模拟算法
 */
int pf::Simulation::verlet(double &enerkin, double &e_pot) {
    enerkin = 0.0;
    RANTERM();

    // Calculate coordinates
    for (int i = 0; i < param.npartM; i++) {
        particle_list[i].x += (param.dt * particle_list[i].vx + pow(param.dt, 2)
                * (particle_list[i].fxo + particle_list[i].frandxo -
                    param.gm * particle_list[i].vx) / 2.0);
        particle_list[i].y += (param.dt * particle_list[i].vy + pow(param.dt, 2)
                * (particle_list[i].fyo + particle_list[i].frandyo -
                    param.gm * particle_list[i].vy) / 2.0);
        particle_list[i].z += (param.dt * particle_list[i].vz + pow(param.dt, 2)
                * (particle_list[i].fzo + particle_list[i].frandzo -
                    param.gm * particle_list[i].vz) / 2.0);
    }
    e_pot = 0.0;
    double e_bond_tot = 0.0;
    double e_bend_tot = 0.0;
    double e_tors_tot = 0.0;
    double e_unbond_tot = 0.0;
    double e_bind_tot = 0.0;
    force.force(particle_list, e_pot, e_unbond_tot, e_bind_tot, e_tors_tot,
            e_bend_tot, e_bond_tot);

    // Calculate verlocity
    for (int i = 0; i < param.npartM; i++) {
        particle_list[i].vx = param.c_1 * particle_list[i].vx + param.c_2 * 
            (particle_list[i].fxo + particle_list[i].frandxo +
             particle_list[i].fx + particle_list[i].frandx);
        particle_list[i].vy = param.c_1 * particle_list[i].vy + param.c_2 * 
            (particle_list[i].fyo + particle_list[i].frandyo +
             particle_list[i].fy + particle_list[i].frandy);
        particle_list[i].vz = param.c_1 * particle_list[i].vz + param.c_2 * 
            (particle_list[i].fzo + particle_list[i].frandzo +
             particle_list[i].fz + particle_list[i].frandz);
        enerkin += pow(particle_list[i].vx, 2) + pow(particle_list[i].vy, 2) +
            pow(particle_list[i].vz, 2);
    }

    for (int i = 0; i < param.npartM; i++) {
        particle_list[i].fxo = particle_list[i].fx;
        particle_list[i].fyo = particle_list[i].fy;
        particle_list[i].fzo = particle_list[i].fz;
        particle_list[i].frandxo = particle_list[i].frandx;
        particle_list[i].frandyo = particle_list[i].frandy;
        particle_list[i].frandzo = particle_list[i].frandz;
    }

    // Calculate kinetic energy
    enerkin = 0.50 * enerkin * param.amass;
    // It seems that tempins is unused.
    // double tempins = 2.0 * enerkin / (3.0 * param.npartM * param.boltz);

    // !  periodic boundary condiction is used based on the coordinate of 
    // residue 15
    if (param.iFlagMov == 2 && param.nadim % 20 == 0) origin_adjust();
    if (param.nadim % 20 == 0) pbc_shift();

    // !  histogram calculation 
    if(param.nadim >= param.nbinsnap0 && param.nadim % param.nbinsnap == 0) {
        int ibin_f = (param.gQ_f - param.vbin0) / param.dbin_f + 1;
        if (ibin_f <= 0) ibin_f = 1;
        if (ibin_f > param.nbin_f) ibin_f= param.nbin_f;

        int ibin_f1 = (param.gQ_f1 - param.vbin0) / param.dbin_f + 1;
        if (ibin_f1 <= 0) ibin_f1 = 1;
        if (ibin_f1 > param.nbin_f) ibin_f1 = param.nbin_f;

        int ibin_f2 = (param.gQ_f2 - param.vbin0) / param.dbin_f + 1;
        if (ibin_f2 <= 0) ibin_f2 = 1;
        if (ibin_f2 > param.nbin_f) ibin_f2 = param.nbin_f;

        int ibin_b = (param.gQ_b - param.vbin0) / param.dbin_b + 1;
        if (ibin_b <= 0) ibin_b = 1;
        if (ibin_b > param.nbin_b) ibin_b = param.nbin_b;

        int iEbin = (e_pot - param.vEbin0) / param.dEbin + 1;
        if (iEbin <= 0) iEbin = 1;
        if (iEbin > param.nEbin) iEbin = param.nEbin;

        int iEbbin = (e_bind_tot - param.vEbbin0) / param.dEbbin + 1;
        if (iEbbin <= 0) iEbbin = 1;
        if (iEbbin > param.nEbbin) iEbbin = param.nEbbin;

        int iRbin = (param.R - param.vRbin0) / param.dRbin + 1;
        if (iRbin <= 0) iRbin = 1;
        if (iRbin > param.nRbin) iRbin = param.nRbin;

        int iWbin = (param.gQ_w - param.vWbin0) / param.dWbin + 1;
        if (iWbin <= 0) iWbin = 1;
        if (iWbin > param.nWbin) iWbin = param.nWbin;
        // fortran code: if(gQ_b.le.cri_Qb) then iQb_cri=1 else iQb_cri=2 endif
        int iQb_cri = 2;
        if(param.gQ_b <= param.cri_Qb) iQb_cri=1;

        // !  for PFBbin:
        statis.PFbin[ibin_f] = statis.PFbin[ibin_f] + 1;
        statis.PBbin[ibin_b] = statis.PBbin[ibin_b] + 1;
        statis.PFBbin[ibin_f][ibin_b] = statis.PFBbin[ibin_f][ibin_b] + 1;
        statis.PFBEbin[ibin_f][ibin_b][iEbin] =
            statis.PFBEbin[ibin_f][ibin_b][iEbin] + 1;
        statis.PQbEbbin[ibin_b][iEbbin] = statis.PQbEbbin[ibin_b][iEbbin] + 1;
        statis.PFBRbin[ibin_f][ibin_b][iRbin] = 
            statis.PFBRbin[ibin_f][ibin_b][iRbin] + 1;
        statis.PFFBbin[ibin_f1][ibin_f2][ibin_b] = 
            statis.PFFBbin[ibin_f1][ibin_f2][ibin_b] + 1;

        statis.PQwbin[iWbin] = statis.PQwbin[iWbin] + 1;
        statis.PQwQbbin[iWbin][ibin_b] = statis.PQwQbbin[iWbin][ibin_b] + 1;
        statis.PQwRbin[iWbin][iRbin][iQb_cri] =
            statis.PQwRbin[iWbin][iRbin][iQb_cri] + 1;
        statis.PQwEbbin[iWbin][iEbbin][iQb_cri] = 
            statis.PQwEbbin[iWbin][iEbbin][iQb_cri] + 1;
        statis.PQwQfbin[iWbin][ibin_f1] = statis.PQwQfbin[iWbin][ibin_f1] + 1;
    }

    // Write the intermediate results
    string msg = "nadim\tgQ_f\tgQ_b\tgQ_w\tE_k\tE_pot\tE_b\tE_bind\teGr\tR\n";
    if (param.nadim == 0)
        log.info(output_filenames[9], msg);
    if (param.nadim == 0)
        cout << msg;
    msg = log.format("%13d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f \
        %10.3f %10.3f \n", param.nadim, param.gQ_f, param.gQ_b, param.gQ_w,
        enerkin, e_pot, e_bond_tot, e_unbond_tot, e_bind_tot, param.eGr,
        param.R);
    if (param.nadim % param.nsnap == 0)
        log.info(output_filenames[9], msg);
    if (param.nadim % param.nsnap == 0)
        cout << msg;
    return 0;
}

int pf::Simulation::pbc_shift() {
    // fortran code: if(npart1.le.0) goto 400
    if(param.npart1 <= 0)
        return 0;

    // ! note that the chain 2 is located at the center
    double x00 = 0.0;
    double y00 = 0.0;
    double z00 = 0.0;
    for (int i = 0; i < param.npart1; i++) {
        x00 += particle_list[i].x;
        y00 += particle_list[i].y;
        z00 += particle_list[i].z;
    }
    x00 = x00 / param.npart1;
    y00 = y00 / param.npart1;
    z00 = z00 / param.npart1;

    // !  R is the distance to the center (0,0,0) of the box to center of
    // chain A (x00,y00,z00)
    param.R = sqrt(pow(x00, 2) + pow(y00, 2) + pow(z00, 2));

    // x coordinate setting
    int ixflag = 0;
    if (x00 > ( param.pL / 2.0 + param.dl)) {
        ixflag = 1;
        for (int i = 0; i < param.npart1; i++) {
            particle_list[i].x -= param.pL;
        }
    }
    // fortran code: if ( ixflag == 1 ) goto 111
    if ( ixflag != 1 && x00 < (-1 * param.pL / 2.0 - param.dl)) {
        for (int i = 0; i < param.npart1; i++) {
            particle_list[i].x += param.pL;
        }
    }

    // y coordinate setting
    int iyflag = 0;
    if (y00 > (param.pL / 2.0 + param.dl)) {
        iyflag = 1;
        for (int i = 0; i < param.npart1; i++) {
            particle_list[i].y -= param.pL;
        }
    }
    if (iyflag != 1 && y00 < (-1 * param.pL / 2.0 - param.dl)) {
        for (int i = 0; i < param.npart1; i++) {
            particle_list[i].y += param.pL;
        }
    }

    // z coordinate setting
    int izflag = 0;
    if (z00 > (param.pL / 2.0 + param.dl)) {
        izflag = 1;
        for (int i = 0; i < param.npart1; i++) {
            particle_list[i].z -= param.pL;
        }
    }
    if (izflag != 1 && z00 < (-1 * param.pL / 2.0 - param.dl)) {
        for (int i = 0; i < param.npart1; i++) {
            particle_list[i].z += param.pL;
        }
    }
    return 0;
}

int pf::Simulation::output_conformation(int &nOutputGap, int &nOutputCount) {
    string msg = "";
    if (param.nadim == 0 ) {
        msg = to_string(nOutputCount) + '\t' + to_string(param.nadim) + '\t' +
            to_string(param.gQ_f) + '\t' + to_string(param.gQ_b) + '\t' +
            to_string(param.R);
        log.info(output_filenames[27],
                "nOutputCount,  nadim,   gQ_f,     gQ_b,     R");
    }

    if (nOutputCount < param.nConformOutput && param.nadim >= param.nOutput0
            && nOutputGap >= param.ndOutput && param.gQ_f1 >= param.outQf1_i
            && param.gQ_f1 <= param.outQf1_f && param.gQ_f2 >= param.outQf2_i
            && param.gQ_f2 <= param.outQf2_f && param.gQ_b >= param.outQb_i
            && param.gQ_b <= param.outQb_f) {
        nOutputGap = 0;

        write_mol(nOutputCount);
        msg = log.format("%7d %13d %10.3f %10.3f %10.3f\n", nOutputCount,
                param.nadim, param.gQ_f, param.gQ_b, param.R);
        log.info(output_filenames[27], msg);

        for (int j = 0; j < param.nparttol; j++) {
            msg = to_string(particle_list[j].x) + '\t' +
                to_string(particle_list[j].y) + '\t' +
                to_string(particle_list[j].z);
            log.info(output_filenames[26], msg + '\n');
        }

        // fortran code: write(26,'('' '')')
        log.info(output_filenames[26], " ");
        // write(*,'(''n='',i3,'', nadim='',i8,'', gQ='',f7.3, '',
        // old gQ='', 1i3)')nOutputCount,nadim,gQ,natcont
        // gQ and natcont don't define, so we define them as follos
        double gQ = 0.0;
        int natcont = 0;
        msg = log.format("n=%3d nadim=%4d gQ=%7.3f, old gQ=%3d\n",
                nOutputCount, param.nadim, gQ, natcont);
        cout << msg;
    }

    return 0;
}

int pf::Simulation::write_mol(int &nOutputCount) {
    nOutputCount = nOutputCount + 1;

    string msg = log.format("Molecule-%7d\n", nOutputCount);
    log.info(output_filenames[22], msg);
    log.info(output_filenames[22],
            " ViewerPro         3D                             0\n");

	if(param.iFlagMov == 1 || param.iFlagMov == 0) {
        msg = log.format("%3d %3d  0  0  0  0  0  0  0  0999 V2000\n",
                param.npart1, param.npart1 - 1);
    } else {
        msg = log.format("%3d %3d  0  0  0  0  0  0  0  0999 V2000\n",
                param.nparttol, param.nparttol - 1);
    } 
    log.info(output_filenames[22], msg);

    for (int i = 0; i < param.npartM; i++) {
        msg = log.format("%10.4f %10.4f %10.4f C   0  0  0  0  0  0  0  0  \
                0  0\n", particle_list[i].x, particle_list[i].y,
                particle_list[i].z);
    }
    for (int i = 0; i < param.npartM - 1; i++) {
        int j = i + 1; 
        if (j <= param.npart1 || i > param.npart1) {
            msg = log.format("%3d %3d  1  0  0  0\n", i, j);
            log.info(output_filenames[22], msg);
        }
    }
    log.info(output_filenames[22], "M  END\n");
	log.info(output_filenames[22], "$$$$\n");

    return 0;
}

// Viewer software: Pymol
int pf::Simulation::write_target(int &nOutputCount) {
    nOutputCount += 1;
	string msg = log.format("Molecule-%7d\n", nOutputCount);
    log.info(output_filenames[22], msg);
	msg = log.format("Molecule-%7d\n",	nOutputCount);
    log.info(output_filenames[22], msg);
	msg = log.format("%3d %3d  0  0  0  0  0  0  0  0999 V2000\n",
            nOutputCount);
    log.info(output_filenames[22], msg);

    for (int i = param.npart1 + 1; i < param.nparttol; i++) {
        int j = i + 1;
        msg = log.format("%3d %3d  1  0  0  0", i, j);
        log.info(output_filenames[22], msg);
    }
    log.info(output_filenames[22], "M  END");
    log.info(output_filenames[22], "$$$$");

    return 0;
}

int pf::Simulation::write_histogram() {
    string msg = "";
    // !  histogram result saving to file
    // BUG: i starts from 1 in original code, but now it starts from 0, so we 
    // add 1 to the i-related compuation
    for (int i = 0; i < param.nbin_f; i++) {
        string msg = log.format("%6.2f, %12.1f\n",
                param.vbin0 + (i + 1 - 0.5) * param.dbin_f, statis.PFbin[i]);
        log.info(output_filenames[15], msg);
    }
    
    // fortran code: if (iFlagMov.le.0) goto 25
    // BUG: i starts from 1 in original code, but now it starts from 0
    if (param.iFlagMov > 0) {
        for (int i = 0; i < param.nbin_b; i++) {
            double sumQ_b = 0;
            for (int j = 0; j < param.nbin_f; j++) {
		        sumQ_b += statis.PFBbin[j][i];
            }
            msg = log.format("%6.2f, %12.1f\n",
                    param.vbin0 + (i + 1 - 0.5) * param.dbin_b, sumQ_b);
            log.info(output_filenames[14], msg);
        }
    }

    // do 1003 i=1, nbin_f
    //   do 1004 j=1, nbin_b        
    for (int i = 0; i < param.nbin_f; i++) {
        for (int j = 0; j < param.nbin_b; j++) {
            if (param.iFlagMov > 0 && statis.PFBbin[i][j] > 1e-5) {
                msg = log.format("%7.2f %7.2f, %12.1f\n",
                        param.vbin0 + (i + 1 - 0.5) * param.dbin_f,
                        param.vbin0 + (j + 1 - 0.5) * param.dbin_b,
                        statis.PFBbin[i][j]);
                log.info(output_filenames[19], msg);
            }
            // fortran code: if ( IsEbin == 0 ) goto 1007
            if (param.IsEbin != 0 ) {
                for (int k = 0; k < param.nEbin; k++) {
                    if (statis.PFBEbin[i][j][k] > 1e-5) {
                        msg = log.format("%7.2f %7.2f, %8.2f, %12.1f",
                                param.vbin0 + (i + 1 - 0.5) * param.dbin_f,
                                param.vbin0 + (j + 1 - 0.5) * param.dbin_b,
                                param.vEbin0 + (k + 1 - 0.5) * param.dEbin,
                                statis.PFBEbin[i][j][k]);
                        log.info(output_filenames[23], msg);
                    } 
                } // 1005 continue
            } // 1007 continue

            // fortran code: if ( IsRbin == 0 ) goto 1008
            if (param.IsRbin != 0 ) {
                for (int l = 0; l < param.nRbin; l++) {
                    if(statis.PFBRbin[i][j][l] > 1e-5) {
                        msg = log.format("%7.2f %7.2f, %8.2f, %12.1f",
                                param.vbin0 + (i + 1 - 0.5) * param.dbin_f,
                                param.vbin0 + (j + 1 - 0.5) * param.dbin_b,
                                param.vRbin0 + (l + 1 - 0.5) * param.dRbin,
                                statis.PFBRbin[i][j][l]);
                        log.info(output_filenames[24], msg);
                    } 
                } // 1006 continue
            } // 1008 continue

            // fortran code: if ( iFlagMov .ne. 2 ) goto 1010
            if (param.iFlagMov == 2) {
                for (int l = 0; l < param.nbin_f; l++) {
                    if (statis.PFFBbin[i][l][j] > 1e-5) {
                        msg = log.format("%7.2f %7.2f, %8.2f, %12.1f",
                                param.vbin0 + (i + 1 - 0.5) * param.dbin_f,
                                param.vbin0 + (l + 1 - 0.5) * param.dbin_f,
                                param.vbin0 + (j +1 - 0.5) * param.dbin_b,
                                statis.PFFBbin[i][l][j]);
                        log.info(output_filenames[25], msg);
                    }
                } // 1009 continue  
            } // 1010 continue

        } // 1004 continue
    } // 1003 continue

    // if ( iFlagMov.le.0 .or. IsEbbin.eq.0 ) goto 40
    if (param.iFlagMov > 0 && param.IsEbbin != 0 ) {
        for (int j = 0; j < param.nbin_b; j++) {
            for (int k = 0; k < param.nEbbin; k++) {
                if(statis.PQbEbbin[j][k] > 1e-5) {
                    msg = log.format("%8.2f %8.2f, %12.1f",
                            param.vbin0 + (j + 1 - 0.5) * param.dbin_b,
                            param.vEbbin0 + (k + 1 - 0.5) * param.dEbbin,
                            statis.PQbEbbin[j][k]);
                    log.info(output_filenames[29], msg);
                } 
            } // 30 continue
        } // 35 continue
    } // 40 continue

    // fortran code: if ( iFlagMov.le.0 .or. IsWbin.eq.0 ) goto 180
    if (param.iFlagMov > 0 && param.IsWbin != 0 ) {
        for (int i = 0; i < param.nWbin; i++) {
            msg = log.format("%8.2f, %12.1f", param.vWbin0 +
                    (i + 1 - 0.5) * param.dWbin, statis.PQwbin[i]);
            log.info(output_filenames[31], msg);
            
            for (int j = 0; j < param.nbin_b; j++) {
		        if (statis.PQwQbbin[i][j] > 1e-5) { 
                    msg = log.format("%8.2f %8.2f, %12.1f",
                            param.vWbin0 + (i + 1 - 0.5) * param.dWbin,
                            param.vbin0 + (j + 1 - 0.5) * param.dbin_b,
                            statis.PQwQbbin[i][j]);
                    log.info(output_filenames[31], msg);
                }
            }
            for (int l = 0; l < param.nRbin; l++) {
		        if (statis.PQwRbin[i][l][0] > 1e-5 || 
                        statis.PQwRbin[i][l][1] > 1e-5) {
                    // write(32, '(2f8.2, 2f12.1)') vWbin0+(i-0.5)*dWbin,
                    // vRbin0+(l-0.5)*dRbin, PQwRbin(i,l,1), PQwRbin(i,l,2)
                    // PQwRbin(i,l,1) --> PQwRbin(i,l,0)
                    msg = log.format("%8.2f %8.2f, %12.1f %12.1f",
                            param.vWbin0 + (i + 1 - 0.5) * param.dWbin,
                            param.vRbin0 + (l + 1 - 0.5) * param.dRbin,
                            statis.PQwRbin[i][l][0], statis.PQwRbin[i][l][1]);
                    log.info(output_filenames[31], msg);
                }
            }
            for (int k = 0; k < param.nEbbin; k++) {
		        if (statis.PQwEbbin[i][k][0] > 1e-5 ||
                        statis.PQwEbbin[i][k][1] > 1e-5) {
                    msg = log.format("%8.2f %8.2f, %12.1f %12.1f",
                            param.vWbin0 + (i + 1 - 0.5) * param.dWbin,
                            param.vEbbin0 + (k + 1 - 0.5) * param.dEbbin,
                            statis.PQwEbbin[i][k][0], statis.PQwEbbin[i][k][1]);
                    log.info(output_filenames[33], msg);
                }
            }
            for (int j = 0; j < param.nbin_f; j++) {
                if (statis.PQwQfbin[i][j] > 1e-5) {
                    msg = log.format("%8.2f %8.2f, %12.1f %12.1f",
                            param.vWbin0 + (i + 1 - 0.5) * param.dWbin,
                            param.vbin0 + (j + 1 - 0.5) * param.dbin_f,
                            statis.PQwQfbin[i][j]);
                    log.info(output_filenames[34], msg);
                }
            }
        } // 150 continue
    } // 180 continue
    return 0;
}

/**
 * conformation to coordinate system
 * 对体系进行对齐
 */
int pf::Simulation::origin_adjust() {
    if (param.nparttol - param.npart1 <= 0) return 0;

    double x00 = 0.0;
    double y00 = 0.0;
    double z00 = 0.0;

    // fortran code: do 281 j=npart1+1,nparttol
    //double xparticle[param.nparttol];
    for (int j = param.npart1; j < param.nparttol; j++) {
        // xparticle[j] is never be used
        //x00 = x00/* + xparticle[j] */;
        x00 = x00 + particle_list[j].x;
        y00 = y00 + particle_list[j].y;
        z00 = z00 + particle_list[j].z;
    }

    x00 /= (param.nparttol - param.npart1);
    y00 /= (param.nparttol - param.npart1);
    z00 /= (param.nparttol - param.npart1);


    for (int j = 0; j < param.nparttol; j++) {
        particle_list[j].x -= x00;
        particle_list[j].y -= y00;
        particle_list[j].z -= z00;
    }

    return 0;
}

int pf::Simulation::InitVel(double &enerkin) {
    double sumvx = 0.0;
    double sumvy = 0.0;
    double sumvz = 0.0;
    double vm = sqrt(param.temp * param.boltz / param.amass);

    for (int i = 0; i < param.npartM; i++) {
        particle_list[i].vx = vm * gauss(xsi);
        particle_list[i].vy = vm * gauss(xsi);
        particle_list[i].vz = vm * gauss(xsi);
        sumvx += particle_list[i].vx;
        sumvy += particle_list[i].vy;
        sumvz += particle_list[i].vz;
    }

    sumvx /= param.npartM;
    sumvy /= param.npartM;
    sumvz /= param.npartM;
    enerkin = 0.0;

    for (int i = 0; i < param.npartM; i++) {
        particle_list[i].vx -= sumvx;
        particle_list[i].vy -= sumvy;
        particle_list[i].vz -= sumvz;
        enerkin += (pow(particle_list[i].vx, 2) + pow(particle_list[i].vy, 2)
                    + pow(particle_list[i].vz, 2));
    }
    enerkin =  0.5 * enerkin * param.amass;
    // tempins is unused.
    // double tempins = 2.0 * enerkin / (3.0 * param.npartM * param.boltz);
    return 0;
}

double pf::Simulation::gauss(double xsi) {
    double r = 0.0;
    for (int j = 1; j < 13; j++) {
        r += rannyu();
    }

    r = (r - 6.0) / 4.0;
    double r2 = pow(r, 2);
    double gauss = ((((a9 * r2 + a7) * r2 + a5) * r2 + a3) * r2 + a1) * r;

    return gauss;
}

// I directly move the body of this function into the rannyu()
int pf::Simulation::rnyubd() {
    m[0] = 502, m[1] = 1521, m[2] = 4071, m[3] = 2107;
    l[0] = 0, l[1] = 0, l[2] = 0, l[3] = 1;
    return 0;
}

double pf::Simulation::rannyu() {
    // fortran code: common /rnyucm/ m(4),l(4), TL: when do they execute
    // TL: Is there any relationship between m/l(4) with m/l1-4?
    double ooto12 = 1.0 / 4096.0;
    int itwo12 = 4096;
    // Initialize m1-4 l1-4 in setrn()
    // int m1 = 0, m2 = 0, m3 = 0, m4 = 0, l1 = 0, l2 = 0, l3 = 0, l4 = 0;
    //////////////////Body of rnyubd()/////////////////////
    m[0] = 502, m[1] = 1521, m[2] = 4071, m[3] = 2107;
    // int l1 = 0, l2 = 0, l3 = 0, l4 = 1;
    // ///////////////////////////////////////////////////////
    // int i1 = l1 * m4 + l2 * m3 + l3 * m2 + l4 * m1;
    // int i2 = l2 * m4 + l3 * m3 + l4 * m2;
    // int i3 = l3 * m4 + l4 * m3;
    // int i4 = l4 * m4;
    // l4 = fmod(i4, itwo12);
    // i3 = i3 + i4 / itwo12;
    // l3 = fmod(i3, itwo12);
    // i2 = i2 + i3 / itwo12;
    // l2 = fmod(i2, itwo12);
    // l1 = fmod((i1 + i2 / itwo12), itwo12);
    // double rannyu = ooto12 * (float(l1) +
    //         ooto12 * (float(l2) + ooto12 * (float(l3) + ooto12 * (float(l4)))));
    int i1 = l[0] * m[3] + l[1] * m[2] + l[2] * m[1] + l[3] * m[0];
    int i2 = l[1] * m[3] + l[2] * m[2] + l[3] * m[1];
    int i3 = l[2] * m[3] + l[3] * m[2];
    int i4 = l[3] * m[3];
    l[3] = i4 % itwo12;
    i3 = i3 + i4 / itwo12;
    l[2] = i3 % itwo12;
    i2 = i2 + i3 / itwo12;
    l[1] = i2 % itwo12;
    l[0] = (i1 + i2 / itwo12) % itwo12;
    double rannyu = ooto12 * (float(l[0]) +
            ooto12 * (float(l[1]) + ooto12 * (float(l[2]) + ooto12 * (float(l[3])))));
    return  rannyu;
}

int pf::Simulation::setrn(int iseed[4]) {
    int isn, ipe, ipd, id;
    for (int j = 0; j < 4; j++) {
        isn = 0;
        for (int i = 0; i < 4; i++) {
            // fortran code: ipe = 4 - i
            ipe = 4 - (i + 1);
            ipd = pow(10, ipe);
            id = iseed[j] / ipd;
            isn += id * pow(8, ipe);
            iseed[j] -= id * ipd;
        }
        iseed[j] = isn;
    }

    for (int i = 0; i < 4; i++) {
        l[i] = iseed[i];
    }
    // fortran code: l(4)=2*(l(4)/2)+1
    l[3] = 2 * (l[3] / 2) + 1;
    return 0;
}

// Never be used.
int pf::Simulation::savern(int iseed[4]) {
    for (int i = 0; i < 4; i++) {
        iseed[i] = l[i];
    }
    for (int j = 0; j < 4; j++) {
        int isn = 0;
        for (int i = 0; i < 4; i++) {
            // fortran code: ipe = 4 - i
            int ipe = 4 - (i + 1);
            int ipo = pow(8, ipe);
            int id = iseed[j] / ipo;
            isn += id * pow(10, ipe);
            iseed[j] -= ipo * id;
        }
        iseed[j] = isn; 
    }
    return 0;
}

// Unused function
double pf::Simulation::calculate_gyrationradius() {
    double ctx = 0;
    double cty = 0;
    double ctz = 0;
    double anr = 0;
    for (int i = 0; i < param.nparttol; i++) {
        ctx += particle_list[i].x;
        cty += particle_list[i].y;
        ctz += particle_list[i].z;
        anr += (pow(particle_list[i].x, 2) + pow(particle_list[i].y, 2) +
            pow(particle_list[i].z, 2));
    }
    ctx /= param.nparttol;
    cty /= param.nparttol;
    ctz /= param.nparttol;
    double gyr1 = anr / param.nparttol;
    double gyr = sqrt(gyr1 - (pow(ctx, 2) + pow(cty, 2) + pow(ctz, 2)));       

    return gyr;
}

int pf::Simulation::nativeinformation() {
    param.nQnative_f1 = 0;
    param.nQnative_f2 = 0;
    param.nQnative_b = 0;

    // kunbond[i]: indicates whether exits interaction between tow residues
    // (标记两个残基之间是否有相互作用)
    double xij, yij, zij;
    for (int k = 0; k < param.nunbond; k++) {
        int i = param.iun[k] - 1;
        int j = param.jun[k] - 1;
        // fortran code: xij=(x(i)-x(j)) yij=(y(i)-y(j)) zij=(z(i)-z(j))
        // Because appNCS.dat' data starts from 1, i-->i-1, j-->j-1
        xij = particle_list[i].x - particle_list[j].x;
        yij = particle_list[i].y - particle_list[j].y;
        zij = particle_list[i].z - particle_list[j].z;
        param.runbond_nat[k] = sqrt(pow(xij, 2) + pow(yij, 2) +
                pow(zij, 2));
        if (param.kunbond[k] == 1) {
            if (i <= param.npart1 && j <= param.npart1) {
                param.nQnative_f1 = param.nQnative_f1 + 1;
            } else if (i > param.npart1 && j > param.npart1) {
                param.nQnative_f2 = param.nQnative_f2 + 1;
            } else {
                param.nQnative_b = param.nQnative_b + 1;
            }
        }
    }

    // bond length (键长)
    for (int i = 0; i < param.nparttol - 1; i++) {
        int j = i + 1;
        xij = particle_list[j].x - particle_list[i].x;
        yij = particle_list[j].y - particle_list[i].y;
        zij = particle_list[j].z - particle_list[i].z;
        param.rbond_nat[i] = sqrt(pow(xij, 2) + pow(yij, 2) +
                pow(zij, 2));
    }

    // bond angle (键角)
    double xkj, ykj, zkj, rij, rkj, costheta;
    for (int i = 0; i < param.nparttol - 2; i++) {
        int j = i + 1;
        int k = i + 2;
        xij = particle_list[i].x - particle_list[j].x;
        yij = particle_list[i].y - particle_list[j].y;
        zij = particle_list[i].z - particle_list[j].z;

        xkj = particle_list[k].x - particle_list[j].x;
        ykj = particle_list[k].y - particle_list[j].y;
        zkj = particle_list[k].z - particle_list[j].z;

        rij = sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2));
        rkj = sqrt(pow(xkj, 2) + pow(ykj, 2) + pow(zkj, 2));

        if(rij < eps || rkj < eps) break;

        costheta = (xij * xkj + yij * ykj + zij * zkj) / (rij * rkj);
        // if(abs(costheta).gt.epstht) costheta=sign(epstht,costheta)
        if(abs(costheta) > epstht) {
            costheta = costheta > 0 ? epstht : (-1 * epstht);
        }

        param.theta_nat[i] = acos(costheta);
    }

    // dihedral angle (二面角)
    double xkl, ykl, zkl, rkl;
    for (int i = 0; i < param.nparttol - 3; i++) {
        int j = i + 1;
        int k = i + 2;
        int l = i + 3;

        xij = particle_list[i].x - particle_list[j].x;
        yij = particle_list[i].y - particle_list[j].y;
        zij = particle_list[i].z - particle_list[j].z;

        xkj = particle_list[k].x - particle_list[j].x;
        ykj = particle_list[k].y - particle_list[j].y;
        zkj = particle_list[k].z - particle_list[j].z;

        xkl = particle_list[k].x - particle_list[l].x;
        ykl = particle_list[k].y - particle_list[l].y;
        zkl = particle_list[k].z - particle_list[l].z;

        rij = sqrt(pow(xij, 2) + pow(yij, 2) + pow(zij, 2));
        rkj = sqrt(pow(xkj, 2) + pow(ykj, 2) + pow(zkj, 2));
        rkl = sqrt(pow(xkl, 2) + pow(ykl, 2) + pow(zkl, 2)); 

        double xmj = yij * zkj - ykj * zij;
        double ymj = zij * xkj - zkj * xij;
        double zmj = xij * ykj - xkj * yij;

        double xnk = ykj * zkl - ykl * zkj;
        double ynk = zkj * xkl - zkl * xkj;
        double znk = xkj * ykl - xkl * ykj;

        double xil = ymj * znk - ynk * zmj;
        double yil = zmj * xnk - znk * xmj;
        double zil = xmj * ynk - xnk * ymj;

        double rnk=sqrt(pow(xnk, 2) + pow(ynk, 2) + pow(znk, 2));
        double rmj=sqrt(pow(xmj, 2) + pow(ymj, 2) + pow(zmj, 2));
        if(pow(rnk, 2) < eps || pow(rmj, 2) < eps) break;

        double phi = (xnk * xmj + ynk * ymj + znk * zmj) / (rnk * rmj);

        phi = acos(phi);

        //fortran code: dihedral_nat(i)=sign(phi,xkj*xil+ykj*yil+zkj*zil)
        double tempval = xkj * xil + ykj * yil + zkj * zil;
        param.dihedral_nat[i] = tempval > 0 ? phi : (-1.0 * phi);
    }
    return 0;
}
/////////////////////////////////////////////////////////////////////////
/******************End of definition of class Simulation****************/
/////////////////////////////////////////////////////////////////////////

int main()
{
    pf::Simulation sim("input.1BE9.test.dat");
    sim.start_simulation();
    return 0;
}
