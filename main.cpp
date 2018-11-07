#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>       /* exp */
#include <algorithm>    // std::min
#include <complex>

#include "mt64.h"

using namespace std;

// useful information obtained here
// http://mathfaculty.fullerton.edu/mathews//n2003/montecarlo/MonteCarloMod/Links/MonteCarloMod_lnk_4.html
// and here
// http://mathfaculty.fullerton.edu/mathews//n2003/MonteCarloMod.html

double pi = M_PI;

class ZPrimeModel {
	
	public:
	
	// default constructor
	ZPrimeModel(string _model, double _mZP, int _gmZmode) {
		
		// mathemtical constants
		pi = M_PI;
		
		// using Hewett, Rizzo angle convention...
		double thetaE6 = 0;
		if (_model=="SQ")        thetaE6 = asin(3.0/8.0*sqrt(6.0));  // SQ model // this is wrong
		else if (_model=="I")    thetaE6 = -1.0*asin(sqrt(5./8.));   // I model
		else if (_model=="N")    thetaE6 = asin(-1.0/4.0);	         // N model
		else if (_model=="Psi")  thetaE6 = 0.0;	                     // Psi model
		else if (_model=="Chi")  thetaE6 = pi/2.;                    // Chi model
		else if (_model=="Eta")  thetaE6 = asin(sqrt(3./8.));        // Eta model
		else {
			cout << "error in model selection" << endl;
			exit(0);
		}

		// E6 charges are linear combination of Psi and Chi states
		double A = cos(thetaE6)/(2.0*sqrt(6.0));  // psi basis
		double B = sin(thetaE6)/(2.0*sqrt(10.0)); // chi basis
  
		// set E6 fermion left chiral charges
		gLPu = (1.0)*A-(-1.0)*B;
		gLPd = (1.0)*A-(-1.0)*B;
		gLPv = (1.0)*A-(3.0)*B;
		//gLPv =  (1.0)*A-(3.0)*B;  // with RH nu on
		gLPl = (1.0)*A-(3.0)*B;

		// set E6 fermion right chiral charges
		gRPu = (-1.0)*A-(1.0)*B;
		gRPd = (-1.0)*A-(-3.0)*B;
		gRPv = 0.0;
		//gRPv =  (-1.0)*A-(5.0)*B; // with RH nu on
		gRPl = (-1.0)*A-(1.0)*B;
		
		// set fermion masses
		mtop = 171.0;

		// electoweak parameters
		sw2 = 0.2312;             // sin^2 theta_W       
		alpha = 1.0/137.036;
		//alpha = alphaEM(pow(_mZP,2));

		// derived quantities
		sw  = sqrt(sw2);
		cw2 = 1-sw2;
		cw  = sqrt(1-sw2);	
		
		// set Z' coupling constants
		e   = sqrt(4*pi*alpha);
		gZ  = e/(sw*cw);
		gZP = e*sqrt(5.0/(3.0*(1-sw2)));
		
		// set Z' mass in GeV
		mZP = _mZP;
		
		// set Z' decay width in GeV
		gammaZP = ZPrimeWidth(_mZP);
		
		// set a/Z/Z' interference
		gmZmode = _gmZmode;
		
	}
	
	// calculate the Z' decay width
	double ZPrimeWidth(double mZP) {
  
  		// decay width prefactor
		double gamma0 = pow(gZP,2)*mZP/(24.0*pi);

		// color factor to account for three quark colors
		double CF = 3;

		// effect of QCD radiative corrections
		//double rQCD = 1.0 + alphaS(mZP*mZP)/TMath::Pi();
		double rQCD  = 1.0;
  
		// compute Z' partial widths
		double gamma_v = gamma0*(pow(gLPv,2) + pow(gRPv,2));          // Z' partial width to neutrinos
		double gamma_l = gamma0*(pow(gLPl,2) + pow(gRPl,2));          // Z' partial width to charged leptons
		double gamma_d = gamma0*CF*rQCD*(pow(gLPd,2) + pow(gRPd,2));  // Z' partial width to down-quarks
		double gamma_u = gamma0*CF*rQCD*(pow(gLPu,2) + pow(gRPu,2));  // Z' partial width to light up-quarks
		//cout << "     Z' width to v: " << gamma_v << endl;
		//cout << "     Z' width to l: " << gamma_l << endl;
		//cout << "     Z' width to d: " << gamma_d << endl;
		//cout << "     Z' width to u: " << gamma_u << endl;

		// include threshold factors for top decays (top is an up-type quark)
		// in the small r limit, the top partial width becomes that of the up-quark
		double r = pow(mtop/mZP,2);
		double gamma_t = gamma0*CF*rQCD*(0.5*pow(gLPu+gRPu,2)*(1+2.0*r) + 0.5*pow(gLPu-gRPu,2)*(1-4.0*r))*sqrt(1-4.0*r);
		//cout << "     Z' width to t: " << gamma_t << endl;
  
		// compute Z' total width by summing over fermion partial widths (3 generations per fermion)
		gammaZP = 3.0*gamma_v + 3.0*gamma_l + 3.0*gamma_d + 2.0*gamma_u;
		if(mZP > 2.0*mtop) gammaZP += gamma_t;
		//cout << "     Z' Width: " << gammaZP << endl;
	
		return gammaZP;
	}

	// calculate running of alphaEM
	double alphaEM(const double &Q2) {
		double alphaEMmZ = 0.00781751;  // Pythia aEM value
		double b         = 0.725;       // Pythia mystery value
		double mZ  = 91.1876;           // GeV, Pythia value
		return alphaEMmZ/(1.0-b*alphaEMmZ*log(Q2/(mZ*mZ)));
	}
	
	// mathematical constants
	double pi;
	
	// fermion left-handed couplings 	
	double gLPu;
	double gLPd;
	double gLPv;
	double gLPl;

	// fermion right-handed couplings 		
	double gRPu;
	double gRPd;
	double gRPv;
	double gRPl;
	
	// fermion masses
	double mtop;
	
	// EW parameters
	double sw2;
	double alpha;
	
	// dervied 
	double sw;
	double cw2;
	double cw;
	
	// Z' parameters
	double e;
	double gZ;
	double gZP;
	double mZP;
	double gammaZP;
	
	// how Z' interfaces with Standard Model photon and Z boson
	int gmZmode;
	
};

double sigma(const double &x, const double &s, const ZPrimeModel &model) {
	
	// cross section parameteres
	// 4 pi alpha^2 (ħc)^2 / 3 s =  86.8 nB*GeV / s [GeV], where we used (ħc/GeV)^2 = 0.3894 mb = 0.3894E6 nB
	double hbarc2 = 0.3894*pow(10,6); // GeV nB  
	double pi = M_PI;

	// define electroweak scheme
	//double GF  = 1.16637*pow(10,-5); // GeV^-2
	double mZ  = 91.1876;            // GeV, Pythia value
	double sw2 = 0.2312;             // sin^2 theta_W       
	double alpha = 1.0/137.036; // uhhhhhhhhhhhhhhhhhh

	// set derived quantities from scheme choice
	double sw  = sqrt(sw2);
	double cw2 = 1-sw2;
	double cw  = sqrt(1-sw2);

	// ----------------------------------------

	// Drell-Yan Parameters
	
	// photon and Z gauge coupling constants
	double e   = sqrt(4*pi*alpha);
	double gZ  = e/(sw*cw);
	
	// Z boson width in GeV
	double gammaZ = 2.4952;
	
	// fermion electric charge Q
	int Qe = -1;
	int Qmu = -1;

	// fermion weak isospin t3
	double t3Le = -1./2.;
	double t3Lmu = -1./2.;

	// fermion left chiral charge g_{L}
	double gLe  = t3Le  - Qe*sw2;
	double gLmu = t3Lmu - Qmu*sw2;

	// fermion right chiral charge g_{R}
	double gRe  = 0     - Qe*sw2;
	double gRmu = 0     - Qmu*sw2;

	// ----------------------------------------
	
	// Z' Parameters
	
	// Z' gauge coupling constant
	double gZP = model.gZP;
	
	// Z' mass in GeV
	double mZP = model.mZP;
	
	// Z' decay width in GeV
	double gammaZP = model.gammaZP;
	
	// fermion left prime chiral charge g'_{L}
	double gLPe  = model.gLPl;
	double gLPmu = model.gLPl;

	// fermion right prime chiral charge g'_{R}
	double gRPe  = model.gRPl;
	double gRPmu = model.gRPl;

	// ----------------------------------------

	// cross section calculation
	
	// compute the mandalstam invariants
	double t = s/2.0*(1.0+x);
	double u = s/2.0*(1.0-x);

	// ...and square them for later use
	double s2 = pow(s,2);
	double t2 = pow(t,2);
	double u2 = pow(u,2);

	// initialize the photon propagator
	complex<double> aPropagator(s, 0.0);
	aPropagator = 1.0/aPropagator;

	// initialize the Z boson propagator
	//complex<double> zPropagator(s-mZ*mZ, s*gammaZ/mZ);
	complex<double> zPropagator(s-mZ*mZ, gammaZ*mZ);
	zPropagator = 1.0/zPropagator;
	
	// initialize Z' propagator
	//complex<double> zPrimePropagator(s-mZP*mZP, s*gammaZP/mZP);
	complex<double> zPrimePropagator(s-mZP*mZP, gammaZP*mZP);
	zPrimePropagator = 1.0/zPrimePropagator;
  
	// construct chiral amplitudes
	
	complex<double> ameLL  = (e*Qe)*(e*Qmu)*aPropagator;
	complex<double> zmeLL  = (gZ*gLe)*(gZ*gLmu)*zPropagator;
	complex<double> zpmeLL = (gZP*gLPe)*(gZP*gLPmu)*zPrimePropagator;
	
	complex<double> ameRR  = (e*Qe)*(e*Qmu)*aPropagator;
	complex<double> zmeRR  = (gZ*gRe)*(gZ*gRmu)*zPropagator;
	complex<double> zpmeRR = (gZP*gRPe)*(gZP*gRPmu)*zPrimePropagator;
	
	complex<double> ameLR  = (e*Qe)*(e*Qmu)*aPropagator;
	complex<double> zmeLR  = (gZ*gLe)*(gZ*gRmu)*zPropagator;
	complex<double> zpmeLR = (gZP*gLPe)*(gZP*gRPmu)*zPrimePropagator;
	
	complex<double> ameRL  = (e*Qe)*(e*Qmu)*aPropagator;
	complex<double> zmeRL  = (gZ*gRe)*(gZ*gLmu)*zPropagator;
	complex<double> zpmeRL = (gZP*gRPe)*(gZP*gLPmu)*zPrimePropagator;
	
	// construct total amplitude based on user set gmZmode
	int gmZmode = model.gmZmode;
	
	complex<double> meLL, meLR, meRL, meRR;
	if(gmZmode==0) {
		meLL = ameLL + zmeLL + zpmeLL;
		meLR = ameLR + zmeLR + zpmeLR;
		meRL = ameRL + zmeRL + zpmeRL;
		meRR = ameRR + zmeRR + zpmeRR;
	} else if(gmZmode==1) {
		meLL = ameLL;
		meLR = ameLR;
		meRL = ameRL;
		meRR = ameRR;
	} else if(gmZmode==2) {
		meLL = zmeLL;
		meLR = zmeLR;
		meRL = zmeRL;
		meRR = zmeRR;
	} else if(gmZmode==3) {
		meLL = zpmeLL;
		meLR = zpmeLR;
		meRL = zpmeRL;
		meRR = zpmeRR;
	} else if(gmZmode==4) {
		meLL = ameLL + zmeLL;
		meLR = ameLR + zmeLR;
		meRL = ameRL + zmeRL;
		meRR = ameRR + zmeRR;
	} else if(gmZmode==5) {
		meLL = ameLL + zpmeLL;
		meLR = ameLR + zpmeLR;
		meRL = ameRL + zpmeRL;
		meRR = ameRR + zpmeRR;
	} else if(gmZmode==6) {
		meLL = zmeLL + zpmeLL;
		meLR = zmeLR + zpmeLR;
		meRL = zmeRL + zpmeRL;
		meRR = zmeRR + zpmeRR;
	} else {
		cout << "failure" << endl;
	}
	
	
	// differential cross section prefactor: dsigma/dcos = 1/(32*pi*s)*|M|^2
	double sigma0 = 1.0/(32*pi*s);
	
	// compute the differential cross section dsigma/d\cos\theta
	double sigma = sigma0 * u2 * real(meLL*conj(meLL));
	sigma += sigma0 * u2 * real(meRR*conj(meRR));
	sigma += sigma0 * t2 * real(meLR*conj(meLR));
	sigma += sigma0 * t2 * real(meRL*conj(meRL));

	// finally give sigma dimensions
	sigma *= hbarc2;

	// cross check that the calculated cross section at s=mZ is consistent with the following
	/*
	double g0 = pow(gZ,2)*mZ/(24.0*pi);
	double gamma_e  = g0*(0.5*pow(gLe+gRe,2)*1 + 0.5*pow(gLe-gRe,2)*1);   // Z partial width to electrons
	double gamma_mu = g0*(0.5*pow(gLmu+gRmu,2)*1 + 0.5*pow(gLmu-gRmu,2)*1);  // Z partial width to muons
	double sigmaZ = 12*pi/pow(mZ,2)*gamma_e*gamma_mu/pow(gammaZ,2);
	sigmaZ *= hbarc2;
	cout << sigmaZ << endl;
	*/
	
	
	// cross check that the calculated cross section at s=mZ' is consistent with the following
	/*
	double gamma0 = pow(gZP,2)*mZP/(24.0*pi);
	double gamma_e = gamma0*(pow(gLPe,2) + pow(gRPe,2));          // Z' partial width to electrons
	double gamma_mu = gamma0*(pow(gLPmu,2) + pow(gRPmu,2));          // Z' partial width to muons
	double sigmaZP = 12*pi/pow(mZP,2)*gamma_e*gamma_mu/pow(gammaZP,2);
	sigmaZP *= hbarc2;
	cout << sigmaZP << endl;
	*/
	
	return sigma;
	
}

int main(int argc, char* argv[]) {
	
  	// seed random number generator
	int seed = 0x123;
	init_genrand64(seed);
	
	// initialize beam parameteres
	double sqrts = atof(argv[1]); // collider center mass (CM) energy in GeV
	double s = pow(sqrts,2);
	double E = sqrts/2; // energy of beams 1&2 in the CM frame
	
	// set the initial state four vectors in the CM frame
	double e_E=0, e_Px=0, e_Py=0, e_Pz=0;
	double p_E=0, p_Px=0, p_Py=0, p_Pz=0;
	
	double mu1_E=0, mu1_Px=0, mu1_Py=0, mu1_Pz=0;
	double mu2_E=0, mu2_Px=0, mu2_Py=0, mu2_Pz=0;
	
	// initialize Z'chi model with mZ' = 3 TeV and full interference: gamma/Z/Z'
	ZPrimeModel model("Chi",3000,0);
	
	// number of events to generate
	int nEvents = 1e5;
	
	// integration parameters
	double bounds = (1.0-(-1.0)); // bounds of integration for dcost (not integrating over phi)
	
	// cross section to compute with MC integration
	double f_avg = 0;
	double f2_avg = 0;
	double mcError = 0;
	
	
	
	
	
	//cout << "mu1_E, mu1_Px, mu1_Py, mu1_Pz, mu2_E, mu2_Px, mu2_Py, mu2_Pz, weight" << endl;
	cout << "costheta, weight" << endl;
	
	// do monte-carlo integration
	for(int i=1; i<nEvents; i++) {
		
		//if(i%1000==0) cout << "iteration " << i << endl;
		
		// random number generation
		double r1 = genrand64_real1();
		double r2 = genrand64_real2();
		//cout << r1 << " " << r2 << endl;

		// generate \cos(\theta) in interval of [-1,1]
		double cost = 1-2*r1;
		double sint = sqrt(1-pow(cost,2));
		
		// generate \phi in interval of [0,2pi)
		double phi = 2*pi*r2;
			
		// compute the cross section
		double sig = sigma(cost, s, model);
		
		// integrating dsigma/dcos over dcos
		f_avg += sig;
		f2_avg += pow(sig,2);
		
		// compute the event weight
		double wt = sig*bounds;
		
		//cout << i << " " << f_avg << " " << f2_avg << " " << wt << endl;
	
		// set initial state four vectors in the CM frame
		double e_E = E, e_Px = 0, e_Py = 0, e_Pz = E;
		double p_E = E, p_Px = 0, p_Py = 0, p_Pz = -E;
	
		// set final state four vectors in the CM frame
		double mu1_Px = E*sint*sin(phi);
		double mu1_Py = E*sint*cos(phi);
		double mu1_Pz = E*cost;
		double mu1_E  = E;
		
		double mu2_Px = -mu1_Px;
		double mu2_Py = -mu1_Py;
		double mu2_Pz = -mu1_Pz;
		double mu2_E  = E;
		
		// compute the running mc error
		mcError = bounds*sqrt( (f2_avg/i - pow(f_avg/i,2) ) / i);
		//cout << i << " " << mcError << endl;

		// write event kinematics to file
		//cout << mu1_E  << ", " << mu1_Px << ", " << mu1_Py << ", " << mu1_Pz << ", " << 
		//	mu2_E  << ", " << mu2_Px << ", " << mu2_Py << ", " << mu2_Pz << ", " << wt << endl;
		
		cout << cost << " " << wt << endl;
		
	}
	
	
	// compute final values
	
	f_avg = f_avg/nEvents;
	cout << "f_avg:     " << f_avg << endl;
	
	f2_avg=f2_avg/nEvents;
	cout << "f2_avg:    " << f2_avg << endl;
	
	double f = bounds*f_avg;
	
	mcError = bounds*sqrt( (f2_avg - pow(f_avg,2) ) / nEvents );
	cout << "cross section: " << f << " +/- " << mcError << " nb" << endl;
  
	return 0;
	
}

