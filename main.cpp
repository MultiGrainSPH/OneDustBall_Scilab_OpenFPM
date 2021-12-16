
#include <math.h>
#include <fstream>
#include "Vector/vector_dist.hpp"
#include "Draw/DrawParticles.hpp"


#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"
#include "NN/CellList/CellListM.hpp"
#include "Vector/vector_dist_multiphase_functions.hpp"

//#include "Grid/grid_dist_id.hpp"
//#include <map_vector.hpp>


// A constant to indicate boundary particles
#define BOUNDARY 0

// A constant to indicate Gas particles
#define Gas 1

// A constant to indicate Dust particles
#define Dust 2

#define X  0
#define Y  1
#define Z  2

// gamma in the formulas
double gamma_ = 4.0/3.0;

// alpha in the formula
const double viscoalpha = -1.0; // was -1.0
const double viscobeta  = 0.0; // was 0.0

double sizeX;// = 1.0; // [cm]

size_t countX;// = 4096;
size_t removed = 0;

double H;// = 1.0e-2;//sqrt(3.0)*dp;

double tstop;		// =rhod/tstop
double eps0; 		// rho_gas/rho_dust
double rho0; 		// start value of medium density
double E0; 			// start value of specific internal energy of gas
double u0gas, u0dust;// start value of the velocity on the edge of the ball

// Eta in the formulas
double Eta2;// = 0.01 * H*H;

double t_end = 0.001;



const double RhoMin = 0;
const double RhoMax = 1e+10;
const double Energymin = 1.0e-20;
double DtMin = 1e-9;
double DtMax = 5e-3;


// Constant used to define time integration
const double CFLnumber = 0.2; // was 0.4

double MassGas, MassDust;


// Properties
const size_t type = 0; 	// Gas, Dust or BOUNDARY
const int rho = 1; 		// Density
const int Pressure = 2; // Pressure
const int Energy = 3; 	// Internal energy
const int newrho = 4; 	// new rho calculated in the force calculation
const int force = 5; 	// calculated force
const int velocity = 6; // velocity
const int Mass = 7; 	// Mass of the SPH-particle
const int DeltaE = 8; 	// changing energy
const int v_prev = 9; // velocity

openfpm::vector<std::string> namesgas( {"type",  "rho",  "pressure","energy","newrho",  "force",  "V_gas", "Mass","de","v_prev"});
openfpm::vector<std::string> namesdust({"type_d","rho_d","2",       "E_d",   "newrho_d","force_d","U_dust","Mass","8", "u_prev"});
//openfpm::vector<std::string> namesanalytic( {"0","arho",  "apres","aenergy","4",  "5",  "aV_gas", "7","8","9"});

// Type of the vector containing particles
typedef vector_dist<3, double, aggregate<size_t,double,  double,    double,  double,  double[3], double[3],  double, double, double[3]>> particles;
//                  |   |                 |      |          |        |        |        |           |           |       |        |
//                  L___|                 |      |          |        |        |        |           |          mass   DeltaE    v_prev
//                  x,y,z                type   density  Pressure  energy   newrho    force     velocity
//

double Rgas = 1.0; 	// Radius of the gaseous ball
double Rdust = 1.0; // Radius of dusty ball
double zeta = 1.0; 	// factor of the individual force: vnew = v_old + (zeta*force_ind + (1 - zeta)*force_avg)*dt

openfpm::vector< particles > phases;

inline void EqState(particles & vd)
{
	auto it = vd.getDomainIterator();
	while (it.isNext())
	{
		auto a = it.get();
		const double rho_a = vd.template getProp<rho>(a);
        const double e_a = vd.template getProp<Energy>(a);

		vd.template getProp<Pressure>(a) = rho_a*(gamma_ - 1.0)*e_a;
		++it;
	}
}


double a2; //a2 = 1.0/M_PI/H/H/H; // 3D

inline double Wab(double r) // SPH-kernel
{
	r /= H; // doing q

	if (r < 1.0)
		return a2*(1.0 - 1.5*r*r + 0.75*r*r*r);
	else if (r < 2.0)
		return a2*(0.25*(2.0 - r)*(2.0 - r)*(2.0 - r));
	else
		return 0.0;
}


inline void DWab(Point<3, double> & dx, Point<3, double> & DW, double r)
{

    const double qq=r/H;
    const double fac1 = (-3.0 + 2.25*qq)/(H * H);
    const double b1 = (qq < 1.0)?1.0e0:0.0e0;
    const double wqq = (2.0 - qq);
    const double fac2 = -3.0 * wqq * wqq / (4.0 * qq * H * H);
    const double b2 = (qq >= 1.0 && qq < 2.0)?1.0e0:0.0e0;
    const double factor = (b1*fac1 + b2*fac2);

    DW.get(X) = a2 * factor * dx.get(X);
    DW.get(Y) = a2 * factor * dx.get(Y);
    DW.get(Z) = a2 * factor * dx.get(Z);
}




inline double Pi(const Point<3,double> & dr, double rr2, const Point<3,double> & dv, const double rhoa, const double rhob, const double Pa, const double Pb, double & visc)
{
	const double dot = dr.get(0)*dv.get(0) + dr.get(1)*dv.get(1) + dr.get(2)*dv.get(2);
	const double dot_rr2 = dot/(rr2 + Eta2);
	visc = std::max(dot_rr2, visc);

	if (dot < 0)
	{
		const double ca2 = Pa > 0.0 ? gamma_*Pa/rhoa : -gamma_*Pa/rhoa;
    	const double cb2 = Pb > 0.0 ? gamma_*Pb/rhob : -gamma_*Pb/rhob;
    	const double cbar = 0.5e0*(sqrt(ca2) + sqrt(cb2));
		const double amubar = H*dot_rr2; // was float
		const double rhobar = (rhoa + rhob)*0.5e0;
		const double pi_visc = (-viscoalpha*cbar*amubar + viscobeta*amubar*amubar)/rhobar;

		if (pi_visc != pi_visc){
		std::cout << "Pi_visc: " << ca2 << " "<< cb2 <<" " << amubar << " " << Pa <<" " << Pb << " " << rhoa << " " << rhob << std::endl;
		}

		return pi_visc;
    }
	else
		return 0.0;
}



void max_acceleration_and_velocity(particles & vd, double & max_acc, double & max_vel)
{
    // Calculate the maximum acceleration
    auto part = vd.getDomainIterator();
    double Pmin = 1000.0;
    while (part.isNext())
    {
        auto a = part.get();
        Point<3,double> acc(vd.getProp<force>(a));
        double acc2 = norm(acc);
        Point<3,double> vel(vd.getProp<velocity>(a));
        double vel2 = norm(vel);
        if (vel2 >= max_vel)
            max_vel = vel2;
        if (acc2 >= max_acc)
            max_acc = acc2;
        //const double Pa = vd.getProp<Pressure>(a);
		//if (Pmin > abs(Pa)) {
		//	Pmin = abs(Pa);
		//	Point<3,double> xa = vd.getPos(a);
		//	Rball = norm(xa);
		//}

        ++part;
    }
    //max_acc = sqrt(max_acc);
    //max_vel = sqrt(max_vel);
    Vcluster<> & v_cl = create_vcluster();
    v_cl.max(max_acc);
    v_cl.max(max_vel);
    //v_cl.min(Rball);
    v_cl.execute();
}


void getBall_radii(particles & gas, particles & dust, double & Rgas, double & Rdust)
{
    // Calculate the maximum acceleration
    auto gaspart = gas.getDomainIterator();
    double Pmin = 1000.0;
    while (gaspart.isNext())
    {
        auto a = gaspart.get();
        const double Pa = abs(gas.getProp<Pressure>(a));
		if (Pmin > Pa)	{
			Pmin = Pa;
			Point<3,double> xa = gas.getPos(a);
			Rgas = norm(xa);
		}
        ++gaspart;
    }

    auto dustpart = dust.getDomainIterator();
    while (dustpart.isNext())
    {
        auto a = dustpart.get();
		if (dust.template getProp<type>(a) == Dust) {
			Point<3,double> xa = dust.getPos(a);
			const double r = norm(xa);
			if ( r > Rdust) {
				Rdust = r;
			}
		}
        ++dustpart;
    }

    Vcluster<> & v_cl = create_vcluster();
    v_cl.min(Rgas);
    v_cl.max(Rdust);
    v_cl.execute();
}

double calc_deltaT(particles & vd, double ViscDtMax)
{
    double Maxacc = 0.0;
    double Maxvel = 0.0;
    max_acceleration_and_velocity(vd,Maxacc,Maxvel);
    //-dt1 depends on force per unit mass.
    const double dt_f = (Maxacc)?sqrt(H/Maxacc):std::numeric_limits<int>::max();
    //-dt2 combines the Courant and the viscous time-step controls.
    const double dt_cv = H/(Maxvel*5.0 + H*ViscDtMax);
    //-dt new value of time step.
    double dt = double(CFLnumber)*std::min(dt_f,dt_cv);
    if (dt < double(DtMin))
        dt = double(DtMin);
    if (dt > double(DtMax))
    	dt = double(DtMax);
    return dt;
}

template<typename CellList> inline void calc_IDIC(particles & gas, particles & dust,  CellList & NNGas, CellList & NNDust, CellList & NNIDIC, double & max_visc, double & dt) //CellsData & Cells, particles Cells,
{
	gas.updateCellList(NNGas); // Update the cell-lists
	dust.updateCellList(NNDust);
    gas.updateCellList(NNIDIC);

	auto CellsCount = NNIDIC.size();

	openfpm::vector <size_t> GasInCell(CellsCount), DustInCell(CellsCount);
	openfpm::vector <double> eps(CellsCount);
	openfpm::vector <double[3]> Vavg(CellsCount), Uavg(CellsCount), Aavg(CellsCount), PP(CellsCount), QQ(CellsCount);

    for (size_t i = 0; i < CellsCount; ++i){
        GasInCell.get(i) = 0;
        DustInCell.get(i) = 0;
      	eps.get(i) = 0.0;

        for (size_t q = X; q <= Z; ++q){
        	Vavg.get(i)[q] = 0.0;
        	Aavg.get(i)[q] = 0.0;
        	Uavg.get(i)[q] = 0.0;
			PP.get(i)[q] = 0.0;
        	QQ.get(i)[q] = 0.0;
        }

    }

    auto gaspart = gas.getDomainIterator();


	while (gaspart.isNext()) // For each particle ...
	{
		auto a = gaspart.get(); // ... a
		//if (gas.template getProp<type>(a) == Gas){
		Point<3,double> xa 	= gas.getPos(a); // Get the position xp of the particle
		const double massa 	= gas.template getProp<Mass>(a); // Take the mass of the particle dependently if it is Gas or BOUNDARY
		const double rhoa 	= gas.template getProp<rho>(a); // Get the density of the of the particle a
		const double Pa 	= gas.template getProp<Pressure>(a); // Get the pressure of the particle a
		Point<3,double> va 	= gas.template getProp<velocity>(a); // Get the Velocity of the particle a

		gas.template getProp<force>(a)[X] 	= 0.0; // Reset the force counter
		gas.template getProp<force>(a)[Y] 	= 0.0; // Reset the force counter
		gas.template getProp<force>(a)[Z] 	= 0.0; // Reset the force counter
		gas.template getProp<DeltaE>(a) = 0.0;

        auto Np = NNGas.template getNNIterator<NO_CHECK>(NNGas.getCell(xa)); // the neighborhood iterator

        // If it is a Gas particle calculate based on equation 1 and 2
		// Get an iterator over the neighborhood particles of a
		while (Np.isNext() == true) // For each neighborhood particle
		{
			auto b = Np.get(); // ... b
			if (a.getKey() == b) {++Np; continue; } // if (a == b) skip this particle
			Point<3,double> xb 	= gas.getPos(b); // Get the position xb of the particle
			const double massb 	= gas.template getProp<Mass>(b);
			Point<3,double> vb 	= gas.template getProp<velocity>(b);
			const double Pb 	= gas.template getProp<Pressure>(b);
			const double rhob 	= gas.template getProp<rho>(b);

			Point<3,double> dr = xa - xb; // Get the distance between a and b
			const double r = norm(dr); // take the norm of this vector, was norm2


			if (r < 2.0*H) // if they interact
			{

				const double r2 = r*r;//sqrt(r2);
				const Point<3, double> v_rel = va - vb;

    			Point<3, double> DW;
    			DWab(dr, DW, r);

    			const double factor = massb*( (Pa/(rhoa*rhoa) + Pb/(rhob*rhob)) + Pi(dr, r2, v_rel, rhoa, rhob, Pa, Pb, max_visc));

    			gas.template getProp<force>(a)[X] += -factor * DW.get(X);
    			gas.template getProp<force>(a)[Y] += -factor * DW.get(Y);
    			gas.template getProp<force>(a)[Z] += -factor * DW.get(Z);

    			const double factorE = massb* (Pa/(rhoa*rhoa)  +  0.5e0 * Pi(dr, r2, v_rel, rhoa, rhob, Pa, Pb, max_visc)) ; // + Pb/(rhob*rhob))
                gas.template getProp<DeltaE>(a) += factorE * (v_rel.get(X) * DW.get(X) + v_rel.get(Y) * DW.get(Y) + v_rel.get(Z) * DW.get(Z));


			}
			++Np;
		}


		//if (gas.template getProp<type>(a) = Gas){
			auto CellID = NNIDIC.getCell(xa);
			const Point<3,double> fa = gas.template getProp<force>(a);

        	++GasInCell.get(CellID);
        	Vavg.get(CellID)[X] += va.get(X);
        	Vavg.get(CellID)[Y] += va.get(Y);
        	Vavg.get(CellID)[Z] += va.get(Z);

        	Aavg.get(CellID)[X] += fa.get(X);
        	Aavg.get(CellID)[Y] += fa.get(Y);
        	Aavg.get(CellID)[Z] += fa.get(Z);
       // }

		//}
		++gaspart;
	}


    auto dustpart = dust.getDomainIterator();

	while (dustpart.isNext()) // For each dust particle ...
	{
		auto a = dustpart.get(); // ... a
		//if (dust.template getProp<type>(a) == Dust){
			Point<3,double> xa = dust.getPos(a); // Get the position xa of the particle
			const Point<3,double> va = dust.template getProp<velocity>(a);
			//const double massa = dust.template getProp<Mass>(a); // Take the mass of the particle

			dust.template getProp<force>(a)[X]  = 0.0; // Reset the force counter
			dust.template getProp<force>(a)[Y]  = 0.0; // Reset the force counter
			dust.template getProp<force>(a)[Z]  = 0.0; // Reset the force counter


			//if (dust.template getProp<type>(a) = Dust){
				auto CellID = NNIDIC.getCell(xa);
				++DustInCell.get(CellID);
	        	Uavg.get(CellID)[X] += va.get(X);
	        	Uavg.get(CellID)[Y] += va.get(Y);
	        	Uavg.get(CellID)[Z] += va.get(Z);
	        //}

			//}
		++dustpart;
	}

	// Calculate delta t integration
	//dt = calc_deltaT(gas, max_visc);
	dt = DtMin;

    // calculate parameters for the next sub-step with known average parameters

    for (size_t i = 0; i < CellsCount; ++i){

    	const double Ngas = double(GasInCell.get(i));
    	const double Ndust = double(DustInCell.get(i));
    	if (Ngas > 0) {
    		eps.get(i) = (Ndust*MassDust)/(Ngas*MassGas);
    	} else {
    		eps.get(i) = 1e200;
    	}

		for (size_t q = X; q <= Z; ++q){

    		if (Ngas > 0) {
    			Vavg.get(i)[q] /= Ngas;
    			Aavg.get(i)[q] /= Ngas;
    		} else {
    			Vavg.get(i)[q] = 0.0;
    			Aavg.get(i)[q] = 0.0;
    		}

    		if (Ndust > 0) {
    			Uavg.get(i)[q] /= Ndust;
			} else {
    			Uavg.get(i)[q] = 0.0;
    		}

    	} // q

    	for (size_t q = X; q <= Z; ++q){
    		PP.get(i)[q] = Vavg.get(i)[q] - Uavg.get(i)[q];
    		QQ.get(i)[q] = Vavg.get(i)[q] + eps.get(i)*Uavg.get(i)[q];
		}

		for (size_t q = X; q <= Z; ++q){
    		PP.get(i)[q] = (PP.get(i)[q] + Aavg.get(i)[q]*dt )/(1.0 + (eps.get(i) + 1.0)*dt/tstop);
    		QQ.get(i)[q] = QQ.get(i)[q] + Aavg.get(i)[q]*dt;

    		Vavg.get(i)[q] = (QQ.get(i)[q] + eps.get(i)*PP.get(i)[q])/( 1.0 + eps.get(i) );
    		Uavg.get(i)[q] = (QQ.get(i)[q] - PP.get(i)[q])/( 1.0 + eps.get(i) );
    	}
    }



    // calculate new velocities
    gaspart = gas.getDomainIterator();
    while (gaspart.isNext()){
		auto a = gaspart.get();
       	auto CellID = NNIDIC.getCell(gas.getPos(a));

       	if (gas.template getProp<type>(a) == Gas){ // not BOUNDARY!
			//gas.template getProp<velocity>(a)[X] = gas.template getProp<v_prev>(a)[X];
			//gas.template getProp<velocity>(a)[Y] = gas.template getProp<v_prev>(a)[Y];
			//gas.template getProp<velocity>(a)[Z] = gas.template getProp<v_prev>(a)[Z];
    	//} else {*/

		for (size_t q = X; q <= Z; ++q){

			gas.template getProp<v_prev>(a)[q] = gas.template getProp<velocity>(a)[q];
       		//gas.template getProp<force>(a)[q] += eps.get(CellID)*Uavg.get(CellID)[q]/tstop;
       		if (eps.get(CellID) != 0.0 && eps.get(CellID) < 1e200 ){
       	    	gas.template getProp<velocity>(a)[q] = (gas.template getProp<v_prev>(a)[q] + (zeta*gas.template getProp<force>(a)[q] + (1.0 - zeta)*Aavg.get(CellID)[q])*dt + eps.get(CellID)*Uavg.get(CellID)[q]*dt/tstop)/(1.0 + dt*eps.get(CellID)/tstop);
       	    	//factor of the individual force: vnew = v_old + (zeta*force_ind + (1 - zeta)force_avg)*dt

       	    	//gas.template getProp<velocity>(a)[q] = (gas.template getProp<v_prev>(a)[q] + Aavg.get(CellID)[q]*dt + eps.get(CellID)*Uavg.get(CellID)[q]*dt/tstop)/(1.0 + dt*eps.get(CellID)/tstop); // AVERAGE FORCE
       	    	//gas.template getProp<velocity>(a)[q] = (gas.template getProp<v_prev>(a)[q] + gas.template getProp<force>(a)[q]*dt + eps.get(CellID)*Uavg.get(CellID)[q]*dt/tstop)/(1.0 + dt*eps.get(CellID)/tstop); // INDIVIDUAL FORCE
       		}else{
       	    	gas.template getProp<velocity>(a)[q] = gas.template getProp<v_prev>(a)[q] + gas.template getProp<force>(a)[q]*dt;
       		}

       		if (gas.template getProp<velocity>(a)[q] != gas.template getProp<velocity>(a)[q] )	{
       			Point<3, double> xa = gas.getPos(a);
				std::cout << "q=" << q << " "<<gas.template getProp<v_prev>(a)[q] << " " << gas.template getProp<force>(a)[q]*dt << " " << eps.get(CellID)*Uavg.get(CellID)[q];
       			std::cout << " [" << xa.get(X) << "," << xa.get(Y) << "," << xa.get(Z) << "]" <<std::endl;
       		}


       	}
       	}//not boundary

		++gaspart;
	}


    dustpart = dust.getDomainIterator();
    while (dustpart.isNext()){
		auto a = dustpart.get();
       	auto CellID = NNIDIC.getCell(dust.getPos(a));

		if (dust.template getProp<type>(a) == Dust){
		for (size_t q = X; q <= Z; ++q){
			dust.template getProp<v_prev>(a)[q] = dust.template getProp<velocity>(a)[q];
        	dust.template getProp<force>(a)[q] =  Vavg.get(CellID)[q]/tstop;
        	//if (eps.get(CellID) != 0.0 && eps.get(CellID) < 1e200){
       	    	dust.template getProp<velocity>(a)[q] = (dust.template getProp<v_prev>(a)[q] + dust.template getProp<force>(a)[q]*dt)/(1.0 + dt/tstop);
       		//}
       	}
       	}// not boundary
		++dustpart;
	}


}

template<typename CellList> inline void calc_df(particles & gas, particles & dust,  CellList & NNGas, CellList & NNDust, CellList & NNIDIC, double & max_visc, double & dt) //CellsData & Cells, particles Cells,
{
	gas.updateCellList(NNGas); // Update the cell-lists
	dust.updateCellList(NNDust);

    auto gaspart = gas.getDomainIterator();

	double vrbound_Gas = 0.0;
	double sumr2_Gas = 0.0;

	while (gaspart.isNext()) // For each particle ...
	{
		auto a = gaspart.get(); // ... a
		Point<3,double> xa 	= gas.getPos(a); // Get the position xp of the particle
		const double massa 	= gas.template getProp<Mass>(a); // Take the mass of the particle dependently if it is Gas or BOUNDARY
		//const double rhoa 	= gas.template getProp<rho>(a); // Get the density of the of the particle a
		//const double Pa 	= gas.template getProp<Pressure>(a); // Get the pressure of the particle a
		//Point<3,double> va 	= gas.template getProp<velocity>(a); // Get the Velocity of the particle a
		Point<3,double> vanext 	= gas.template getProp<velocity>(a);

		gas.template getProp<newrho>(a) 	= massa*Wab(0.0);
        //gas.template getProp<DeltaE>(a)   	= 0.0;

		if (gas.template getProp<type>(a) == Gas){
    				vrbound_Gas += vanext.get(X)*xa.get(X) + vanext.get(Y)*xa.get(Y) + vanext.get(Z)*xa.get(Z);
    				sumr2_Gas += norm2(xa);
    			}


        auto Np = NNGas.template getNNIterator<NO_CHECK>(NNGas.getCell(xa)); // the neighborhood iterator

 		// Get an iterator over the neighborhood particles of a
		while (Np.isNext() == true) // For each neighborhood particle
		{
			auto b = Np.get(); // ... b
			if (a.getKey() == b) {++Np; continue; } // if (a == b) skip this particle
			Point<3,double> xb 	= gas.getPos(b); // Get the position xb of the particle
			const double massb 	= gas.template getProp<Mass>(b);
			//Point<3,double> vb 	= gas.template getProp<velocity>(b);
			//Point<3,double> vbnext 	= gas.template getProp<velocity>(b);
			//const double Pb 	= gas.template getProp<Pressure>(b);
			//const double rhob 	= gas.template getProp<rho>(b);

			Point<3,double> dr = xa - xb; // Get the distance between a and b
			const double r = norm(dr); // take the norm of this vector, was norm2


			if (r < 2.0*H) // if they interact
			{
				gas.template getProp<newrho>(a) += massb*Wab(r);
				//const double r2 = r*r;//sqrt(r2);
				//const Point<3, double> v_rel = va - vb;

    			//Point<3, double> DW;
    			//DWab(dr, DW, r);

    			//const double factor = massb* (Pa/(rhoa*rhoa)  +  0.5e0 * Pi(dr, r2, v_rel, rhoa, rhob, Pa, Pb, max_visc)) ; // + Pb/(rhob*rhob))

                //gas.template getProp<DeltaE>(a) += factor * (v_rel.get(X) * DW.get(X) + v_rel.get(Y) * DW.get(Y) + v_rel.get(Z) * DW.get(Z));
                //gas.template getProp<DeltaE>(a) += -massb*Pa/(rhoa*rhoa) * (vb.get(X) * DW.get(X) + vb.get(Y) * DW.get(Y) + vb.get(Z) * DW.get(Z));
                //gas.template getProp<DeltaE>(a) += -massb*Pb/(rhob*rhob) * (va.get(X) * DW.get(X) + va.get(Y) * DW.get(Y) + va.get(Z) * DW.get(Z));
                //gas.template getProp<DeltaE>(a) += massa*(Pa/(rhoa*rhoa) + 0.5*Pi(dr, r2, v_rel, rhoa, rhob, Pa, Pb, max_visc)) * (v_rel.get(X) * DW.get(X) + v_rel.get(Y) * DW.get(Y) + v_rel.get(Z) * DW.get(Z));

			}
			++Np;
		}

		++gaspart;
	}


    auto dustpart = dust.getDomainIterator();
    double vrbound_Dust = 0.0;
	double sumr2_Dust = 0.0;

	while (dustpart.isNext()) // For each dust particle ...
	{
		auto a = dustpart.get(); // ... a
		Point<3,double> xa = dust.getPos(a); // Get the position xa of the particle
		const double massa = dust.template getProp<Mass>(a); // Take the mass of the particle
		Point<3,double> vanext 	= dust.template getProp<velocity>(a);

		//dust.template getProp<force>(a)[X]  = 0.0; // Reset the force counter
		//dust.template getProp<force>(a)[Y]  = 0.0; // Reset the force counter
		//dust.template getProp<force>(a)[Z]  = 0.0; // Reset the force counter
		dust.template getProp<newrho>(a)   	= massa*Wab(0.0);
        //dust.template getProp<DeltaE>(a)   	= 0.0;

        if (dust.template getProp<type>(a) == Dust){
    				vrbound_Dust += vanext.get(X)*xa.get(X) + vanext.get(Y)*xa.get(Y) + vanext.get(Z)*xa.get(Z);
    				sumr2_Dust += norm2(xa);
    			}

        auto Np = NNDust.template getNNIterator<NO_CHECK>(NNDust.getCell(dust.getPos(a))); // the neighborhood iterator

        // If it is a Gas particle calculate based on equation 1 and 2
		// Get an iterator over the neighborhood particles of a
		while (Np.isNext() == true) // For each neighborhood particle
		{
			auto b = Np.get(); // ... b
			if (a.getKey() == b) { ++Np;	continue; } // if (a == b) skip this particle
			Point<3,double> xb = dust.getPos(b); // Get the position xp of the particle
			const double massb = dust.template getProp<Mass>(b);

			Point<3,double> dr = xa - xb; // Get the distance between a and b
			const double r = norm(dr); // take the norm of this vector, was norm2

			if (r < 2.0*H) // if they interact
		    {
        	    dust.template getProp<newrho>(a) += massb*Wab(r);
		    }
		    ++Np;
		}


		++dustpart;
	}


	// Get the values of velocities of the boundary particles
	gaspart = gas.getDomainIterator();
 	// Get an iterator over particles
	while (gaspart.isNext() == true) // For each boundary particle
	{
		auto a = gaspart.get();
		Point<3,double> xa 	= gas.getPos(a);
		if (gas.template getProp<type>(a) == BOUNDARY){

			gas.template getProp<velocity>(a)[X] = vrbound_Gas * xa.get(X)/sumr2_Gas;
			gas.template getProp<velocity>(a)[Y] = vrbound_Gas * xa.get(Y)/sumr2_Gas;
			gas.template getProp<velocity>(a)[Z] = vrbound_Gas * xa.get(Z)/sumr2_Gas;

			gas.template getProp<v_prev>(a)[X] = gas.template getProp<velocity>(a)[X];
			gas.template getProp<v_prev>(a)[Y] = gas.template getProp<velocity>(a)[Y];
			gas.template getProp<v_prev>(a)[Z] = gas.template getProp<velocity>(a)[Z];

    	}

    	++gaspart;
	}

	dustpart = dust.getDomainIterator();
	while (dustpart.isNext() == true) // For each boundary particle
	{
		auto a = dustpart.get();
		Point<3,double> xa 	= dust.getPos(a);
		if (dust.template getProp<type>(a) == BOUNDARY){

			dust.template getProp<velocity>(a)[X] = vrbound_Dust * xa.get(X)/sumr2_Dust;
			dust.template getProp<velocity>(a)[Y] = vrbound_Dust * xa.get(Y)/sumr2_Dust;
			dust.template getProp<velocity>(a)[Z] = vrbound_Dust * xa.get(Z)/sumr2_Dust;

			dust.template getProp<v_prev>(a)[X] = dust.template getProp<velocity>(a)[X];
			dust.template getProp<v_prev>(a)[Y] = dust.template getProp<velocity>(a)[Y];
			dust.template getProp<v_prev>(a)[Z] = dust.template getProp<velocity>(a)[Z];

    	}

    	++dustpart;
	}

}




openfpm::vector<size_t> to_remove;

size_t cnt = 0;

void euler_int_position(particles & vd, const double dt, const size_t TheType)
{
	// particle iterator
	auto part = vd.getDomainIterator();

	// For each particle ...
	while (part.isNext())
	{
		// ... a
		auto a = part.get();

		// if the particle is boundary

        	const double vel[3] = {vd.template getProp<v_prev>(a)[X],vd.template getProp<v_prev>(a)[Y], vd.template getProp<v_prev>(a)[Z]};

			//-Calculate displacement and update position
			const double dx = vel[X]*dt;// + acc*dt205;
			const double dy = vel[Y]*dt;
			const double dz = vel[Z]*dt;

	    	vd.getPos(a)[X] += dx;
	    	vd.getPos(a)[Y] += dy;
	    	vd.getPos(a)[Z] += dz; // Now the particles have a new position! Based on this position we should calculate new rho!

	    	Point<3,double> xa 	= vd.getPos(a);
	    	const double r = norm(xa);

	    	//if ((vd.template getProp<type>(a) == BOUNDARY) &&  (r <= Rball)){
	    	//	vd.template getProp<type>(a) = TheType; // if the particle goes out from the ball
	    	//}
	 		++part;
	}

}

void euler_int_rhoE(particles & vd, double dt)
{
	// list of the particle to remove
	to_remove.clear();

	// particle iterator
	auto part = vd.getDomainIterator();

	// For each particle ...
	while (part.isNext())
	{
		auto a = part.get();

        	vd.template getProp<Energy>(a) += vd.template getProp<DeltaE>(a)*dt;
        	vd.template getProp<rho>(a) = vd.template getProp<newrho>(a);

			if (vd.template getProp<Energy>(a) < 0.0 && vd.template getProp<type>(a) != BOUNDARY){
				vd.template getProp<Energy>(a) = Energymin;
			}

	    	// Check if the particle go out of range in space and in density vd.getPos(a)[X] != vd.getPos(a)[X] ||
	    	if ( vd.template getProp<rho>(a) < RhoMin ||
	    		 vd.template getProp<rho>(a) > RhoMax ||
	    		 vd.template getProp<rho>(a) != vd.template getProp<rho>(a) ||
	    		 vd.template getProp<Energy>(a) != vd.template getProp<Energy>(a) ||
	    		 //vd.template getProp<Energy>(a) < 0.0 ||
	    		 vd.getPos(a)[X] != vd.getPos(a)[X] ||
	    		 vd.getPos(a)[Y] != vd.getPos(a)[Y] ||
	    		 vd.getPos(a)[Z] != vd.getPos(a)[Z] )
	    	{
          		removed++;
	       	    to_remove.add(a.getKey());
	    	}

        ++part;
	}

	// remove the particles
	vd.remove(to_remove, 0);
}




int main(int argc, char* argv[])
{
    // initialize the library
	openfpm_init(&argc, &argv);

    std::ifstream fin("input.ini");
    //char buff[50];
    double dt;
    double AllMass ;// Mass of the gaseous ball //= 4.0*M_PI/3.0*sizeX*sizeX*sizeX*rho0;
    double kappa;
    fin >> sizeX >> Rdust >> countX >> H >> kappa >> t_end >> DtMin >> AllMass >> E0 >> tstop >> eps0 >> u0gas >> u0dust >> gamma_ >> zeta;
    fin.close();

    std::cout << " Rgas = sizeX = " << sizeX << " Rdust = " << Rdust << std::endl;
    std::cout << " CountX = " << countX << " H = " << H << " kappa = " << kappa << std::endl;
    std::cout << " t_end = " << t_end << " dt_min = " << DtMin << std::endl;
    std::cout << " E0 = " << E0 << " tstop = " << tstop << " eps0 = " << eps0 << " u0gas = " << u0gas << " u0dust = " << u0dust << std::endl;
    std::cout << " gamma = "<< gamma_ << " zeta = " << zeta << std::endl;

	Rgas = sizeX;

    dt = DtMin;
    const double dp = 2.0*sizeX/(double(countX) );
	rho0 = AllMass/(4.0*M_PI/3.0*sizeX*sizeX*sizeX);
    a2 = 1.0/(M_PI * H * H * H); // 3D   // a2 = 2.0/(3.0*H); // 1D
    std::cout << "GasMass = " << AllMass << " rho0 = " << rho0 << " a2 = " << a2 << " dp = "<< dp << std::endl;
    Eta2 = 0.01 * H * H;

	// Here we define our domain a 3D box
	const double mm = 1.8;
	Box<3,double> domain({-mm*sizeX, -mm*sizeX, -mm*sizeX}, {mm*sizeX, mm*sizeX, mm*sizeX}); // full coordinate domain where everything happens
	//Box<3,double> Box1({-sizeX+dp/2.0, -sizeX+dp/2.0, -sizeX+dp/2.0}, {sizeX-dp/2.0, sizeX-dp/2.0, sizeX-dp/2.0});
	Box<3,double> Box1({-sizeX-2.0*H, -sizeX-2.0*H, -sizeX-2.0*H}, {sizeX+2.0*H, sizeX+2.0*H, sizeX+2.0*H});
	//size_t sz1[3] = {countX, countX, countX}; // count of particles in any directions
	size_t sz2[3] = {int(mm*(countX+1)), int(mm*(countX+1)), int(mm*(countX+1))}; // count of particles in any directions

	// Here we define the boundary conditions of our problem
    size_t bc[3] = {PERIODIC, PERIODIC, PERIODIC};

	// extended boundary around the domain, and the processor domain
	Ghost<3,double> g(2.1 * H); // 3 -- means 3D

	phases.add(particles(0, domain, bc, g)); // add gas , DEC_GRAN(2)
	phases.add(particles(0, domain, bc, g)); // add dust this way , DEC_GRAN(2)
	//phases.add(particles(0, domain, bc, g)); // add particles for analytical solution

	// return an iterator to the Gas particles to add to vd
	auto Gas_it  = DrawParticles::DrawBox(phases.get(0), sz2, domain, Box1);
	auto Dust_it = DrawParticles::DrawBox(phases.get(1), sz2, domain, Box1);
	//auto AGas_it = DrawParticles::DrawBox(phases.get(2), sz2, domain, Box1);


	// for each particle inside the box ...
	while (Gas_it.isNext())
	{

		// ... set its position ...
		const double xx = Gas_it.get().get(X) + dp/2.0;
		const double yy = Gas_it.get().get(Y) + dp/2.0;
		const double zz = Gas_it.get().get(Z) + dp/2.0;
		const double rr = sqrt(xx*xx + yy*yy + zz*zz);

		if (rr <= Rgas + 2.0*H) {

			// ... add a particle ...
			phases.get(0).add();

			phases.get(0).getLastPos()[X] = xx;
			phases.get(0).getLastPos()[Y] = yy;
			phases.get(0).getLastPos()[Z] = zz;

			if (rr <= Rgas) {
			phases.get(0).template getLastProp<type>() = Gas;


			} else {
				phases.get(0).template getLastProp<type>() = BOUNDARY;

			}

			phases.get(0).template getLastProp<rho>() = rho0;

		size_t ii = 0;
		//for (size_t ii = 0; ii < 1; ++ii){ // Gas only
		// and parameters.
            phases.get(ii).template getLastProp<Mass>() = phases.get(ii).template getLastProp<rho>()*dp*dp*dp;
            const double EE = E0*(1.0 - (rr/sizeX)*(rr/sizeX));
            phases.get(ii).template getLastProp<Energy>() = EE;// > 0.0? EE : 0.0;
            phases.get(ii).template getLastProp<Pressure>() = phases.get(ii).template getLastProp<rho>()*(gamma_ - 1.0)*phases.get(ii).template getLastProp<Energy>(); // Eqstate
            phases.get(ii).template getLastProp<velocity>()[X] = u0gas*xx/Rgas;
            phases.get(ii).template getLastProp<velocity>()[Y] = u0gas*yy/Rgas;
            phases.get(ii).template getLastProp<velocity>()[Z] = u0gas*zz/Rgas;
            phases.get(ii).template getLastProp<v_prev>()[X] = u0gas*xx/Rgas;
            phases.get(ii).template getLastProp<v_prev>()[Y] = u0gas*yy/Rgas;
            phases.get(ii).template getLastProp<v_prev>()[Z] = u0gas*zz/Rgas;

		//}

		MassGas = phases.get(0).template getLastProp<Mass>();

		}

		if (rr <= Rdust + 2.0*H) {

			// ... add a particle ...
			phases.get(1).add();
			//phases.get(2).add();

			phases.get(1).getLastPos()[X] = xx + (1.0e-4)*(rand()/RAND_MAX - 0.5);
			phases.get(1).getLastPos()[Y] = yy + (1.0e-4)*(rand()/RAND_MAX - 0.5);
			phases.get(1).getLastPos()[Z] = zz + (1.0e-4)*(rand()/RAND_MAX - 0.5);


			if (rr <= Rdust) {
			phases.get(1).template getLastProp<type>() = Dust;

			} else {
				phases.get(1).template getLastProp<type>() = BOUNDARY;
			}

			phases.get(1).template getLastProp<rho>() = eps0*rho0; // dust to gas ratio

		size_t ii = 1;
		//for (size_t ii = 1; ii < 2; ++ii){ // Dust only
		// and parameters.
            phases.get(ii).template getLastProp<Mass>() = phases.get(ii).template getLastProp<rho>()*dp*dp*dp;
            phases.get(ii).template getLastProp<velocity>()[X] = u0dust*xx/Rdust;
            phases.get(ii).template getLastProp<velocity>()[Y] = u0dust*yy/Rdust;
            phases.get(ii).template getLastProp<velocity>()[Z] = u0dust*zz/Rdust;
            phases.get(ii).template getLastProp<v_prev>()[X] = u0dust*xx/Rdust;
            phases.get(ii).template getLastProp<v_prev>()[Y] = u0dust*yy/Rdust;
            phases.get(ii).template getLastProp<v_prev>()[Z] = u0dust*zz/Rdust;

		//}
		MassDust = phases.get(1).template getLastProp<Mass>();

		}

		// next Gas or Dust particle
		++Gas_it;
		++Dust_it;
		//++AGas_it;
	}

//std::cout <<"Mass of gas sph particle is " << MassGas << " calc: " << AllMass/(M_PI*countX*countX*countX/6.0) << std::endl;

	phases.get(0).map();
	phases.get(1).map();
	EqState(phases.get(0));

	phases.get(0).ghost_get<0,1,2,3,4,5,6,7,8,9>();//
	phases.get(1).ghost_get<0,1,2,3,4,5,6,7,8,9>();

//std::cout <<"Start NNGas" << std::endl;
	auto NNGas  = phases.get(0).getCellList(2.0*H);
//std::cout <<"Start NNDust" << std::endl;
	auto NNDust = phases.get(1).getCellList(2.0*H); // maybe, vector??
//std::cout <<"Start NNIDIC" << std::endl;
    auto NNIDIC = phases.get(0).getCellList(kappa*H);

	phases.get(0).setPropNames(namesgas);
	phases.get(1).setPropNames(namesdust);
	//phases.get(2).setPropNames(namesanalytic);

//std::cout <<"Names setted" << std::endl;
    //tstop = 1.0/K; // density of dust is 1.0 g/cm^3


	// Evolve

	size_t write = 0;
	size_t it = 0;
	//size_t it_reb = 0;
	double t = 0.0;

    double max_visc = 0.0;

    phases.get(0).write_frame("slices/Gas",write,"",CSV_WRITER); // write initial conditions
	phases.get(1).write_frame("slices/Dust",write,"",CSV_WRITER);
	//phases.get(2).write_frame("slices/AGas",write);
	write++;

	std::cout << "TIME: " << t << "  write initial conditions, max_visc = " << max_visc << std::endl;

	std::ofstream foutr("radius.tab");
	foutr << t << " " << Rgas << " " << Rdust << std::endl;

	// ======================= MAIN LOOP =======================
	max_visc = 0.0;
	while (t <= t_end)
	{

		Vcluster<> & v_cl = create_vcluster();
		//timer it_time;


		phases.get(0).map();
		phases.get(1).map();
		//phases.get(2).map();

		//EqState(phases.get(0));

		// Calculate pressure from the density

		max_visc /= 10.0;
        v_cl.max(max_visc);
        v_cl.execute();


        phases.get(0).ghost_get<0,1,2,3,4,5,6,7,8,9>();
		phases.get(1).ghost_get<0,1,2,3,4,5,6,7,8,9>();


		// Calc forces phases.get(2),
		calc_IDIC(phases.get(0), phases.get(1), NNGas, NNDust, NNIDIC, max_visc, dt); //phases.get(2)
		euler_int_position(phases.get(0), dt, Gas);
		euler_int_position(phases.get(1), dt, Dust);
		//euler_int_position(phases.get(2), dt);
		phases.get(0).map();
		phases.get(1).map();
		//phases.get(2).map();
		phases.get(0).ghost_get<0,1,2,3,4,5,6,7,8,9>();
		phases.get(1).ghost_get<0,1,2,3,4,5,6,7,8,9>();
		calc_df(phases.get(0), phases.get(1), NNGas, NNDust, NNIDIC, max_visc, dt);
		//phases.get(0).ghost_get<0,1,2,3,4,5,6,7,8,9>();
		//phases.get(1).ghost_get<0,1,2,3,4,5,6,7,8,9>();
		euler_int_rhoE(phases.get(0), dt);
		euler_int_rhoE(phases.get(1), dt);
		EqState(phases.get(0));
		getBall_radii(phases.get(0), phases.get(1), Rgas, Rdust);

		// Euler step
		it++;

		// Get the maximum viscosity term across processors
		v_cl.max(max_visc);
		v_cl.execute();


		cnt++;


		t += dt;

		if (write <= t*100)
		{
			phases.get(0).write_frame("slices/Gas",write,"",CSV_WRITER);
			phases.get(1).write_frame("slices/Dust",write,"",CSV_WRITER);


			if (v_cl.getProcessUnitID() == 0)
			{std::cout << "TIME: " << t << "  write file " << write << "   procID=" << v_cl.getProcessUnitID() << "   " << cnt << "   Max visc: " << max_visc << " dt = " << dt << " Removed: " << removed << " Rgas|dust = " << Rgas << "|" << Rdust << std::endl;
			foutr << t << " " << Rgas << " " << Rdust << std::endl;
			}
			write++;
		}
		/*else
		{
			if (v_cl.getProcessUnitID() == 0)
			{std::cout << "TIME: " << t << "        " << it_time.getwct() << "   procID=" << v_cl.getProcessUnitID() << "   " << cnt << "    Max visc: " << max_visc << " Removed: " << removed<< std::endl;}
		}*/

	}

	foutr.close();

	openfpm_finalize();
}
