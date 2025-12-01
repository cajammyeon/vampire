//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2022. All rights reserved.
//
//   Email: mara.strungaru@york.ac.uk
//
//   implementation based on the paper Phys. Rev. B 103, 024429, (2021) M.Strungaru, M.O.A. Ellis et al
//------------------------------------------------------------------------------


// C++ standard library headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

// Vampire headers
#include "anisotropy.hpp"
#include "atoms.hpp"
#include "create.hpp"
#include "sld.hpp"
#include "sim.hpp"

// sld module headers
#include "internal.hpp"

namespace sld
{

   	void compute_fields(const int start_index, // first atom for exchange interactions to be calculated
						const int end_index,
						const std::vector<int>& neighbour_list_start_index,
						const std::vector<int>& neighbour_list_end_index,
						const std::vector<int>& type_array, // type for atom
						const std::vector<int>& neighbour_list_array, // list of interactions between atom
						const std::vector<double>& x_coord_array, // coord vectors for atoms
						const std::vector<double>& y_coord_array,
						const std::vector<double>& z_coord_array,
						const std::vector<double>& x_spin_array, // coord vectors for atoms
						const std::vector<double>& y_spin_array,
						const std::vector<double>& z_spin_array,
						std::vector<double>& forces_array_x, //  vectors for forces
						std::vector<double>& forces_array_y,
						std::vector<double>& forces_array_z,
						std::vector<double>& fields_array_x, //  vectors for fields
						std::vector<double>& fields_array_y,
						std::vector<double>& fields_array_z)
	{


		internal::compute_exchange(start_index, end_index,
								neighbour_list_start_index, neighbour_list_end_index,
								type_array, neighbour_list_array,
								x_coord_array, y_coord_array, z_coord_array,
								x_spin_array, y_spin_array, z_spin_array,
								forces_array_x, forces_array_y, forces_array_z,
								fields_array_x, fields_array_y, fields_array_z);

		if(sld::internal::pseudodipolar) internal::compute_sld_coupling(start_index, end_index,
																		neighbour_list_start_index, neighbour_list_end_index,
																		type_array, neighbour_list_array,
																		x_coord_array, y_coord_array, z_coord_array,
																		x_spin_array, y_spin_array, z_spin_array,
																		forces_array_x, forces_array_y, forces_array_z,
																		fields_array_x, fields_array_y, fields_array_z);
		if(sld::internal::full_neel) internal::compute_sld_coupling_neel(start_index, end_index,
																		neighbour_list_start_index, neighbour_list_end_index,
																		type_array, neighbour_list_array,
																		x_coord_array, y_coord_array, z_coord_array,
																		x_spin_array, y_spin_array, z_spin_array,
																		forces_array_x, forces_array_y, forces_array_z,
																		fields_array_x, fields_array_y, fields_array_z);

		if (sim::time > sim::equilibration_time) 
		{

			const double Hx=sim::H_vec[0]*sim::H_applied;
			const double Hy=sim::H_vec[1]*sim::H_applied;
			const double Hz=sim::H_vec[2]*sim::H_applied;

			for(int i=start_index;i<end_index; i++)
			{
				fields_array_x[i]+=Hx;
				fields_array_y[i]+=Hy;
				fields_array_z[i]+=Hz;
			}
		}

          
		anisotropy::fields(atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array, atoms::type_array,
						fields_array_x, fields_array_y, fields_array_z,
						start_index, end_index, sim::temperature);
      	return;

    }

	namespace internal
	{

		void compute_exchange(const int start_index,
							const int end_index, 
							const std::vector<int>& neighbour_list_start_index,
							const std::vector<int>& neighbour_list_end_index,
							const std::vector<int>& type_array, 
							const std::vector<int>& neighbour_list_array, 
							const std::vector<double>& x_coord_array, 
							const std::vector<double>& y_coord_array,
							const std::vector<double>& z_coord_array,
							const std::vector<double>& x_spin_array, 
							const std::vector<double>& y_spin_array,
							const std::vector<double>& z_spin_array,
							std::vector<double>& forces_array_x, 
							std::vector<double>& forces_array_y,
							std::vector<double>& forces_array_z,
							std::vector<double>& fields_array_x,
							std::vector<double>& fields_array_y,
							std::vector<double>& fields_array_z)
		{

			double rx, ry, rz;
			double dx, dy, dz;
			double sx, sy, sz;
			double sjx, sjy,sjz;
			double si_dot_sj;
			double fx = 0.0, fy = 0.0, fz = 0.0;
			double hx = 0.0, hy = 0.0, hz = 0.0;
			double rji_sqr, rji, inv_rji; 
			double y, f_exch,  energy = 0.0;
			double J;
			int j;
			double r_sqr_cut = sld::internal::r_cut_fields*sld::internal::r_cut_fields;
	      	double exch_inv_rcut = 1.0/sld::internal::r_cut_fields;
       		double sumJ = 0.0;

       		for(int i=start_index;i<end_index; ++i)
			{

          		const unsigned int imat = atoms::type_array[i];

          		double exch_J0 = sld::internal::mp[imat].J0_ms.get();          // TODO : change the definition of J0_ms
          		double exch_J0_prime = sld::internal::mp[imat].J0_prime.get(); // TODO : change the definition of J0_prime
				
				int nbr_start = neighbour_list_start_index[i];
				int nbr_end = neighbour_list_end_index[i]+1;

				fx = 0.0;
				fy = 0.0;
				fz = 0.0;
				hx = 0.0;
				hy = 0.0;
				hz = 0.0;
				sumJ = 0.0;
				energy = 0.0;

				rx = x_coord_array[i];
				ry = y_coord_array[i];
				rz = z_coord_array[i];

				sx = x_spin_array[i];
				sy = y_spin_array[i];
				sz = z_spin_array[i];

				

          		for( int n = nbr_start; n < nbr_end; ++n)
				{

            		j = neighbour_list_array[n];

					if (j != i)
					{
						dx = -x_coord_array[j] + rx;
						dy = -y_coord_array[j] + ry;
						dz = -z_coord_array[j] + rz;

						dx = sld::PBC_wrap( dx, cs::system_dimensions[0], cs::pbc[0]);
						dy = sld::PBC_wrap( dy, cs::system_dimensions[1], cs::pbc[1]);
						dz = sld::PBC_wrap( dz, cs::system_dimensions[2], cs::pbc[2]);


						rji_sqr = dx*dx + dy*dy + dz*dz;

             			if( rji_sqr < r_sqr_cut)
            			{   

                 			rji = sqrt(rji_sqr);
                 			inv_rji = 1.0/ rji;

                 			y = (1.0 - rji * exch_inv_rcut);
                 			J = exch_J0 * y * y * y;

							sjx = x_spin_array[j];
							sjy = y_spin_array[j];
							sjz = z_spin_array[j];

							hx += J * sjx ;
							hy += J * sjy ;
							hz += J * sjz ;
							sumJ += J;

                 			si_dot_sj = sx * sjx + sy * sjy + sz * sjz;

							f_exch = -exch_J0_prime * y * y;

							fx += f_exch * dx *  (si_dot_sj)* inv_rji;
							fy += f_exch * dy *  (si_dot_sj)* inv_rji;
							fz += f_exch * dz *  (si_dot_sj) * inv_rji;

							energy += J *  (si_dot_sj);
             			}
          			}
       			}

				forces_array_x[i] += fx;
				forces_array_y[i] += fy;
				forces_array_z[i] += fz;

				fields_array_x[i] += hx;
				fields_array_y[i] += hy;
				fields_array_z[i] += hz;

				sld::internal::sumJ[i]=sumJ;
				sld::internal::exch_eng[i]=-0.5*energy;
   			}
   			return;
		}

		void compute_sld_coupling(const int start_index,
								const int end_index, 
								const std::vector<int>& neighbour_list_start_index,
								const std::vector<int>& neighbour_list_end_index,
								const std::vector<int>& type_array,
								const std::vector<int>& neighbour_list_array, 
								const std::vector<double>& x_coord_array, 
								const std::vector<double>& y_coord_array,
								const std::vector<double>& z_coord_array,
								const std::vector<double>& x_spin_array, 
								const std::vector<double>& y_spin_array,
								const std::vector<double>& z_spin_array,
								std::vector<double>& forces_array_x,
								std::vector<double>& forces_array_y,
								std::vector<double>& forces_array_z,
								std::vector<double>& fields_array_x, 
								std::vector<double>& fields_array_y,
								std::vector<double>& fields_array_z)
		{


			double rx, ry, rz;
			double dx, dy, dz;
			double sx, sy, sz;
			double sjx, sjy,sjz;
			double si_dot_sj;
			double fc_x = 0.0, fc_y = 0.0, fc_z = 0.0;
			double hc_x = 0.0, hc_y = 0.0, hc_z = 0.0;
			double rji_sqr, rji, inv_rji,  inv_rji2, inv_rji4, inv_rji6;
			double sj_dot_rji, si_dot_rji;
			double energy_c;
			int j; 

			double r_sqr_cut=sld::internal::r_cut_fields*sld::internal::r_cut_fields;
			double oneover3 = 1.0/3.0;
			double sumC;

			for(int i=start_index;i<end_index; ++i)
			{

				const unsigned int imat = atoms::type_array[i];
				double fact=sld::internal::mp[imat].C0.get()/1.602176634e-19;
				double fact_ms= sld::internal::mp[imat].C0_ms.get();

				int nbr_start = neighbour_list_start_index[i];
				int nbr_end = neighbour_list_end_index[i]+1;

				fc_x = 0.0;
				fc_y = 0.0;
				fc_z = 0.0;
				hc_x = 0.0;
				hc_y = 0.0;
				hc_z = 0.0;
				energy_c=0.0;
				sumC=0.0;

				rx = x_coord_array[i];
				ry = y_coord_array[i];
				rz = z_coord_array[i];

				sx = x_spin_array[i];
				sy = y_spin_array[i];
				sz = z_spin_array[i];

				for( int n = nbr_start; n < nbr_end; ++n)
				{
					j = neighbour_list_array[n];

					if (j != i)
					{

						dx = -x_coord_array[j] + rx;
						dy = -y_coord_array[j] + ry;
						dz = -z_coord_array[j] + rz;

						dx = sld::PBC_wrap( dx, cs::system_dimensions[0], cs::pbc[0]);
						dy = sld::PBC_wrap( dy, cs::system_dimensions[1], cs::pbc[1]);
						dz = sld::PBC_wrap( dz, cs::system_dimensions[2], cs::pbc[2]);


						rji_sqr = dx*dx + dy*dy + dz*dz;

						if(rji_sqr < r_sqr_cut)
						{

                            rji = sqrt(rji_sqr);
                            inv_rji = 1.0/ rji;

                            sjx = x_spin_array[j];
                            sjy = y_spin_array[j];
                            sjz = z_spin_array[j];

                            si_dot_sj = sx * sjx + sy * sjy + sz * sjz;

                            sj_dot_rji = (dx * sjx + dy * sjy + dz * sjz);
                            si_dot_rji = (dx * sx  + dy * sy  + dz * sz);

                            inv_rji2=inv_rji*inv_rji;
                            inv_rji4=inv_rji2*inv_rji2;
                            inv_rji6=inv_rji4*inv_rji2;

                            hc_x +=fact_ms*inv_rji4*(inv_rji2*dx*sj_dot_rji-oneover3*sjx);
                            hc_y +=fact_ms*inv_rji4*(inv_rji2*dy*sj_dot_rji-oneover3*sjy);
                            hc_z +=fact_ms*inv_rji4*(inv_rji2*dz*sj_dot_rji-oneover3*sjz);

							fc_x += fact*inv_rji6*( sj_dot_rji * sx + si_dot_rji * sjx -6.0* dx* sj_dot_rji* si_dot_rji * inv_rji2+ oneover3*4.0*si_dot_sj*dx);
                            fc_y += fact*inv_rji6*( sj_dot_rji * sy + si_dot_rji * sjy -6.0* dy* sj_dot_rji* si_dot_rji * inv_rji2+ oneover3*4.0*si_dot_sj*dy);
                            fc_z += fact*inv_rji6*( sj_dot_rji * sz + si_dot_rji * sjz -6.0* dz* sj_dot_rji* si_dot_rji * inv_rji2+ oneover3*4.0*si_dot_sj*dz);



                            sumC +=fact_ms*inv_rji4;
                  		}
               		}
            	}

				forces_array_x[i] += fc_x;
				forces_array_y[i] += fc_y;
				forces_array_z[i] += fc_z;

				fields_array_x[i] += hc_x;
				fields_array_y[i] += hc_y;
				fields_array_z[i] += hc_z;

				sld::internal::sumC[i]=sumC;
				sld::internal::coupl_eng[i]=-0.5*energy_c;

			}
            return;
        }


		void compute_sld_coupling_neel (const int start_index,
										const int end_index, // last +1 atom to be calculated
										const std::vector<int>& neighbour_list_start_index,
										const std::vector<int>& neighbour_list_end_index,
										const std::vector<int>& type_array, // type for atom
										const std::vector<int>& neighbour_list_array, // list of interactions between atom
										const std::vector<double>& x_coord_array, // coord vectors for atoms
										const std::vector<double>& y_coord_array,
										const std::vector<double>& z_coord_array,
										const std::vector<double>& x_spin_array, // coord vectors for atoms
										const std::vector<double>& y_spin_array,
										const std::vector<double>& z_spin_array,
										std::vector<double>& forces_array_x, //  vectors for forces
										std::vector<double>& forces_array_y,
										std::vector<double>& forces_array_z,
										std::vector<double>& fields_array_x, //  vectors for fields
										std::vector<double>& fields_array_y,
										std::vector<double>& fields_array_z)
		{
			double rx, ry, rz;
			double dx, dy, dz;
			double sx, sy, sz;
			double sjx, sjy,sjz;
			double si_dot_sj;
			double fc_x = 0.0, fc_y = 0.0, fc_z = 0.0;
			double hc_x = 0.0, hc_y = 0.0, hc_z = 0.0;
			double rji_sqr, rji, inv_rji,  inv_rji2, inv_rji4, inv_rji6;
			double sj_dot_rji, si_dot_rji;
			double energy_c;
			int j; 

			double r_sqr_cut=sld::internal::r_cut_fields*sld::internal::r_cut_fields;
			double oneover3=1.0/3.0;

			for(int i=start_index;i<end_index; ++i)
			{
				const unsigned int imat = atoms::type_array[i];
				double fact=sld::internal::mp[imat].C0.get()/1.602176634e-19;
				double fact_ms= sld::internal::mp[imat].C0_ms.get();

				fc_x = 0.0;
				fc_y = 0.0;
				fc_z = 0.0;
				hc_x = 0.0;
				hc_y = 0.0;
				hc_z = 0.0;
				energy_c=0.0;

				rx = x_coord_array[i];
				ry = y_coord_array[i];
				rz = z_coord_array[i];

				sx = x_spin_array[i];
				sy = y_spin_array[i];
				sz = z_spin_array[i];

				int nbr_start = neighbour_list_start_index[i];
				int nbr_end = neighbour_list_end_index[i]+1;


                for( int n = nbr_start; n < nbr_end; ++n)
				{
					j = neighbour_list_array[n];


					if ( j != i)
					{

						dx = -x_coord_array[j] + rx;
						dy = -y_coord_array[j] + ry;
						dz = -z_coord_array[j] + rz;

						dx = sld::PBC_wrap( dx, cs::system_dimensions[0], cs::pbc[0]);
						dy = sld::PBC_wrap( dy, cs::system_dimensions[1], cs::pbc[1]);
						dz = sld::PBC_wrap( dz, cs::system_dimensions[2], cs::pbc[2]);


						rji_sqr = dx*dx + dy*dy + dz*dz;

						if(rji_sqr < r_sqr_cut)
						{
							rji = sqrt(rji_sqr);
							inv_rji = 1.0/ rji;

							sjx = x_spin_array[j];
							sjy = y_spin_array[j];
							sjz = z_spin_array[j];

							si_dot_sj = sx * sjx + sy * sjy + sz * sjz;

							sj_dot_rji = (dx * sjx + dy * sjy + dz * sjz);
							si_dot_rji = (dx * sx  + dy * sy  + dz * sz);

							inv_rji2=inv_rji*inv_rji;
							inv_rji4=inv_rji2*inv_rji2;
							inv_rji6=inv_rji4*inv_rji2;

							double prod1=inv_rji2*si_dot_rji*si_dot_rji- oneover3*si_dot_sj;
							double prod2=inv_rji2*sj_dot_rji*sj_dot_rji- oneover3*si_dot_sj;
							double prod3=si_dot_rji*sj_dot_rji*sj_dot_rji*sj_dot_rji;
							double prod4=sj_dot_rji*si_dot_rji*si_dot_rji*si_dot_rji;
							double sj3=sj_dot_rji*sj_dot_rji*sj_dot_rji;
							double si3=si_dot_rji*si_dot_rji*si_dot_rji;
							double deriv1=2*inv_rji2*si_dot_rji;



                            hc_x +=12.0/35.0*fact_ms*inv_rji4*(inv_rji2*dx*sj_dot_rji-oneover3*sjx);
                            hc_y +=12.0/35.0*fact_ms*inv_rji4*(inv_rji2*dy*sj_dot_rji-oneover3*sjy);
                            hc_z +=12.0/35.0*fact_ms*inv_rji4*(inv_rji2*dz*sj_dot_rji-oneover3*sjz);

                            energy_c+= 12.0/35.0*fact_ms*inv_rji4*(inv_rji2*sj_dot_rji*si_dot_rji- oneover3*si_dot_sj);

                            hc_x +=9.0/5.0*fact_ms*inv_rji4*((deriv1*dx-oneover3*sjx)*prod2+prod1*(-oneover3*sjx));
                            hc_y +=9.0/5.0*fact_ms*inv_rji4*((deriv1*dy-oneover3*sjy)*prod2+prod1*(-oneover3*sjy));
                            hc_z +=9.0/5.0*fact_ms*inv_rji4*((deriv1*dz-oneover3*sjz)*prod2+prod1*(-oneover3*sjz));

                            energy_c +=9.0/5.0*fact_ms*inv_rji4*prod1*prod2;

							hc_x +=-2.0/5.0*fact_ms*inv_rji4*inv_rji4*(dx*sj3+3.0*dx*sj_dot_rji*si_dot_rji*si_dot_rji);
                            hc_y +=-2.0/5.0*fact_ms*inv_rji4*inv_rji4*(dy*sj3+3.0*dy*sj_dot_rji*si_dot_rji*si_dot_rji);
                            hc_z +=-2.0/5.0*fact_ms*inv_rji4*inv_rji4*(dz*sj3+3.0*dz*sj_dot_rji*si_dot_rji*si_dot_rji);

                            energy_c +=-2.0/5.0*fact_ms*inv_rji4*inv_rji4*(prod3+prod4);


                            fc_x += 12.0/35.0*fact*inv_rji6*( sj_dot_rji * sx + si_dot_rji * sjx -6.0* dx* sj_dot_rji* si_dot_rji * inv_rji2+ oneover3*4.0*si_dot_sj*dx);
                            fc_y += 12.0/35.0*fact*inv_rji6*( sj_dot_rji * sy + si_dot_rji * sjy -6.0* dy* sj_dot_rji* si_dot_rji * inv_rji2+ oneover3*4.0*si_dot_sj*dy);
                            fc_z += 12.0/35.0*fact*inv_rji6*( sj_dot_rji * sz + si_dot_rji * sjz -6.0* dz* sj_dot_rji* si_dot_rji * inv_rji2+ oneover3*4.0*si_dot_sj*dz);

                            fc_x += 9.0/5.0*fact*((-4*dx*inv_rji6)*prod1*prod2+ inv_rji4*prod2*(2*sx*inv_rji2*si_dot_rji-2*dx*si_dot_rji*si_dot_rji*inv_rji4) + inv_rji4*prod1*(2*sjx*inv_rji2*sj_dot_rji-2*dx*sj_dot_rji*sj_dot_rji*inv_rji4));
                            fc_y += 9.0/5.0*fact*((-4*dy*inv_rji6)*prod1*prod2+ inv_rji4*prod2*(2*sy*inv_rji2*si_dot_rji-2*dy*si_dot_rji*si_dot_rji*inv_rji4) + inv_rji4*prod1*(2*sjy*inv_rji2*sj_dot_rji-2*dy*sj_dot_rji*sj_dot_rji*inv_rji4));
                            fc_z += 9.0/5.0*fact*((-4*dz*inv_rji6)*prod1*prod2+ inv_rji4*prod2*(2*sz*inv_rji2*si_dot_rji-2*dz*si_dot_rji*si_dot_rji*inv_rji4) + inv_rji4*prod1*(2*sjz*inv_rji2*sj_dot_rji-2*dz*sj_dot_rji*sj_dot_rji*inv_rji4));

                            fc_x += -2.0/5.0*fact*((-4*dx*inv_rji6)*(inv_rji4*prod3+inv_rji4*prod4)+ inv_rji4*(-4*dx*inv_rji6*prod3 +inv_rji4*sx*sj3+inv_rji4*si_dot_rji*3*sjx*sj_dot_rji*sj_dot_rji)+ inv_rji4*(-4*dx*inv_rji6*prod4 +inv_rji4*sjx*si3+inv_rji4*sj_dot_rji*3*sx*si_dot_rji*si_dot_rji));
                            fc_y += -2.0/5.0*fact*((-4*dy*inv_rji6)*(inv_rji4*prod3+inv_rji4*prod4)+ inv_rji4*(-4*dy*inv_rji6*prod3 +inv_rji4*sy*sj3+inv_rji4*si_dot_rji*3*sjy*sj_dot_rji*sj_dot_rji)+ inv_rji4*(-4*dy*inv_rji6*prod4 +inv_rji4*sjy*si3+inv_rji4*sj_dot_rji*3*sy*si_dot_rji*si_dot_rji));
                            fc_z += -2.0/5.0*fact*((-4*dz*inv_rji6)*(inv_rji4*prod3+inv_rji4*prod4)+ inv_rji4*(-4*dz*inv_rji6*prod3 +inv_rji4*sz*sj3+inv_rji4*si_dot_rji*3*sjz*sj_dot_rji*sj_dot_rji)+ inv_rji4*(-4*dz*inv_rji6*prod4 +inv_rji4*sjz*si3+inv_rji4*sj_dot_rji*3*sz*si_dot_rji*si_dot_rji));
						}
               		}
            	}
				forces_array_x[i] += fc_x;
				forces_array_y[i] += fc_y;
				forces_array_z[i] += fc_z;

				fields_array_x[i] += hc_x;
				fields_array_y[i] += hc_y;
				fields_array_z[i] += hc_z;

				sld::internal::sumC[i]=hc_x*sx+hc_y*sy+hc_z*sz;
				sld::internal::coupl_eng[i]=-0.5*energy_c;
            }
            return;
        }
    } 
} 
