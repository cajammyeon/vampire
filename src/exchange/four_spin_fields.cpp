//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2023. All rights reserved.
//
//   Email: mara.strungaru@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "atoms.hpp" // for exchange list type defs
#include "sim.hpp"
#include "exchange.hpp"
#include "atoms.hpp"
#include "sld.hpp"

// exchange module headers
#include "internal.hpp"

double dot_product(double six, double siy, double siz, double sjx, double sjy, double sjz){
   double dp = six*sjx + siy*sjy + siz*sjz;
   return dp;
}

namespace exchange
{

	namespace internal
	{

		void four_spin_exchange_fields(const int start_index, // first atom for exchange interactions to be calculated
									const int end_index,                                 // field vectors for atoms
									std::vector<double>& field_array_x,
									std::vector<double>& field_array_y,
									std::vector<double>& field_array_z)
		{ 

			// =========================================================================
			// Memory allocation, allocate once before iteration, allow for overwrite
			// =========================================================================
			const int num_four_spin_neighbours = four_spin_neighbour_list_array_l.size();
			int atom, natomj, natomk, natoml;

			double Jij;
			double sjx, sjy, sjz;
			double skx, sky, skz;
			double slx, sly, slz;
			double sk_dot_sl, sj_dot_sk, sj_dot_sl;

			double athird = 1.0 / 3.0;

			double dist_a_x, dist_j_x, dist_k_x, dist_l_x;
			double dist_a_y, dist_j_y, dist_k_y, dist_l_y;
			double dist_a_z, dist_j_z, dist_k_z, dist_l_z;

			// =========================================================================
			// Iterate through the neighbour list for interaction calculation
			// =========================================================================
			for(int nn = 0; nn < num_four_spin_neighbours; ++nn)
			{
				atom   = four_spin_neighbour_list_array_i[nn];
				natomj = four_spin_neighbour_list_array_j[nn];
				natomk = four_spin_neighbour_list_array_k[nn];
				natoml = four_spin_neighbour_list_array_l[nn];

				dist_a_x = 1.0; 
				dist_j_x = 1.0; 
				dist_k_x = 1.0;
				dist_l_x = 1.0;

			 	dist_a_y = 1.0; 
				dist_j_y = 1.0;
				dist_k_y = 1.0;
				dist_l_y = 1.0;

				dist_a_z = 1.0;
				dist_j_z = 1.0; 
				dist_k_z = 1.0;
				dist_l_z = 1.0;

				// TODO : insert distance dependence here !!!!!
				if (internal::enable_fourspin_distance) {

					// ========== resolve x distance ==========
					dist_a_x = atoms::x_coord_array[atom];
					dist_j_x = atoms::x_coord_array[natomj];
					dist_k_x = atoms::x_coord_array[natomk];
					dist_l_x = atoms::x_coord_array[natoml];

					dist_j_x = 1 / sld::PBC_wrap(dist_j_x - dist_a_x, cs::system_dimensions[0], cs::pbc[0]);
					dist_k_x = 1 / sld::PBC_wrap(dist_k_x - dist_a_x, cs::system_dimensions[0], cs::pbc[0]);
					dist_l_x = 1 / sld::PBC_wrap(dist_l_x - dist_a_x, cs::system_dimensions[0], cs::pbc[0]);

					// ========== resolve y distance ==========
					dist_a_y = atoms::y_coord_array[atom];
					dist_j_y = atoms::y_coord_array[natomj];
					dist_k_y = atoms::y_coord_array[natomk];
					dist_l_y = atoms::y_coord_array[natoml];

					dist_j_y = 1 / sld::PBC_wrap(dist_j_y - dist_a_y, cs::system_dimensions[1], cs::pbc[1]);
					dist_k_y = 1 / sld::PBC_wrap(dist_k_y - dist_a_y, cs::system_dimensions[1], cs::pbc[1]);
					dist_l_y = 1 / sld::PBC_wrap(dist_l_y - dist_a_y, cs::system_dimensions[1], cs::pbc[1]);

					// ========== resolve z distance ==========
					dist_a_z = atoms::z_coord_array[atom];
					dist_j_z = atoms::z_coord_array[natomj];
					dist_k_z = atoms::z_coord_array[natomk];
					dist_l_z = atoms::z_coord_array[natoml];

					dist_j_z = 1 / sld::PBC_wrap(dist_j_z - dist_a_z, cs::system_dimensions[2], cs::pbc[2]);
					dist_k_z = 1 / sld::PBC_wrap(dist_k_z - dist_a_z, cs::system_dimensions[2], cs::pbc[2]);
					dist_l_z = 1 / sld::PBC_wrap(dist_l_z - dist_a_z, cs::system_dimensions[2], cs::pbc[2]);
				}

				Jij = four_spin_exchange_list[nn];

				sjx = atoms::x_spin_array[natomj] * dist_j_x;
				sjy = atoms::y_spin_array[natomj] * dist_j_y;
				sjz = atoms::z_spin_array[natomj] * dist_j_z;

				skx = atoms::x_spin_array[natomk] * dist_k_x;
				sky = atoms::y_spin_array[natomk] * dist_k_y;
				skz = atoms::z_spin_array[natomk] * dist_k_z;

				slx = atoms::x_spin_array[natoml] * dist_l_x;
				sly = atoms::y_spin_array[natoml] * dist_l_y;
				slz = atoms::z_spin_array[natoml] * dist_l_z;

				sk_dot_sl = dot_product(skx,sky,skz,slx,sly,slz);
				sj_dot_sk = dot_product(skx,sky,skz,sjx,sjy,sjz);
				sj_dot_sl = dot_product(sjx,sjy,sjz,slx,sly,slz);

				field_array_x[atom] += (Jij*athird) * (sjx * sk_dot_sl + skx * sj_dot_sl + slx * sj_dot_sk);
				field_array_y[atom] += (Jij*athird) * (sjy * sk_dot_sl + sky * sj_dot_sl + sly * sj_dot_sk);
				field_array_z[atom] += (Jij*athird) * (sjz * sk_dot_sl + skz * sj_dot_sl + slz * sj_dot_sk);
			}

			return;

		}

	}
}
