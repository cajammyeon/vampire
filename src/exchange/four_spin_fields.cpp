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

			// =========================================================================
			// Iterate through the neighbour list for interaction calculation
			// =========================================================================
			for(int nn = 0; nn < num_four_spin_neighbours; ++nn)
			{
				atom   = four_spin_neighbour_list_array_i[nn];
				natomj = four_spin_neighbour_list_array_j[nn];
				natomk = four_spin_neighbour_list_array_k[nn];
				natoml = four_spin_neighbour_list_array_l[nn];

				Jij = four_spin_exchange_list[nn];

				sjx = atoms::x_spin_array[natomj];
				sjy = atoms::y_spin_array[natomj];
				sjz = atoms::z_spin_array[natomj];

				skx = atoms::x_spin_array[natomk];
				sky = atoms::y_spin_array[natomk];
				skz = atoms::z_spin_array[natomk];

				slx = atoms::x_spin_array[natoml];
				sly = atoms::y_spin_array[natoml];
				slz = atoms::z_spin_array[natoml];

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
