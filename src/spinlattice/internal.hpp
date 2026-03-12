//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2022. All rights reserved.
//
//   Email: mara.strungaru@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef SLD_INTERNAL_H_
#define SLD_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the sld module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <vector>

// Vampire headers
#include "sld.hpp"

// sld module headers
#include "internal.hpp"


namespace sld
{
 	namespace internal
	{

      	class set_double_t
		{

     	private:
        	double value; 
        	bool setf; 
		public:
			set_double_t() : value(0.0), setf(false) {}

			void set(double in_value)
			{
				value = in_value;
				setf = true;
			};

			double get() {return value;};
			bool is_set(){return setf; };

      };

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-----------------------------------------------------------------------------
      // internal materials class for storing material parameters
      //-----------------------------------------------------------------------------
      	class mp_t
		{

          	private:

          	public:
				//------------------------------
				// material parameter variables
				//------------------------------
				set_double_t mass;
				set_double_t V0;
				set_double_t J0;
				set_double_t C0;
				set_double_t damp_lat;
				set_double_t eq_damp_lat;

				set_double_t J0_ms;
				set_double_t C0_ms;
				set_double_t J0_prime;
				set_double_t F_th_sigma;
				set_double_t F_th_sigma_eq;

				// constructor
				mp_t (const unsigned int max_materials = 100)
				{
					mass.set(5.7915e-3);
					V0.set(0.15);
					J0.set(1);
					J0_prime.set(1/1.602176634e-19);
					J0_ms.set(1/2.04028e-23);
					C0.set(1/1.602176634e-19);
					C0_ms.set(0.5/2.04028e-23);
					F_th_sigma.set(1.0);
					F_th_sigma_eq.set(1.0);
					damp_lat.set(0.06);
					eq_damp_lat.set(0.6);
				};

       	};
		//-------------------------------------------------------------------------
		// Internal shared variables
		//-------------------------------------------------------------------------

		extern bool enabled;                        // bool to enable module
		extern std::vector<sld::internal::mp_t> mp; // array of material properties

		extern double r_cut_pot;    // mechanical potential cutoff
		extern double r_cut_fields; // exchange/coupling cutoff

		extern double dr_init;
		extern double th_velo;

		extern double morse_beta;
		extern double morse_factor;
		extern double alpha_m;
		extern double r0_m;
		extern double morse_D;
		extern bool morse;

		extern bool harmonic; 
		extern bool pseudodipolar;
		extern bool full_neel;

		extern std::vector<double> x0_coord_array;
		extern std::vector<double> y0_coord_array;
		extern std::vector<double> z0_coord_array;

		extern std::vector <double> x_coord_storage_array;
		extern std::vector <double> y_coord_storage_array;
		extern std::vector <double> z_coord_storage_array;

		extern std::vector<double> forces_array_x;
		extern std::vector<double> forces_array_y;
		extern std::vector<double> forces_array_z;

		extern std::vector<double> fields_array_x;
		extern std::vector<double> fields_array_y;
		extern std::vector<double> fields_array_z;

		extern std::vector<double> velo_array_x;
		extern std::vector<double> velo_array_y;
		extern std::vector<double> velo_array_z;
		extern std::vector<double> potential_eng;
		extern std::vector<double> exch_eng;
		extern std::vector<double> coupl_eng;
		extern std::vector<double> sumJ;
		extern std::vector<double> sumC;

		extern std::vector<int> test_atom_list; //Core atoms of each octant

		void initialise_positions(std::vector<double>& x0_coord_array, // coord vectors for atoms
								std::vector<double>& y0_coord_array,
								std::vector<double>& z0_coord_array,
								std::vector<double>& x_coord_array, // coord vectors for atoms
								std::vector<double>& y_coord_array,
								std::vector<double>& z_coord_array,
								const double dr);

		void thermal_velocity(std::vector<double>& x_velo_array, // coord vectors for atoms
							std::vector<double>& y_velo_array,
							std::vector<double>& z_velo_array);

     	void initialise_sld_parameters();

		void compute_forces_harmonic(const int start_index,
									const int end_index, // last +1 atom to be calculated
									const std::vector<int>& neighbour_list_start_index,
									const std::vector<int>& neighbour_list_end_index,
									const std::vector<int>& type_array, // type for atom
									const std::vector<int>& neighbour_list_array, // list of interactions between atom
									const std::vector<double>& x0_coord_array, // coord vectors for atoms
									const std::vector<double>& y0_coord_array,
									const std::vector<double>& z0_coord_array,
									const std::vector<double>& x_coord_array, // coord vectors for atoms
									const std::vector<double>& y_coord_array,
									const std::vector<double>& z_coord_array,
									std::vector<double>& forces_array_x, //  vectors for forces
									std::vector<double>& forces_array_y,
									std::vector<double>& forces_array_z,
									std::vector<double>& potential_eng);

    	void compute_exchange(const int start_index,
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
							std::vector<double>& fields_array_x, //  vectors for forces
							std::vector<double>& fields_array_y,
							std::vector<double>& fields_array_z);

		void compute_forces_morse(const int start_index,
								const int end_index, // last +1 atom to be calculated
								const std::vector<int>& neighbour_list_start_index,
								const std::vector<int>& neighbour_list_end_index,
								const std::vector<int>& type_array, // type for atom
								const std::vector<int>& neighbour_list_array, // list of interactions between atom
								const std::vector<double>& x_coord_array, // coord vectors for atoms
								const std::vector<double>& y_coord_array,
								const std::vector<double>& z_coord_array,
								std::vector<double>& forces_array_x, //  vectors for forces
								std::vector<double>& forces_array_y,
								std::vector<double>& forces_array_z,
								std::vector<double>& potential_eng);

		void compute_sld_coupling(const int start_index,
								const int end_index, // last +1 atom to be calculated
								const std::vector<int>& neighbour_list_start_index,
								const std::vector<int>& neighbour_list_end_index,
								const std::vector<int>& type_array, // type for atom
								const std::vector<int>& neighbour_list_array, // list of interactions between atom
								const std::vector<double>& x_coord_array, // coord vectors for atoms
								const std::vector<double>& y_coord_array,
								const std::vector<double>& z_coord_array,
								const std::vector<double>& x_spin_array, // spin  vectors for atoms
								const std::vector<double>& y_spin_array,
								const std::vector<double>& z_spin_array,
								std::vector<double>& forces_array_x, //  vectors for forces
								std::vector<double>& forces_array_y,
								std::vector<double>& forces_array_z,
								std::vector<double>& fields_array_x, //  vectors for forces
								std::vector<double>& fields_array_y,
								std::vector<double>& fields_array_z);

		void compute_sld_coupling_neel(const int start_index,
									const int end_index, // last +1 atom to be calculated
									const std::vector<int>& neighbour_list_start_index,
									const std::vector<int>& neighbour_list_end_index,
									const std::vector<int>& type_array, // type for atom
									const std::vector<int>& neighbour_list_array, // list of interactions between atom
									const std::vector<double>& x_coord_array, // coord vectors for atoms
									const std::vector<double>& y_coord_array,
									const std::vector<double>& z_coord_array,
									const std::vector<double>& x_spin_array, // spin  vectors for atoms
									const std::vector<double>& y_spin_array,
									const std::vector<double>& z_spin_array,
									std::vector<double>& forces_array_x, //  vectors for forces
									std::vector<double>& forces_array_y,
									std::vector<double>& forces_array_z,
									std::vector<double>& fields_array_x, //  vectors for forces
									std::vector<double>& fields_array_y,
									std::vector<double>& fields_array_z);

		void cayley_update(const int start_index,
						const int end_index,
						double dt,
						std::vector<double>& x_spin_array, // coord vectors for atoms
						std::vector<double>& y_spin_array,
						std::vector<double>& z_spin_array,
						std::vector<double>& fields_array_x, //  vectors for fields
						std::vector<double>& fields_array_y,
						std::vector<double>& fields_array_z);

		void add_spin_noise(const int start_index,
							const int end_index,
							double dt,
							const std::vector<int>& type_array, // type for atom
							const std::vector<double>& x_spin_array, // coord vectors for atoms
							const std::vector<double>& y_spin_array,
							const std::vector<double>& z_spin_array,
							std::vector<double>& fields_array_x, //  vectors for fields
							std::vector<double>& fields_array_y,
							std::vector<double>& fields_array_z,
							std::vector<double>& Hx_th, //  vectors for fields
							std::vector<double>& Hy_th,
							std::vector<double>& Hz_th);

		//MPI variables
		extern std::vector<std::vector<int> > c_octants; //Core atoms of each octant
		extern std::vector<std::vector<int> > b_octants; //Boundary atoms of each octant

		extern std::vector <int> all_atoms_octant_start_index;
		extern std::vector <int> all_atoms_octant_end_index;
		extern std::vector <int> all_atoms_octant;

   	} 
}

#endif //SLD_INTERNAL_H_