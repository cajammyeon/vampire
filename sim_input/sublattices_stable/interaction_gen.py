import numpy as np

# --- CONFIGURATION ---
# Unit cell size (Angstroms) for FeRh (approx doubled Z for AFM)
uc_size = np.array([3, 3, 6]) 

# Atoms: ID, Fractional Pos (x,y,z), Name
# 0=Fe_Up, 1=Fe_Down, 2=Rh
atoms = [
    {'id': 0, 'pos': np.array([0.0, 0.0, 0.0]), 'mat': 'Fe1'},
    {'id': 1, 'pos': np.array([0.0, 0.0, 0.5]), 'mat': 'Fe2'},
    {'id': 2, 'pos': np.array([0.5, 0.5, 0.25]), 'mat': 'Rh'},
]

# Interaction Cutoff (Angstroms)
# 5.5 A captures neighbors sufficient for Fe-Rh and Fe-Fe interactions
cutoff = 7.5

# --- GENERATION ---
interactions = []
int_id = 0

for atom_i in atoms:
    for atom_j in atoms:
        # Search surrounding unit cells (dx, dy, dz from -2 to +2)
        for dx in range(-2, 3):
            for dy in range(-2, 3):
                for dz in range(-2, 3):
                    
                    # Calculate real distance
                    shift = np.array([dx, dy, dz])
                    pos_j_real = (atom_j['pos'] + shift) * uc_size
                    pos_i_real = atom_i['pos'] * uc_size
                    dist = np.linalg.norm(pos_j_real - pos_i_real)
                    
                    # Exclude self-interaction (dist=0) and check cutoff
                    if 0.01 < dist <= cutoff:
                        # SET INTERACTION TO SIMPLE 1
                        J_val = 1.0 
                        
                        interactions.append([int_id, atom_i['id'], atom_j['id'], dx, dy, dz, J_val])
                        int_id += 1

# --- OUTPUT ---
print(f"# Interactions n exctype; id i j dx dy dz Jij")
print(f"{len(interactions)} normalised-isotropic")
for i in interactions:
    # Print formatted columns
    print(f"{i[0]}\t{i[1]}\t{i[2]}\t{i[3]}\t{i[4]}\t{i[5]}\t{i[6]}")