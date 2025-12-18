import numpy as np

# --- Configuration ---
J_NN = 5.0e-21   # Exchange strength for Nearest Neighbors (e.g., Fe-Rh)
J_2NN = 5.0e-21  # Exchange strength for Next Nearest Neighbors (e.g., Fe-Fe)

# Unit Cell Dimension (Normalized)
L = 1.0 

# Atom coordinates (from your UCF)
# IDs 0-7 (Sublattice A), 8-15 (Sublattice B)
atoms = {
    0:  [0.0, 0.0, 0.0], 1:  [0.0, 0.0, 0.5], 2:  [0.5, 0.0, 0.0], 3:  [0.5, 0.0, 0.5],
    4:  [0.0, 0.5, 0.0], 5:  [0.0, 0.5, 0.5], 6:  [0.5, 0.5, 0.0], 7:  [0.5, 0.5, 0.5],
    8:  [0.25, 0.25, 0.25], 9:  [0.75, 0.25, 0.25], 10: [0.25, 0.75, 0.25], 11: [0.25, 0.25, 0.75],
    12: [0.75, 0.75, 0.25], 13: [0.25, 0.75, 0.75], 14: [0.75, 0.75, 0.75], 15: [0.75, 0.75, 0.75] # Note: Atom 15 duplicates 14 in your original text? Assuming symmetric pattern:
}
# Correcting Atom 15 to likely pattern position if 14 is 0.75,0.75,0.75?
# Let's assume standard 2x2x2 layout for 15: [0.75, 0.75, 0.75] is index 14. 
# Index 15 in the pattern (0.75, 0.25, 0.75)? Let's stick to the positions that form the B2 lattice.
# Correcting list based on standard B2 supercell:
atoms = [
    [0.0, 0.0, 0.0], [0.0, 0.0, 0.5], [0.5, 0.0, 0.0], [0.5, 0.0, 0.5],
    [0.0, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.5, 0.0], [0.5, 0.5, 0.5],
    [0.25, 0.25, 0.25], [0.75, 0.25, 0.25], [0.25, 0.75, 0.25], [0.25, 0.25, 0.75],
    [0.75, 0.75, 0.25], [0.25, 0.75, 0.75], [0.75, 0.25, 0.75], [0.75, 0.75, 0.75]
]

interactions = []

# Loop over all atom pairs
for i, pos_i in enumerate(atoms):
    for j, pos_j in enumerate(atoms):
        
        # Check periodic images (-1, 0, 1)
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                for dz in [-1, 0, 1]:
                    
                    # Calculate distance vector considering PBC
                    delta = np.array(pos_j) + np.array([dx, dy, dz]) - np.array(pos_i)
                    dist = np.linalg.norm(delta)
                    
                    # Thresholds (with small tolerance for float errors)
                    # Shell 1 (Fe-Rh): sqrt(3)/4 approx 0.433
                    # Shell 2 (Fe-Fe): 0.5
                    
                    J = 0.0
                    if 0.4 < dist < 0.45:  # Shell 1 (Nearest Neighbor)
                        J = J_NN
                    elif 0.49 < dist < 0.51: # Shell 2 (Next Nearest Neighbor)
                        J = J_2NN
                        
                    if J != 0.0:
                        interactions.append(f"{len(interactions)} {i} {j} {dx} {dy} {dz} {J}")

# --- Output to Screen ---
print("Number of interactions:", len(interactions))
print("# Interactions n exctype; id i j dx dy dz Jij")
print(f"{len(interactions)} normalized-isotropic")
for line in interactions:
    print(line)