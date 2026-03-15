find . -maxdepth 1 \( -name "*.pdb" -o -name "*.cif" \) -exec basename {} \; | sed 's/\.pdb$//; s/\.cif$//' | sort
