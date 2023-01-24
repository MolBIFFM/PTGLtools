/*
 * This file is part of the PTGLtools software package.
 *
 * Copyright Tim Schäfer 2012. PTGLtools is free software, see the LICENSE and README files for details.
 *
 * @author ts
 */
package tools;

import proteinstructure.Atom;
import java.util.ArrayList;
import tools.DP;

/**
 * Some helper functions for computing various 3D stuff. Used for visualization scripts and other stuff.
 * @author ts
 */
public class Comp3DTools {
    
    /**
     * Computes the center atoms of a list of atoms, i.e., the atom with the minimal max distance to any other atom in the list. Uses the 3D coords from the PDB file.
     * @param atoms the atom list
     * @return the resulting center atom
     */
    public static Atom getCenterAtomOf(ArrayList<Atom> atoms) {
        Atom a, b, center;
        a = b = center = null;
        int maxDistForAtom, dist = 0; // just assign a small start value
        int MAXDIST = Integer.MAX_VALUE;   // just assign a *very* large start value
        int totalMinMaxDist = MAXDIST;
        //Integer atomRadius = Settings.getInteger("PTGLgraphComputation_I_aa_atom_radius");

        if(atoms.size() < 1) {
            DP.getInstance().w("empty atom list");
            return(null);
        }
                
        for(Integer i = 0; i < atoms.size(); i++) {

            a = atoms.get(i);
            maxDistForAtom = 0;

            for(Integer j = 0; j < atoms.size(); j++) {     // we need to compare the atom to itself (distance = 0 then) if there only is a single atom in this residue (which holds for ligands like 'MG'). So j=i, not j=i+1.

                b = atoms.get(j);
                dist = a.distToAtom(b);

                if(dist > maxDistForAtom) {
                    maxDistForAtom = dist;
                }
            }

            // We determined the maximal distance of this atom to any other atom of this residue.
            // Now check whether this maxDist is smaller than the smallest current maxDist.

            if(maxDistForAtom < totalMinMaxDist) {
                totalMinMaxDist = maxDistForAtom;

                // Also update the current center atom. We can't break here though because
                //  this may still get improved/overwritten during the rest of the loop.
                center = a;
            }

        }
             
        return(center);
    }
    
}
