import sys, os
from pymol import cmd, stored, editing

@cmd.extend
def scanFactor( nucleosome, probe, cutoff=2, clashKeep=1, writeModels=False ):
  """
USAGE

  scanFactor nucleosome object,
             probe object,
             probe DNA chain ID,
             distance in Ångstrom below which atoms are considered clashing,
             keep models with fewer than the specified number of residue clashes

  """

  results=[]
  if writeModels:
    dirName =  None
    dirName = ("%s_%s_models" % (nucleosome, probe) )
    if not os.path.exists(dirName):
      os.mkdir("%s" % (dirName) )


  if not cmd.select("polymer.nucleic and %s" % (nucleosome) ) :
    print("\nThe nucleosome object \'%s\' must contain contain DNA – are you sure this is a nucleosome model?\n" % (nucleosome) )
    return

  if not cmd.select("polymer.nucleic and %s" % (probe) ) :
    print("\nThe probe object \'%s\' must contain at least two DNA bases for superposition, the first two of which should correspond\nto the position of the central two bases of the factor's recognition sequence when bound to DNA\n" % (probe) )
    return
      

  # Determine probe DNA chain ID
  cmd.select("probeDNA","polymer.nucleic and %s" % (probe) )
  probeDNAChain = cmd.get_model("probeDNA").atom[0].chain

  # Renumber probe DNA chain to start from 1
  firstBase = []
  firstBase = int(cmd.get_model("%s and chain %s" % (probe, probeDNAChain)).atom[0].resi)
  secondBase = firstBase+1
  offset = firstBase-1

  print ("Renumbering %s chain %s to start at residue 1\n" % (probe, probeDNAChain) )
  cmd.alter("%s and chain %s" % (probe, probeDNAChain), "resi=str(int(resi)-%d)" % (offset) )

  # Remove hydrogens
  print ("\nRemoving hydrogens\n")
  cmd.remove("hydrogens")
  
  # Determine nucleosomal DNA chain IDs
  chains = []
  chains = cmd.get_chains("polymer.nucleic and %s" % (nucleosome) )

  # Loop over nuclesomal DNA
  for chain in chains:

     firstBase = int(cmd.get_model("%s and chain %s" % (nucleosome, chain)).atom[0].resi)
     chainLength = (len(cmd.get_model("polymer and %s and chain %s" % (nucleosome, chain)).get_residues()))
     
     for i in range(firstBase, (firstBase + chainLength) ):

        j=(i+1)
        
        # Handle lack of 5' phosphates (such as on the first base)
        if not cmd.select("%s and chain %s and resi %d and name P" % (nucleosome, chain, i) ) :
           cmd.select("moving", "%s and chain %s and resi 2 and backbone and not (n. P or n. OP1 or n. OP2)" % (probe, probeDNAChain) )
           cmd.select("target", "%s and chain %s and resi %d and backbone and not (n. P or n. OP1 or n. OP2)" % (nucleosome, chain, j) )
        # Handle overhanging superposition for final base
        elif i == (firstBase + chainLength - 1):
           cmd.select("moving", "%s and chain %s and resi 1 and backbone" % (probe, probeDNAChain) )
           cmd.select("target", "%s and chain %s and resi %d and backbone" % (nucleosome, chain, i) )
        else:
           cmd.select("moving", "%s and chain %s and resi 1-2 and backbone" % (probe, probeDNAChain) )
           cmd.select("target", "%s and chain %s and resi %d-%d and backbone" % (nucleosome, chain, i, j) )

        #numAtomMoving=(cmd.count_atoms("moving"))
        #print ("Number of moving atoms: %d" % (numAtomMoving) )
        #numAtomTarget=(cmd.count_atoms("target"))
        #print ("Number of target atoms: %d" % (numAtomTarget) ) 

        alignResult=[]
        alignResult = cmd.pair_fit("moving", "target")
        
        if alignResult:
        
          # alignResult = cmd.super("moving", "target")
          # print ("RMSD after refinement %.4f with %d atoms in alignhment after %d refinement cycles" % (alignResult[0],alignResult[1],alignResult[2]) )
          # print(" ".join('%s' % x for x in alignResult))

          clashes=("clashes")
          # Calculate clashes excluding DNA portion of probe
          cmd.select("%s" % (clashes), "(%s and not chain %s) within %f of %s" % (probe, probeDNAChain, float(cutoff), nucleosome) )
          scoreAtoms = (cmd.count_atoms("%s" % (clashes) ))
          cmd.select("%s" % (clashes), "clashes byres clashes")
          scoreRes = (cmd.count_atoms("%s" % (clashes) + " and name CA"))

          # Write out models (nuclesome + current superimposed probe) for Rosetta scoring
          if  writeModels :
            modelName = None
            modelName = ("%s_probe_%s_%d-%d" % (nucleosome, chain, i, j) )

            editing.copy_to("%s" % (modelName), "(%s or %s)" % (nucleosome, probe) )
            cmd.save("%s/%s.pdb" % (dirName,modelName), "(%s)" % (modelName) )
            cmd.delete("%s" % (modelName) )

         # Retain models with no. residue clashes less than or equal to 'clashKeep'
          if scoreRes < (int(clashKeep)+1):
             editing.copy_to("probe_%s_%d-%d" % (chain, i, j), "%s" % (probe) )
             cmd.group("ClashKeep %s)" % (clashKeep),"nonclashing + probe_%s_%d-%d" % (chain, i, j) )
             cmd.order("probe_*", "yes")

          print ("DNA chain %s bases %d - %d : %s clashes with %d non-hydrogen atoms in %d residues" % (chain, i, j, probe, scoreAtoms, scoreRes) )
          clashes=("%s,%d,%d,%d,%d" % (chain, i, j, scoreAtoms, scoreRes) )
          results.append((clashes))
        # Catch superpostion failures
        else:
            
          print ("DNA chain %s bases %d - %d : %s superpositon failed and clashes were not calculated – check for missing atoms, alternate conformations, or unsual bases at these positions" % (chain, i, j, probe) )
          clashes=("%s,%d,%d,-,-" % (chain, i, j) )
          results.append((clashes))


  with open('%s_scanFactor.csv' % (nucleosome), 'w') as output:
    output.write("DNAchain,start,end,atomClashes,residueClashes\n")
    for i in results:
      output.write("%s\n" % i)
  print ("\nOutput saved to %s_scanFactor.csv\n" % (nucleosome) )
