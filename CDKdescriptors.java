package main;

import org.openscience.cdk.RingSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.isomorphism.VentoFoggia;
import org.openscience.cdk.qsar.descriptors.molecular.JPlogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.VABCDescriptor;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

public class CDKdescriptors {
	
	public SmilesParser   sp  = new SmilesParser(SilentChemObjectBuilder.getInstance());
	
	/**
	 * Calculating volume decriptor.
	 * 
	 * @param ac	IAtomContainer
	 * @return		boolean
	 */
	
	public String volumeDescriptor(IAtomContainer ac) {
		VABCDescriptor volumeDesc = new VABCDescriptor();
		return volumeDesc.calculate(ac).getValue().toString();
	}
	
	/**
	 * Calculating log p with Mannhold method.
	 * 
	 * @param ac	IAtomContainer	
	 * @return		boolean
	 */
	
	public String MannholdLogPDescriptor(IAtomContainer ac) {
		JPlogPDescriptor logp = new JPlogPDescriptor();
		return logp.calculate(ac).getValue().toString();
	}

	/**
	 * Searching for triple bonds in input atom container.
	 * 
	 * @param ac	IAtomContainer	
	 * @return		boolean
	 */
	
	public boolean detectTripleBonds(IAtomContainer ac) {
		boolean check = false;
		for(IBond b: ac.bonds()) {
			if(b.getOrder().ordinal()==2) {
				check=true;
				break;
			}
		}
		return check;
	}
	
	/**
	 * Searching for allkenes in input atom container.
	 * 
	 * @param ac	IAtomContainer 
	 * @return		boolean
	 * @throws InvalidSmilesException
	 */
	
	public boolean detectAllenes(IAtomContainer ac) throws InvalidSmilesException {
		boolean check = false;
		if(substructureSearch(sp.parseSmiles("C=C=C"),ac)) check=true;
		return check;
	}
	
	/**
	 * Substructure search in an input atom container.
	 * 
	 * @param sub		IAtomContainer
	 * @param ac		IAtomContainer
	 * @return			boolean
	 * @throws InvalidSmilesException
	 */
	
	public boolean substructureSearch(IAtomContainer sub, IAtomContainer ac) throws InvalidSmilesException {
		boolean contains= false;
		Pattern	pattern = VentoFoggia.findSubstructure(sub);
		int[] match = pattern.match(ac);
		if (match.length > 0) contains=true;
		return contains;
	}
	
	/**
	 * Checking whether the molecule includes 3 or 4-membered rings.
	 * 
	 * @param ac	IAtomContainer
	 * @return 		boolean 
	 */
	
	public boolean includes3or4rings(IAtomContainer ac) {
		boolean check=false;
		AllRingsFinder ringsFinder = new AllRingsFinder();
		RingSet search = new RingSet();
		IRingSet ringSet = ringsFinder.findAllRings(ac);
		search.add(ringSet);
		
		int atomCount=0;
	    for (IAtomContainer container : ringSet.atomContainers()) {
	    	atomCount=container.getAtomCount();
	    	if(atomCount==3 || atomCount==4) {
	    		check=true;
	        	break;
	        }   
	    }
	    return check;
	}
}
