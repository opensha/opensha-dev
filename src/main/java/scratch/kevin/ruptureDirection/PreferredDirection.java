package scratch.kevin.ruptureDirection;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

import org.opensha.commons.data.Named;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RuptureTreeNavigator;
import org.opensha.sha.faultSurface.EvenlyGriddedSurface;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;

import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;

import scratch.UCERF3.erf.ETAS.SeisDepthDistribution;

public class PreferredDirection implements Named {

	private String name;
	private FaultSystemSolution sol;
	private FaultSystemRupSet rupSet;
	private List<FaultSection> targetSects;
	private List<EvenlyGriddedSurface> targetSectSurfs;
	private List<Boolean> alongStrikes;
	private double targetFractInPreferred;
	private double maxBilateralFract;
	
	private HashSet<Integer> targetSectIDs;
	private ImmutableSet<Integer> affectedRuptures;
	
	private List<EvenlyDiscretizedFunc> targetDepthCmlWeights;
	private ClusterRuptures cRups;
	
	private Random rand;

	/**
	 * 
	 * @param name
	 * @param sol fault system solution
	 * @param targetSects sections sharing this preferred rupture direction
	 * @param alongStrikes boolean for each section stating if the preferred direction is along strike (true) or not (false)
	 * @param targetFractInPreferred the target fraction of ruptures nucleating on these sections that should propagate
	 * in the preferred dirction
	 * @param maxBilateralFract maximum fraction of a rupture that can go bilaterally in the non-preferred direction for
	 * it to still be called a preferred direction hypocenter
	 * @param aseisReducesArea true if aseismicity reduces area of fault sections
	 */
	public PreferredDirection(String name, FaultSystemSolution sol, List<FaultSection> targetSects, List<Boolean> alongStrikes,
			double targetFractInPreferred, double maxBilateralFract, boolean aseisReducesArea) {
		this.name = name;
		this.sol = sol;
		this.rupSet = sol.getRupSet();
		this.cRups = rupSet.requireModule(ClusterRuptures.class);
		Preconditions.checkState(!targetSects.isEmpty(), "No sections supplied");
		this.targetSects = targetSects;
		Preconditions.checkState(alongStrikes.size() == targetSects.size(), "Along strikes list is of different size");
		this.alongStrikes = alongStrikes;
		Preconditions.checkState(targetFractInPreferred > 0d && targetFractInPreferred <= 1d,
				"Bad target fract (must be >0 and <=1): %s", targetFractInPreferred);
		this.targetFractInPreferred = targetFractInPreferred;
		Preconditions.checkState(maxBilateralFract > 0d && maxBilateralFract < 1d,
				"Bad max bilateral fract (must be >0 and <1): %s", maxBilateralFract);
		this.maxBilateralFract = maxBilateralFract;
		
		ImmutableSet.Builder<Integer> affectedRuptures = ImmutableSet.builder();
		targetSectIDs = new HashSet<>();
		targetSectSurfs = new ArrayList<>();
		for (FaultSection sect : targetSects) {
			affectedRuptures.addAll(rupSet.getRupturesForSection(sect.getSectionId()));
			targetSectIDs.add(sect.getSectionId());
			targetSectSurfs.add((EvenlyGriddedSurface)sect.getFaultSurface(0.1d, false, aseisReducesArea));
		}
		this.affectedRuptures = affectedRuptures.build();
		
		SeisDepthDistribution seisDepthDistribution = new SeisDepthDistribution();
		double delta=2;
		HistogramFunction binnedDepthDistFunc = new HistogramFunction(1d, 12,delta);
		for(int i=0;i<binnedDepthDistFunc.size();i++) {
			double prob = seisDepthDistribution.getProbBetweenDepths(binnedDepthDistFunc.getX(i)-delta/2d,binnedDepthDistFunc.getX(i)+delta/2d);
			binnedDepthDistFunc.set(i,prob);
		}
		binnedDepthDistFunc.normalizeBySumOfY_Vals();
		targetDepthCmlWeights = new ArrayList<>();
		for (EvenlyGriddedSurface surf : targetSectSurfs) {
			int midCol = surf.getNumCols()/2;
			double topDepth = surf.getLocation(0, midCol).depth;
			double botDepth = surf.getLocation(surf.getNumRows()-1, midCol).depth;
			EvenlyDiscretizedFunc depthWeights = new EvenlyDiscretizedFunc(topDepth, botDepth, surf.getNumRows());
			for (int i=0; i<depthWeights.size(); i++) {
				double depth = depthWeights.getX(i);
				double minDepth = Math.max(0d, depth-depthWeights.getDelta()*0.5);
				double maxDepth = Math.min(24d, depth+depthWeights.getDelta()*0.5);
				if (minDepth < maxDepth)
					depthWeights.set(i, seisDepthDistribution.getProbBetweenDepths(minDepth, maxDepth));
			}
			double sumWeights = depthWeights.calcSumOfY_Vals();
			if (sumWeights == 0d) {
				// zero weight, probably too deep, just distribute uniformly
				for (int i=0; i<depthWeights.size(); i++)
					depthWeights.set(i, 1d);
				sumWeights = depthWeights.calcSumOfY_Vals();
			}
			depthWeights.scale(1d/sumWeights);
			// now change to cumulative
			double cmlSum = 0d;
			for (int i=0; i<depthWeights.size(); i++) {
				cmlSum += depthWeights.getY(i);
				depthWeights.set(i, cmlSum);
			}
			Preconditions.checkState((float)depthWeights.getY(depthWeights.size()-1) == 1f);
			targetDepthCmlWeights.add(depthWeights);
		}
		
		rand = new Random((long)rupSet.getNumRuptures()*(long)targetSects.size() + targetSects.get(0).getSectionId());
	}
	
	public ImmutableSet<Integer> getAffectedRuptures() {
		return affectedRuptures;
	}
	
	public boolean affectsRupture(int rupIndex) {
		return affectedRuptures.contains(rupIndex);
	}
	
	public double getTargetFractInPreferred() {
		return targetFractInPreferred;
	}
	
	public double getFractNuclOnTargetSects(int rupIndex) {
		double areaOn = 0d;
		double areaOff = 0d;
		
		for (int sectIndex : rupSet.getSectionsIndicesForRup(rupIndex)) {
			double sectArea = rupSet.getAreaForSection(sectIndex);
			if (targetSectIDs.contains(sectIndex))
				areaOn += sectArea;
			else
				areaOff += sectArea;
		}
		
		return areaOn/(areaOn+areaOff);
	}
	
	public List<Location> drawPreferredDirectionHypos(int rupIndex, int num) {
		return drawPreferredDirectionHypos(rupIndex, num, false);
	}
	
	public List<Location> drawPreferredDirectionHypos(int rupIndex, int num, boolean debug) {
		if (debug) System.out.println("Debugging "+name+" preferred hypocenters for rupture "+rupIndex);
		// these are indexes in our target array, not sect IDs
		List<Integer> affectedSectIndexes = new ArrayList<>();
		String affectedStr = null;
		for (int sectID : rupSet.getSectionsIndicesForRup(rupIndex)) {
			if (targetSectIDs.contains(sectID)) {
				// find the index
				int index = -1;
				for (int i=0; i<targetSects.size(); i++) {
					if (targetSects.get(i).getSectionId() == sectID) {
						index = i;
						break;
					}
				}
				Preconditions.checkState(index >= 0);
				affectedSectIndexes.add(index);
				if (debug) {
					if (affectedStr == null)
						affectedStr = "";
					else
						affectedStr += ",";
					affectedStr += sectID;
				}
			}
		}
		
		ClusterRupture cRup = cRups.get(rupIndex);
		Preconditions.checkState(!affectedSectIndexes.isEmpty(),
				"Rupture %s doesn't affect any sections for %s?\n\tRupture: %s\n\tMy Sections: %s",
				rupIndex, name, cRup, targetSectIDs);
		if (debug) System.out.println("\tRupture: "+cRup);
		if (debug) System.out.println("\tAffected sections: "+affectedStr); 
		
		List<List<Integer>> validColumns = new ArrayList<>();
		int numValidColumns = 0;
		for (int index=0; index<affectedSectIndexes.size(); index++) {
			int listIndex = affectedSectIndexes.get(index);
			FaultSection nuclSect = targetSects.get(listIndex);
			Preconditions.checkState(cRup.contains(nuclSect), "ClusterRupture for %s doesn't contain nuclSect=%s? %s",
					rupIndex, nuclSect.getSectionId(), cRup);
			EvenlyGriddedSurface nuclSurf = targetSectSurfs.get(listIndex);
			List<Integer> subValidColumns = new ArrayList<>();
			
			List<Location> jumpFromLocs = new ArrayList<>();
			for (int testIndex=0; testIndex<affectedSectIndexes.size(); testIndex++) {
				int testListIndex = affectedSectIndexes.get(testIndex);
				if (testListIndex == listIndex)
					jumpFromLocs.add(null);
				else
					jumpFromLocs.add(findJumpFromLoc(cRup, nuclSect, targetSects.get(testListIndex)));
			}
			
			if (debug) System.out.println("\tTesting nucleation from "+nuclSect.getSectionId()+". "+nuclSect.getSectionName());
			for (int col=0; col<nuclSurf.getNumCols(); col++) {
				double areaInPreferredDirection = 0d;
				double areaAgainstPreferredDirection = 0d;
				for (int testIndex=0; testIndex<affectedSectIndexes.size(); testIndex++) {
					int testListIndex = affectedSectIndexes.get(testIndex);
					FaultSection sect = targetSects.get(testListIndex);
					Preconditions.checkState(cRup.contains(nuclSect), "ClusterRupture for %s doesn't contain targetSect=%s? %s",
							rupIndex, sect.getSectionId(), cRup);
					FaultTrace trace = sect.getFaultTrace();
					double sectArea = rupSet.getAreaForSection(sect.getSectionId());
					boolean alongStrike = alongStrikes.get(testListIndex);
					if (listIndex == testListIndex) {
						// nucleates on this section
						Location surfFirstLoc = nuclSurf.get(0, 0);
						double firstLocDist = LocationUtils.horzDistanceFast(surfFirstLoc, trace.first());
						double lastLocDist = LocationUtils.horzDistanceFast(surfFirstLoc, trace.last());
						boolean surfForwards = firstLocDist < lastLocDist;
						boolean preferBefore = (!surfForwards && alongStrike) || (surfForwards && !alongStrike);
						double fractBefore = (col+1d)/(double)nuclSurf.getNumCols();
						if (preferBefore) {
							areaInPreferredDirection += sectArea*fractBefore;
							areaAgainstPreferredDirection += sectArea*(1-fractBefore);
						} else {
							areaInPreferredDirection += sectArea*(1-fractBefore);
							areaAgainstPreferredDirection += sectArea*fractBefore;
						}
					} else {
						// nucleates on another section, figure out the direction
						Location fromLoc = jumpFromLocs.get(testIndex);
						double firstLocDist = LocationUtils.horzDistanceFast(fromLoc, trace.first());
						double lastLocDist = LocationUtils.horzDistanceFast(fromLoc, trace.last());
						boolean isAlongStrike = firstLocDist < lastLocDist;
						if (isAlongStrike == alongStrike)
							// preferred direction
							areaInPreferredDirection += sectArea;
						else
							// against the grain
							areaAgainstPreferredDirection += sectArea;
					}
				}
				double fractBilateral = areaAgainstPreferredDirection/(areaAgainstPreferredDirection+areaInPreferredDirection);
				boolean valid = fractBilateral <= maxBilateralFract;
				if (debug) {
					Location loc = nuclSurf.get(0, col);
					System.out.println("\t\tColumn "+col+", fractBilateral="+(float)areaAgainstPreferredDirection
							+"/"+(float)(areaAgainstPreferredDirection+areaInPreferredDirection)+" = "+(float)fractBilateral
							+" (valid? "+valid+", surfLoc=["+(float)loc.lat+","+(float)loc.lon+"])");
				}
				if (valid) {
					// it's a match
					subValidColumns.add(col);
					numValidColumns++;
				}
			}
			validColumns.add(subValidColumns);
		}
		
		if (debug) System.out.println("\tFound "+numValidColumns+" valid columns");
//		if (numValidColumns == 0 && !debug) {
//			// do it again with debug on
//			synchronized (this) {
//				System.err.println("FAILED for "+rupIndex+", redoing");
//				drawPreferredDirectionHypos(rupIndex, num, true);
//			}
//		}
//		Preconditions.checkState(numValidColumns > 0, "No valid columns found for rupture %s: %s", rupIndex, cRup);
		if (numValidColumns == 0)
			// can happen in rare casese with weird jumps, for example between parallel strands that are both targets
			return null;
		
		List<Integer> randIndexes;
		if (num < 0) {
			// all of them
			randIndexes = new ArrayList<>(numValidColumns);
			for (int i=0; i<numValidColumns; i++)
				randIndexes.add(i);
		} else {
			// random sample
			randIndexes = new ArrayList<>(num);
			for (int i=0; i<num; i++)
				randIndexes.add(rand.nextInt(numValidColumns));
		}
		
		List<Location> ret = new ArrayList<>();
		for (int colIndex : randIndexes) {
			if (debug) System.out.println("\tDrawing random hypocenter using colIndex="+colIndex);
			// find that column
			int curIndex = 0;
			EvenlyGriddedSurface surf = null;
			EvenlyDiscretizedFunc cmlDepthWeights = null;
			int col = -1;
			for (int index=0; index<validColumns.size(); index++) {
				List<Integer> subCols = validColumns.get(index);
				for (int subCol : subCols) {
					if (curIndex == colIndex) {
						int listIndex = affectedSectIndexes.get(index);
						surf = targetSectSurfs.get(listIndex);
						cmlDepthWeights = targetDepthCmlWeights.get(listIndex);
						col = subCol;
						if (debug) System.out.println("\t\tFound it for sect "+targetSects.get(listIndex).getSectionId()+", col="+col);
						break;
					}
					curIndex++;
				}
				if (surf != null)
					break;
			}
			Preconditions.checkState(surf != null);
			
			// now need to choose a random depth
			double depth = cmlDepthWeights.getClosestXtoY(rand.nextDouble());
			int row = cmlDepthWeights.getClosestXIndex(depth);
			Location hypo = surf.get(row, col);
			if (debug) System.out.println("\t\tRand depth="+(float)depth+" is at row "+row);
			if (debug) System.out.println("\t\tHypocenter "+ret.size()+": "+hypo);
			ret.add(hypo);
		}
		
		return ret;
	}
	
	private Location findJumpFromLoc(ClusterRupture cRup, FaultSection nuclSect, FaultSection targetSect) {
		Preconditions.checkState(!nuclSect.equals(targetSect));
		RuptureTreeNavigator nav = cRup.getTreeNavigator();
		Location loc = findJumpFromLocRecursive(nav, nuclSect, targetSect, true);
		if (loc == null)
			loc = findJumpFromLocRecursive(nav, nuclSect, targetSect, false);
		Preconditions.checkNotNull(loc, "Didn't find path from nuclSect=%s to targetSect=%s for rup: %s",
				nuclSect.getSectionId(), targetSect.getSectionId(), cRup);
		return loc;
	}
	
	private Location findJumpFromLocRecursive(RuptureTreeNavigator nav, FaultSection curSect,
			FaultSection targetSect, boolean forwards) {
		Preconditions.checkState(!curSect.equals(targetSect));
		List<FaultSection> dests = new ArrayList<>();
		if (forwards) {
			for (FaultSection sect : nav.getDescendants(curSect))
				dests.add(sect);
		} else {
			FaultSection predecessor = nav.getPredecessor(curSect);
			if (predecessor != null)
				dests.add(predecessor);
		}
		Location ret = null;
		for (FaultSection dest : dests) {
			if (dest.getSectionId() == targetSect.getSectionId()) {
				// found it
				FaultTrace curTrace = curSect.getFaultTrace();
				Location curFirstLoc = curTrace.first();
				Location curLastLoc = curTrace.last();
				FaultTrace destTrace = targetSect.getFaultTrace();
				Location destFirstLoc = destTrace.first();
				Location destLastLoc = destTrace.last();
				double firstDist = Math.min(LocationUtils.horzDistanceFast(curFirstLoc, destFirstLoc),
						LocationUtils.horzDistanceFast(curFirstLoc, destLastLoc));
				double lastDist = Math.min(LocationUtils.horzDistanceFast(curLastLoc, destFirstLoc),
						LocationUtils.horzDistanceFast(curLastLoc, destLastLoc));
				if (firstDist < lastDist)
					ret = curFirstLoc;
				else
					ret = curLastLoc;
			} else {
				// keep searching
				ret = findJumpFromLocRecursive(nav, dest, targetSect, forwards);
			}
			if (ret != null)
				return ret;
		}
		return ret;
	}

	@Override
	public String getName() {
		return name;
	}
	
	public List<FaultSection> getTargetSects() {
		return ImmutableList.copyOf(targetSects);
	}
	
	public List<Boolean> getTargetAlongStrikes() {
		return ImmutableList.copyOf(alongStrikes);
	}

}
