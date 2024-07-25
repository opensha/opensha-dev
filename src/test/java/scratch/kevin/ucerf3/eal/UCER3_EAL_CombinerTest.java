package scratch.kevin.ucerf3.eal;

import static org.junit.Assert.*;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.dom4j.DocumentException;
import org.junit.BeforeClass;
import org.junit.Test;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.IDPairing;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.param.BackgroundRupParam;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.util.FocalMech;
import org.opensha.sra.calc.parallel.MPJ_CondLossCalc;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.U3CompoundFaultSystemSolution;
import scratch.UCERF3.U3FaultSystemRupSet;
import scratch.UCERF3.U3FaultSystemSolution;
import scratch.UCERF3.U3FaultSystemSolutionFetcher;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.enumTreeBranches.MaxMagOffFault;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.griddedSeismicity.AbstractGridSourceProvider;
import scratch.UCERF3.griddedSeismicity.GridSourceFileReader;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.inversion.SectionCluster;
import scratch.UCERF3.inversion.SectionClusterList;
import scratch.UCERF3.inversion.OldSectionConnectionStrategy;
import scratch.UCERF3.inversion.UCERF3SectionConnectionStrategy;
import scratch.UCERF3.inversion.laughTest.UCERF3PlausibilityConfig;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.logicTree.U3LogicTreeBranchNode;
import scratch.UCERF3.utils.DeformationModelFetcher;

public class UCER3_EAL_CombinerTest {
	
	private static InversionFaultSystemRupSet rupSet;
	private static U3FaultSystemSolutionFetcher fetch;
	private static U3FaultSystemSolution trueMeanSol;
	private static double[][] fssLosses;
	private static Map<U3LogicTreeBranch, List<Integer>> mappings;
	
	// mech => [node][mag, loss]
	private static Map<FocalMech, DiscretizedFunc[]> griddedLossFuncs;
	
	private static final int numSols = 50;
	
	private static final Random r = new Random();
	
	public static InversionFaultSystemRupSet buildSmallTestRupSet() {
		U3LogicTreeBranch branch = U3LogicTreeBranch.UCERF2;
		// this list will store our subsections
		List<FaultSection> subSections = Lists.newArrayList();
		
		FaultModels fm = branch.getValue(FaultModels.class);
		List<FaultSection> fsd = fm.getFaultSections();
		double maxSubSectionLength = 0.5;
		double maxDistance = 5d;
		
		// build the subsections
		int sectIndex = 0;
		for (FaultSection parentSect : fsd) {
			if (parentSect.getSectionId() != 301)
				// only one fault, Mojave S
				continue;
			double ddw = parentSect.getOrigDownDipWidth();
			double maxSectLength = ddw*maxSubSectionLength;
			// the "2" here sets a minimum number of sub sections
			List<? extends FaultSection> newSubSects = parentSect.getSubSectionsList(maxSectLength, sectIndex, 2);
			subSections.addAll(newSubSects);
			sectIndex += newSubSects.size();
		}
				
		UCERF3PlausibilityConfig laughTest = UCERF3PlausibilityConfig.getDefault();
		laughTest.setCoulombFilter(null);
		
		// calculate distances between each subsection
		Map<IDPairing, Double> subSectionDistances = DeformationModelFetcher.calculateDistances(maxDistance, subSections);
		Map<IDPairing, Double> reversed = Maps.newHashMap();
		// now add the reverse distance
		for (IDPairing pair : subSectionDistances.keySet()) {
			IDPairing reverse = pair.getReversed();
			reversed.put(reverse, subSectionDistances.get(pair));
		}
		subSectionDistances.putAll(reversed);
		Map<IDPairing, Double> subSectionAzimuths = DeformationModelFetcher.getSubSectionAzimuthMap(
				subSectionDistances.keySet(), subSections);
		
		OldSectionConnectionStrategy connectionStrategy = new UCERF3SectionConnectionStrategy(
				laughTest.getMaxJumpDist(), null);
		
		SectionClusterList clusters = new SectionClusterList(connectionStrategy,
				laughTest, subSections, subSectionDistances, subSectionAzimuths);
		
		List<List<Integer>> ruptures = Lists.newArrayList();
		for (SectionCluster cluster : clusters) {
			ruptures.addAll(cluster.getSectionIndicesForRuptures());
		}
		
		System.out.println("Created "+ruptures.size()+" ruptures");
		
		return new InversionFaultSystemRupSet(branch, clusters, subSections);
	}

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		rupSet = buildSmallTestRupSet();
		List<Integer> mapping = Lists.newArrayList();
		// simple 1 to 1 mapping, only mag variation
		for (int i=0; i<rupSet.getNumRuptures(); i++)
			mapping.add(i);
		mappings = Maps.newHashMap();
		
		MinMaxAveTracker latTrack = new MinMaxAveTracker();
		MinMaxAveTracker lonTrack = new MinMaxAveTracker();
		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			for (Location loc : sect.getFaultTrace()) {
				latTrack.addValue(loc.getLatitude());
				lonTrack.addValue(loc.getLongitude());
			}
		}
		GriddedRegion reg = new GriddedRegion(new Location(latTrack.getMax()+1d, lonTrack.getMax()+1d),
				new Location(latTrack.getMin()-1d, lonTrack.getMin()-1d), 0.1, null);
		Preconditions.checkState(reg.getNodeCount() > 5);
		
		HashSet<U3LogicTreeBranch> branches = new HashSet<U3LogicTreeBranch>();
		
		DiscretizedFunc[] rupMFDs = new DiscretizedFunc[rupSet.getNumRuptures()];
		for (int i=0; i<rupMFDs.length; i++)
			rupMFDs[i] = new ArbitrarilyDiscretizedFunc();
		
		Map<U3LogicTreeBranch, InversionFaultSystemSolution> map = Maps.newHashMap();
		
		Map<Integer, IncrementalMagFreqDist> meanNodeSubSeisMFDs = Maps.newHashMap();
		Map<Integer, IncrementalMagFreqDist> meanNodeUnassociatedMFDs = Maps.newHashMap();
		
		for (int i=0; i<numSols; i++) {
			// get unique branch
			U3LogicTreeBranch branch = null;
			while (branch == null || branches.contains(branch))
				branch = getRandomBranch();
			branches.add(branch);
			
			double[] mags = Arrays.copyOf(rupSet.getMagForAllRups(), rupSet.getNumRuptures());
			double[] rates = new double[mags.length];
			for (int r=0; r<mags.length; r++) {
				if (Math.random()<0.25)
					// perturb +/- 0.2
					mags[r] += (Math.random()-0.5)*0.2;
				
				rates[r] = Math.random()*1e-4;
				
				int index = rupMFDs[r].getXIndex(mags[r]);
				if (index >= 0)
					rupMFDs[r].set(index, rupMFDs[r].getY(index)+rates[r]);
				else
					rupMFDs[r].set(mags[r], rates[r]);
			}
			
			U3FaultSystemRupSet subRupSet = new U3FaultSystemRupSet(rupSet.getFaultSectionDataList(),
					rupSet.getSlipRateForAllSections(), rupSet.getSlipRateStdDevForAllSections(), rupSet.getAreaForAllSections(),
					rupSet.getSectionIndicesForAllRups(), mags, rupSet.getAveRakeForAllRups(), rupSet.getAreaForAllRups(),
					rupSet.getLengthForAllRups(), "");
			
			InversionFaultSystemRupSet subInvRupSet = new InversionFaultSystemRupSet(subRupSet, branch, rupSet.getOldPlausibilityConfiguration(),
					rupSet.getAveSlipForAllRups(), rupSet.getCloseSectionsListList(),
					rupSet.getRupturesForClusters(), rupSet.getSectionsForClusters());
			
			InversionFaultSystemSolution sol = new InversionFaultSystemSolution(subInvRupSet, rates);
			
			GridSourceProvider gridProv = buildRandomGridProv(reg, branch.getValue(MaxMagOffFault.class));
			sol.setGridSourceProvider(gridProv);
			
			int numAboveZero = 0;
			
			for (int n=0; n<reg.getNodeCount(); n++) {
				IncrementalMagFreqDist unMFD = gridProv.getMFD_Unassociated(n);
				IncrementalMagFreqDist subMFD = gridProv.getMFD_SubSeisOnFault(n);
				
				if (unMFD != null) {
					IncrementalMagFreqDist meanUnMFD = meanNodeUnassociatedMFDs.get(n);
					if (meanUnMFD == null) {
						meanUnMFD = new SummedMagFreqDist(mfdMin, mfdMax, mfdNum);
						meanNodeUnassociatedMFDs.put(n, meanUnMFD);
					}
					Preconditions.checkState(unMFD.calcSumOfY_Vals() > 0);
					((SummedMagFreqDist)meanUnMFD).addIncrementalMagFreqDist(unMFD);
				}
				
				if (subMFD != null) {
					IncrementalMagFreqDist meanSubMFD = meanNodeSubSeisMFDs.get(n);
					if (meanSubMFD == null) {
						meanSubMFD = new SummedMagFreqDist(mfdMin, mfdMax, mfdNum);
						meanNodeSubSeisMFDs.put(n, meanSubMFD);
					}
					Preconditions.checkState(subMFD.calcSumOfY_Vals() > 0);
					((SummedMagFreqDist)meanSubMFD).addIncrementalMagFreqDist(subMFD);
				}
				
				if (gridProv.getMFD(n, AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF).calcSumOfY_Vals()>0)
					numAboveZero++;
				Preconditions.checkState(gridProv.getMFD(n, AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF).size() > 0);
			}
			Preconditions.checkState(numAboveZero>0, "Sol "+i+" has all zero mfd nodes!");
			
			map.put(branch, sol);
			
			mappings.put(branch, mapping);
		}
		
		Preconditions.checkState(mappings.size() == numSols);
		
		fetch = new InMemFSSFetch(map);
		double[] rates = new double[rupMFDs.length];
		for (int r=0; r<rupMFDs.length; r++) {
			double sumY = 0;
			for (Point2D pt : rupMFDs[r])
				sumY += pt.getY();
			rates[r] = sumY/(double)numSols;
		}
		trueMeanSol = new U3FaultSystemSolution(rupSet, rates);
		trueMeanSol.setRupMagDists(rupMFDs);
		
		GridSourceFileReader meanGridProv = new GridSourceFileReader(reg, meanNodeSubSeisMFDs, meanNodeUnassociatedMFDs);
		trueMeanSol.setGridSourceProvider(meanGridProv);
		// now create losses for each rup type
		griddedLossFuncs = Maps.newHashMap();
		griddedLossFuncs.put(FocalMech.NORMAL, new DiscretizedFunc[meanGridProv.size()]);
		griddedLossFuncs.put(FocalMech.REVERSE, new DiscretizedFunc[meanGridProv.size()]);
		griddedLossFuncs.put(FocalMech.STRIKE_SLIP, new DiscretizedFunc[meanGridProv.size()]);
		for (int n=0; n<meanGridProv.size(); n++) {
			// make sure not nan
			for (Point2D pt : meanGridProv.getMFD(n))
				Preconditions.checkState(!Double.isNaN(pt.getY()));
			griddedLossFuncs.get(FocalMech.NORMAL)[n] = getRandMFD(mfdMax, 100d);
			griddedLossFuncs.get(FocalMech.REVERSE)[n] = getRandMFD(mfdMax, 100d);
			griddedLossFuncs.get(FocalMech.STRIKE_SLIP)[n] = getRandMFD(mfdMax, 100d);
		}
		
		// now calculate fake fault EALs
		// choose one random empty node
		int randEmptyRup = r.nextInt(rupMFDs.length);
		fssLosses = new double[rupSet.getNumRuptures()][];
		for (int i=0; i<fssLosses.length; i++) {
			if (i == randEmptyRup)
				fssLosses[i] = new double[0];
			else
				fssLosses[i] = new double[rupMFDs[i].size()];
			for (int j=0; j<fssLosses[i].length; j++)
				fssLosses[i][j] = r.nextDouble() * 1000d;
		}
	}
	
	private static double mfdMin = 5.05;
	private static double mfdMax = 8.45;
	private static int mfdNum = 35;
	
	public static GridSourceProvider buildRandomGridProv(GriddedRegion reg, MaxMagOffFault mmaxOff) {
		Map<Integer, IncrementalMagFreqDist> nodeSubSeisMFDs = Maps.newHashMap();
		Map<Integer, IncrementalMagFreqDist> nodeUnassociatedMFDs = Maps.newHashMap();
		
		double mmax = mmaxOff.getMaxMagOffFault();
		
		for (int i=0; i<reg.getNodeCount(); i++) {
			if (r.nextDouble() < 0.5)
				nodeUnassociatedMFDs.put(i, getRandMFD(mmax, 1e-6));
			if (r.nextDouble() < 0.1)
				nodeSubSeisMFDs.put(i, getRandMFD(mmax, 1e-6));
		}
		
		return new GridSourceFileReader(reg, nodeSubSeisMFDs, nodeUnassociatedMFDs);
	}
	
	private static IncrementalMagFreqDist getRandMFD(double mmax, double scale) {
		IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(mfdMin, mfdMax, mfdNum);
		for (int j=0; j<mfdNum; j++) {
			if (mfd.getX(j)>mmax)
				break;
			mfd.set(j, scale*r.nextDouble());
		}
		return mfd;
	}
	
	private static class InMemFSSFetch extends U3FaultSystemSolutionFetcher {
		
		private Map<U3LogicTreeBranch, InversionFaultSystemSolution> map;
		
		public InMemFSSFetch(Map<U3LogicTreeBranch, InversionFaultSystemSolution> map) {
			this.map = map;
		}

		@Override
		public Collection<U3LogicTreeBranch> getBranches() {
			return map.keySet();
		}

		@Override
		protected InversionFaultSystemSolution fetchSolution(
				U3LogicTreeBranch branch) {
			return map.get(branch);
		}
		
	}
	
	private static U3LogicTreeBranch getRandomBranch() {
		List<U3LogicTreeBranchNode<?>> vals = Lists.newArrayList();
		for (Class<? extends U3LogicTreeBranchNode<?>> clazz : U3LogicTreeBranch.getLogicTreeNodeClasses()) {
			if (clazz.equals(FaultModels.class))
				vals.add(FaultModels.FM3_1);
			else
				vals.add(getRandomElem(clazz));
		}
		return U3LogicTreeBranch.fromValues(vals);
	}
	
	private static U3LogicTreeBranchNode<?> getRandomElem(Class<? extends U3LogicTreeBranchNode<?>> clazz) {
		List<U3LogicTreeBranchNode<?>> vals = Lists.newArrayList();
		for (U3LogicTreeBranchNode<?> val : clazz.getEnumConstants()) {
			if (val.getRelativeWeight(InversionModels.CHAR_CONSTRAINED) > 0)
				vals.add(val);
		}
		return vals.get(r.nextInt(vals.size()));
	}

	@Test
	public void testFaultBased() throws DocumentException, IOException {
		UCERF3_EAL_Combiner comb = new UCERF3_EAL_Combiner(fetch, mappings, trueMeanSol, fssLosses, null);
		
		DiscretizedFunc[] rupMFDs = trueMeanSol.getRupMagDists();
		
		double[] faultEALs = comb.getFaultEALs();
		List<U3LogicTreeBranch> branches = comb.getBranches();
		
		assertEquals(branches.size(), numSols);
		assertEquals(branches.size(), faultEALs.length);
		
		for (int i=0; i<numSols; i++) {
			InversionFaultSystemSolution sol = fetch.getSolution(branches.get(i));
			// calculate my own fault eal
			U3FaultSystemRupSet rupSet = sol.getRupSet();
			double eal = 0d;
			for (int r=0; r<rupSet.getNumRuptures(); r++) {
				if (fssLosses[r].length == 0)
					continue;
				double mag = rupSet.getMagForRup(r);
				int mfdIndex = rupMFDs[r].getXIndex(mag);
				double loss = fssLosses[r][mfdIndex];
				eal += loss * sol.getRateForRup(r);
			}
			
			double calcEAL = faultEALs[i];
			double pDiff = DataUtils.getPercentDiff(calcEAL, eal);
			assertEquals(0, pDiff, 1e-4);
		}
	}
	
	private static double lossForGridRup(ProbEqkRupture rup, int node) {
		FocalMech mech;
		if ((float)rup.getAveRake() == 90f)
			mech = FocalMech.REVERSE;
		else if ((float)rup.getAveRake() == -90f)
			mech = FocalMech.NORMAL;
		else
			mech = FocalMech.STRIKE_SLIP;
		double mag = rup.getMag();
		return griddedLossFuncs.get(mech)[node].getY(mag);
	}

	public void doTestGridded(BackgroundRupType bgType) throws DocumentException, IOException {
		GridSourceProvider meanProv = trueMeanSol.getGridSourceProvider();
		GriddedRegion reg = meanProv.getGriddedRegion();
		// first create random data
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(trueMeanSol);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.ONLY);
		erf.setParameter(BackgroundRupParam.NAME, bgType);
		erf.updateForecast();
		Preconditions.checkState(erf.getNumSources() == reg.getNodeCount());
		double[][] origLosses = new double[erf.getNumSources()][];
		Preconditions.checkState(erf.getNumSources() == reg.getNodeCount());
		for (int n=0; n<origLosses.length; n++) {
			ProbEqkSource src = erf.getSource(n);
			origLosses[n] = new double[src.getNumRuptures()];
			for (int r=0; r<src.getNumRuptures(); r++) {
				ProbEqkRupture rup = src.getRupture(r);
				origLosses[n][r] = lossForGridRup(rup, n);
			}
		}
		File tempFile = File.createTempFile("openSHA", "test.bin");
		MPJ_CondLossCalc.writeFSSGridSourcesFile(erf, origLosses, tempFile);
		DiscretizedFunc[] griddedLosses = MPJ_CondLossCalc.loadGridSourcesFile(tempFile, reg);
		tempFile.delete();
		UCERF3_EAL_Combiner comb = new UCERF3_EAL_Combiner(fetch, mappings, trueMeanSol, fssLosses, griddedLosses);
		
		double[] griddedEALs = comb.getGriddedEALs();
		List<U3LogicTreeBranch> branches = comb.getBranches();
		
		assertEquals(branches.size(), numSols);
		assertEquals(branches.size(), griddedEALs.length);
		
		for (int i=0; i<numSols; i++) {
			InversionFaultSystemSolution sol = fetch.getSolution(branches.get(i));
			// calculate my own gridded eal
			GridSourceProvider gridProv = sol.getGridSourceProvider();
			double eal = 0;
			int numNonZeroRate = 0;
			int numNonZeroLoss = 0;
			int numMFDNonZero = 0;
			for (int n=0; n<reg.getNodeCount(); n++) {
				ProbEqkSource src = gridProv.getSource(n, 1d, null, bgType);
				// make sure not nan
				for (Point2D pt : gridProv.getMFD(n, AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF))
					Preconditions.checkState(!Double.isNaN(pt.getY()));
				if (gridProv.getMFD(n, AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF).calcSumOfY_Vals() > 0)
					numMFDNonZero++;
				for (ProbEqkRupture rup : src) {
					double loss = lossForGridRup(rup, n);
					double rate = rup.getMeanAnnualRate(1d);
					if (rate > 0)
						numNonZeroRate++;
					if (loss > 0)
						numNonZeroLoss++;
					Preconditions.checkState(!Double.isNaN(loss));
					Preconditions.checkState(!Double.isNaN(rate), "NaN rate. Prob was "+rup.getProbability());
					eal += loss*rate;
				}
			}
			
			double calcEAL = griddedEALs[i];
			
			System.out.println("Calc: "+calcEAL);
			System.out.println("Actual: "+eal);
			System.out.println("Num Non Zero: rate="+numNonZeroRate+" loss="+numNonZeroLoss+" mfd="+numMFDNonZero);
			
			double pDiff = DataUtils.getPercentDiff(calcEAL, eal);
			assertEquals(0, pDiff, 1e-4);
		}
	}
	
	@Test
	public void testGriddedCrosshair() throws DocumentException, IOException {
		doTestGridded(BackgroundRupType.CROSSHAIR);
	}
	
	@Test
	public void testGriddedFininte() throws DocumentException, IOException {
		doTestGridded(BackgroundRupType.FINITE);
	}
	
	@Test
	public void testGriddedPoint() throws DocumentException, IOException {
		doTestGridded(BackgroundRupType.POINT);
	}

}
