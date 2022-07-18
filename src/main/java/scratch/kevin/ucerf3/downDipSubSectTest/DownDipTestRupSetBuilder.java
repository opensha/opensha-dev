package scratch.kevin.ucerf3.downDipSubSectTest;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.geo.Location;
import org.opensha.commons.util.FaultUtils;
import org.opensha.commons.util.IDPairing;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionInputGenerator;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SlipRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ConstraintRange;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ThreadedSimulatedAnnealing;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.CompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.ProgressTrackingCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.TimeCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRuptureBuilder;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.ClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.RuptureGrowingStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.DistCutoffClosestSectClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.SimpleFaultData;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import scratch.UCERF3.U3FaultSystemRupSet;
import scratch.UCERF3.U3FaultSystemSolution;
import scratch.UCERF3.U3SlipAlongRuptureModelRupSet;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;
import scratch.UCERF3.inversion.SectionCluster;
import scratch.UCERF3.inversion.SectionClusterList;
import scratch.UCERF3.inversion.OldSectionConnectionStrategy;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.U3FaultSystemIO;
import scratch.UCERF3.utils.MFD_InversionConstraint;

public class DownDipTestRupSetBuilder {

	public static void main(String[] args) throws IOException {
		boolean runInversion = true;
		long inversionMins = 1; // run it for this many minutes
		
		File rupSetFile = new File("/tmp/down_dip_sub_sect_rup_set.zip");
		File solFile = new File("/tmp/down_dip_sub_sect_sol.zip");
		
		// we're going to manually build faults here
		// first, build a big fault with down-dip subsections
		String sectName = "Test SubSect Down-Dip Fault";
		int sectID = 0;
		int startID = 0;
		double upperDepth = 0d;
		double lowerDepth = 30d;
		double dip = 20d;
//		int numDownDip = 5;
//		int numAlongStrike = 6;
		int numDownDip = 20;
		int numAlongStrike = 30;
		FaultTrace trace = new FaultTrace(sectName);
		trace.add(new Location(34, -119, upperDepth));
		trace.add(new Location(34.1, -118.75, upperDepth));
		trace.add(new Location(34.15, -118.5, upperDepth));
		trace.add(new Location(34.1, -118.25, upperDepth));
		trace.add(new Location(34, -118, upperDepth));
		double slipRate = 50d; // mm/yr
		
		SimpleFaultData faultData = new SimpleFaultData(dip, lowerDepth, upperDepth, trace);
		double aveRake = 90d;
		
		DownDipSubSectBuilder downDipBuilder = new DownDipSubSectBuilder(sectName, sectID, startID,
				faultData, aveRake, numAlongStrike, numDownDip, slipRate);
		
		List<FaultSection> subSections = new ArrayList<>();
		subSections.addAll(downDipBuilder.getSubSectsList());
		System.out.println("Have "+subSections.size()+" sub-sections for "+sectName);
		startID = subSections.size();
		
		// now lets add a nearby crustal fault that it can interact with
		sectName = "Test Crustal Fault";
		sectID = 1;
		upperDepth = 0d;
		lowerDepth = 14d;
		dip = 90d;
		slipRate = 20d; // mm/yr
		trace = new FaultTrace(sectName);
		trace.add(new Location(33.5, -117.5, upperDepth));
		trace.add(new Location(34.1, -118.25, upperDepth));
		FaultSectionPrefData crustalSect = new FaultSectionPrefData();
		crustalSect.setSectionId(sectID);
		crustalSect.setSectionName(sectName);
		crustalSect.setFaultTrace(trace);
		crustalSect.setAveUpperDepth(upperDepth);
		crustalSect.setAveLowerDepth(lowerDepth);
		crustalSect.setAseismicSlipFactor(0d);
		crustalSect.setAveSlipRate(slipRate);
		crustalSect.setAveDip(dip);
		double maxSectLength = 0.5*crustalSect.getOrigDownDipWidth();
		
		subSections.addAll(crustalSect.getSubSectionsList(maxSectLength, startID, 2));
		System.out.println("Have "+subSections.size()+" sub-sections in total");
		
		for (int s=0; s<subSections.size(); s++)
			Preconditions.checkState(subSections.get(s).getSectionId() == s,
				"section at index %s has ID %s", s, subSections.get(s).getSectionId());
		
		// instantiate plausibility filters
		List<PlausibilityFilter> filters = new ArrayList<>();
		int minDimension = 1; // minimum numer of rows or columns
		double maxAspectRatio = 5d; // max aspect ratio of rows/cols or cols/rows
		filters.add(new RectangularityFilter(downDipBuilder, minDimension, maxAspectRatio));
		
		SectionDistanceAzimuthCalculator distAzCalc = new SectionDistanceAzimuthCalculator(subSections);
		
		// this creates rectangular permutations only for our down-dip fault to speed up rupture building
		RuptureGrowingStrategy permutationStrategy = new DownDipTestPermutationStrategy(downDipBuilder);
		// connection strategy: parent faults connect at closest point, and only when dist <=5 km
		ClusterConnectionStrategy connectionStrategy = new DistCutoffClosestSectClusterConnectionStrategy(
				subSections, distAzCalc, 5d);
		int maxNumSplays = 0; // don't allow any splays
		
//		ClusterRuptureBuilder builder = new ClusterRuptureBuilder(subSections, connectionStrategy,
//				distAzCalc, filters, maxNumSplays);
		ClusterRuptureBuilder builder = new ClusterRuptureBuilder(
				connectionStrategy.getClusters(), filters, maxNumSplays, distAzCalc);
		
		List<ClusterRupture> ruptures = builder.build(permutationStrategy);
		
		System.out.println("Built "+ruptures.size()+" total ruptures");
		
		MySlipEnabledRupSet rupSet = new MySlipEnabledRupSet(ruptures, subSections,
				ScalingRelationships.SHAW_2009_MOD, SlipAlongRuptureModels.UNIFORM);
		
		U3FaultSystemIO.writeRupSet(rupSet, rupSetFile);
		
		if (runInversion) {
			List<InversionConstraint> constraints = new ArrayList<>();
			
			/*
			 * Slip rate constraints
			 */
			double slipRateConstraintWt = 1;
			ConstraintWeightingType slipRateWeighting = ConstraintWeightingType.UNNORMALIZED;
			constraints.add(new SlipRateInversionConstraint(slipRateConstraintWt, slipRateWeighting, rupSet));
			
			/*
			 * MFD constraints
			 */
			double totalRateM5 = 5d; // expected number of M>=5's per year
			double bValue = 1d; // G-R b-value
			// magnitude to switch from MFD equality to MFD inequality
			double mfdTransitionMag = 7.85;
			double mfdEqualityConstraintWt = 10;
			double mfdInequalityConstraintWt = 1000;
			int mfdNum = 40;
			double mfdMin = 5.05d;
			double mfdMax = 8.95;
			GutenbergRichterMagFreqDist mfd = new GutenbergRichterMagFreqDist(
					bValue, totalRateM5, mfdMin, mfdMax, mfdNum);
			int transitionIndex = mfd.getClosestXIndex(mfdTransitionMag);
			// snap it to the discretization if it wasn't already
			mfdTransitionMag = mfd.getX(transitionIndex);
			Preconditions.checkState(transitionIndex >= 0);
			GutenbergRichterMagFreqDist equalityMFD = new GutenbergRichterMagFreqDist(
					bValue, totalRateM5, mfdMin, mfdTransitionMag, transitionIndex);
			MFD_InversionConstraint equalityConstr = new MFD_InversionConstraint(equalityMFD, null);
			GutenbergRichterMagFreqDist inequalityMFD = new GutenbergRichterMagFreqDist(
					bValue, totalRateM5, mfdTransitionMag, mfdMax, mfd.size()-equalityMFD.size());
			MFD_InversionConstraint inequalityConstr = new MFD_InversionConstraint(inequalityMFD, null);

//			constraints.add(new MFDInversionConstraint(rupSet, mfdEqualityConstraintWt, false,
//					Lists.newArrayList(equalityConstr), null));
//			constraints.add(new MFDInversionConstraint(rupSet, mfdInequalityConstraintWt, true,
//					Lists.newArrayList(inequalityConstr), null));
			
			// weight of entropy-maximization constraint (not used in UCERF3)
			double smoothnessWt = 0;
			
			/*
			 * Build inversion inputs
			 */
			InversionInputGenerator inputGen = new InversionInputGenerator(rupSet, constraints);
			
			inputGen.generateInputs(true);
			// column compress it for fast annealing
			inputGen.columnCompress();
			
			// inversion completion criteria (how long it will run)
			CompletionCriteria criteria = TimeCompletionCriteria.getInMinutes(inversionMins);
			
			// Bring up window to track progress
			criteria = new ProgressTrackingCompletionCriteria(criteria, 0.1);

			// this will use all available processors
			int numThreads = Runtime.getRuntime().availableProcessors();

			// this is the "sub completion criteria" - the amount of time (or iterations) between synchronization
			CompletionCriteria subCompetionCriteria = TimeCompletionCriteria.getInSeconds(1); // 1 second;
			
			ThreadedSimulatedAnnealing tsa = new ThreadedSimulatedAnnealing(inputGen.getA(), inputGen.getD(),
					inputGen.getInitialSolution(), smoothnessWt, inputGen.getA_ineq(), inputGen.getD_ineq(),
					numThreads, subCompetionCriteria);
			tsa.setConstraintRanges(inputGen.getConstraintRowRanges());
			
			tsa.iterate(criteria);
			
			// now assemble the solution
			double[] solution_raw = tsa.getBestSolution();
			
			// adjust for minimum rates if applicable
			double[] solution_adjusted = inputGen.adjustSolutionForWaterLevel(solution_raw);
			
			Map<ConstraintRange, Double> energies = tsa.getEnergies();
			if (energies != null) {
				System.out.println("Final energies:");
				for (ConstraintRange range : energies.keySet())
					System.out.println("\t"+range.name+": "+energies.get(range).floatValue());
			}
			
			// now write out the solution
			U3FaultSystemSolution sol = new U3FaultSystemSolution(rupSet, solution_adjusted);
			U3FaultSystemIO.writeSol(sol, solFile);
		}
		
		System.exit(0);
	}
	
	private static class MySlipEnabledRupSet extends U3SlipAlongRuptureModelRupSet {
		
		private double[] rupAveSlips;
		
		public MySlipEnabledRupSet(List<ClusterRupture> ruptures, List<FaultSection> subSections,
				ScalingRelationships scale, SlipAlongRuptureModels slipAlongModel) {
			super(slipAlongModel);
			
			// build a rupture set (doing this manually instead of creating an inversion fault system rup set,
			// mostly as a demonstration)
			double[] sectSlipRates = new double[subSections.size()];
			double[] sectAreasReduced = new double[subSections.size()];
			double[] sectAreasOrig = new double[subSections.size()];
			for (int s=0; s<sectSlipRates.length; s++) {
				FaultSection sect = subSections.get(s);
				sectAreasReduced[s] = sect.getArea(true);
				sectAreasOrig[s] = sect.getArea(false);
				sectSlipRates[s] = sect.getReducedAveSlipRate()*1e-3; // mm/yr => m/yr
			}
			
			double[] rupMags = new double[ruptures.size()];
			double[] rupRakes = new double[ruptures.size()];
			double[] rupAreas = new double[ruptures.size()];
			double[] rupLengths = new double[ruptures.size()];
			rupAveSlips = new double[ruptures.size()];
			List<List<Integer>> rupsIDsList = new ArrayList<>();
			for (int r=0; r<ruptures.size(); r++) {
				ClusterRupture rup = ruptures.get(r);
				List<FaultSection> rupSects = rup.buildOrderedSectionList();
				List<Integer> sectIDs = new ArrayList<>();
				double totLength = 0d;
				double totArea = 0d;
				double totOrigArea = 0d; // not reduced for aseismicity
				List<Double> sectAreas = new ArrayList<>();
				List<Double> sectRakes = new ArrayList<>();
				for (FaultSection sect : rupSects) {
					sectIDs.add(sect.getSectionId());
					double length = sect.getTraceLength()*1e3;	// km --> m
					totLength += length;
					double area = sectAreasReduced[sect.getSectionId()];	// sq-m
					totArea += area;
					totOrigArea += sectAreasOrig[sect.getSectionId()];	// sq-m
					sectAreas.add(area);
					sectRakes.add(sect.getAveRake());
				}
				rupAreas[r] = totArea;
				rupLengths[r] = totLength;
				rupRakes[r] = FaultUtils.getInRakeRange(FaultUtils.getScaledAngleAverage(sectAreas, sectRakes));
				double origDDW = totOrigArea/totLength;
				rupMags[r] = scale.getMag(totArea, totLength, totArea/totLength, origDDW, rupRakes[r]);
				rupsIDsList.add(sectIDs);
				rupAveSlips[r] = scale.getAveSlip(totArea, totLength, totArea/totLength, origDDW, rupRakes[r]);
			}
			
			String info = "Test down-dip subsectioning rup set";
			
			init(subSections, sectSlipRates, null, sectAreasReduced,
					rupsIDsList, rupMags, rupRakes, rupAreas, rupLengths, info);
		}

		@Override
		public double getAveSlipForRup(int rupIndex) {
			return rupAveSlips[rupIndex];
		}

		@Override
		public double[] getAveSlipForAllRups() {
			return rupAveSlips;
		}
		
	}

}
