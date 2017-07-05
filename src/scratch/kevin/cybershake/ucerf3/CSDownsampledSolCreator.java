package scratch.kevin.cybershake.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.hazardMap.HazardCurveSetCalculator;
import org.opensha.sha.cybershake.db.CybershakeSite;
import org.opensha.sha.cybershake.db.CybershakeSiteInfo2DB;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.plot.HazardCurvePlotter;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenuationRelationship;
import org.opensha.sha.imr.attenRelImpl.NGAWest_2014_Averaged_AttenRel;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.SiteTranslator;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.FaultSystemIO;

public class CSDownsampledSolCreator {
	
	private enum RupInRegionStatus {
		FULLY_INSIDE,
		PARTIAL,
		FULLY_OUTSIDE
	}
	
	/**
	 * This will take a FaultSystemSolution and a region and genrate a new FSS that:<br>
	 * <br>
	 * 1) Removes all ruptures that are entirely outside of the region<br>
	 * 2) Removes all ruptures that are partially outside of the region, reassigning their rate
	 * to the largest subset rupture which is entirely contained inside the region.
	 * @param inputSol
	 * @param region
	 * @return
	 */
	public static FaultSystemSolution getDownsampledSol(FaultSystemSolution inputSol, Region region) {
		return getDownsampledSol(inputSol, region, false);
	}
	
	public static FaultSystemSolution getDownsampledSol(FaultSystemSolution inputSol, Region region, boolean excludePartial) {
		FaultSystemRupSet inputRupSet = inputSol.getRupSet();
		
		// set of all sections inside
		HashSet<Integer> sectsInsideSet = new HashSet<Integer>();
		for (int s=0; s<inputRupSet.getNumSections(); s++) {
			// count as inside if any point in the subsection is inside the region
			for (Location loc : inputRupSet.getFaultSectionData(s).getFaultTrace()) {
				if (region.contains(loc)) {
					sectsInsideSet.add(s);
					break;
				}
			}
		}
		
		// list of status for each rupture to be used later
		List<RupInRegionStatus> rupStatuses = Lists.newArrayList();
		// mapping from sections to ID for all ruptures that are fully inside
		Map<String, Integer> insideRupSectsToIndexMap = Maps.newHashMap();
		int numRupsInside = 0;
		
		for (int r=0; r<inputRupSet.getNumRuptures(); r++) {
			boolean partiallyInside = false;
			boolean fullyInside = true;
			
			for (int s : inputRupSet.getSectionsIndicesForRup(r)) {
				boolean inside = sectsInsideSet.contains(s);
				partiallyInside = partiallyInside || inside;
				fullyInside = fullyInside && inside;
			}
			
			RupInRegionStatus status;
			if (fullyInside)
				status = RupInRegionStatus.FULLY_INSIDE;
			else if (partiallyInside && !excludePartial)
				status = RupInRegionStatus.PARTIAL;
			else
				status = RupInRegionStatus.FULLY_OUTSIDE;
			rupStatuses.add(status);
			
			if (status == RupInRegionStatus.FULLY_INSIDE) {
				insideRupSectsToIndexMap.put(sectsHash(inputRupSet.getSectionsIndicesForRup(r)), r);
				numRupsInside++;
			}
		}
		
		// map from rupture ID kept to list of rupture IDs whose rates will be added to it
		Map<Integer, List<Integer>> rupAssignments = Maps.newHashMap();
		// set of ruptures to be removed
		HashSet<Integer> rupsRemoved = new HashSet<Integer>();
		
		int numPartials = 0;
		int numComplicatedPartials = 0;
		
		// now do assignments
		for (int r=0; r<inputRupSet.getNumRuptures(); r++) {
			RupInRegionStatus status = rupStatuses.get(r);
			
			if (status == RupInRegionStatus.FULLY_INSIDE)
				// easy case
				continue;
			if (status == RupInRegionStatus.FULLY_OUTSIDE) {
				// also easy
				rupsRemoved.add(r);
				continue;
			}
			numPartials++;
			// partial case. trickier
			rupsRemoved.add(r);
			
			// now we must find the largest subset which is fully inside
			List<Integer> fsdIndexes = inputRupSet.getSectionsIndicesForRup(r);
			List<Integer> indexesInside = Lists.newArrayList();
			for (int s : fsdIndexes)
				if (sectsInsideSet.contains(s))
					indexesInside.add(s);
			if (indexesInside.size() < 2) {
				// not enough sections inside for a real rupture, don't reassign
				continue;
			}
			Integer mappedRupIndex = insideRupSectsToIndexMap.get(sectsHash(indexesInside));
			if (mappedRupIndex == null) {
				numComplicatedPartials++;
				int sectsToRemove = 1;
				while (mappedRupIndex == null) {
					// more complicated case
					System.out.println("Complicated case for mapping with "
							+indexesInside.size()+"/"+fsdIndexes.size()+" sects inside");
					// we should be able to remove 1 section from either end and get a real rupture
					if (indexesInside.size() < (2+sectsToRemove)) {
						System.out.println("Resolved by removal, not enough sections for complicated solution");
						break;
					}
					List<Integer> newIndexesInside = indexesInside.subList(sectsToRemove, indexesInside.size());
					mappedRupIndex = insideRupSectsToIndexMap.get(sectsHash(newIndexesInside));
					if (mappedRupIndex != null) {
						System.out.println("Resolved by removing "+sectsToRemove+" from start");
					} else {
						newIndexesInside = indexesInside.subList(0, indexesInside.size()-sectsToRemove);
						mappedRupIndex = insideRupSectsToIndexMap.get(sectsHash(newIndexesInside));
						if (mappedRupIndex == null) {
							System.out.println("Couldn't resolve by removing "+sectsToRemove+" from either end, trying again");
							sectsToRemove++;
//							String str = "Couldn't resolve by removing a sect from either side!";
//							str += "\n\tOrig sects: "+Joiner.on(",").join(fsdIndexes);
//							str += "\n\tSects inside: "+Joiner.on(",").join(indexesInside);
//							throw new IllegalStateException(str);
						}
						System.out.println("Resolved by removing "+sectsToRemove+" from end");
					}
				}
			}
			List<Integer> currentAssignments = rupAssignments.get(mappedRupIndex);
			if (currentAssignments == null) {
				currentAssignments = Lists.newArrayList();
				rupAssignments.put(mappedRupIndex, currentAssignments);
			}
			currentAssignments.add(r);
		}
		
		String info = "Downsampled FSS for "+region.getName()+" with "
				+numRupsInside+"/"+inputRupSet.getNumRuptures()+" ruptures retained. "
				+numPartials+" partial ruptures, "+numComplicatedPartials+" of which were complicated.";
		System.out.println(info);
		info += "\n\nOriginal metadata:\n"+inputRupSet.getInfoString();
		
		List<List<Integer>> sectionForRups = Lists.newArrayList();
		double[] mags = new double[numRupsInside];
		double[] rakes = new double[numRupsInside];
		double[] rupAreas = new double[numRupsInside];
//		double[] rupLengths = new double[numRupsInside];
		double[] rupLengths;
		if (inputRupSet.getLengthForAllRups() == null)
			rupLengths = null;
		else
			rupLengths = new double[numRupsInside];
		double[] rates = new double[numRupsInside];
		
		int numAbove6p5 = 0;
		
		int newIndex = 0;
		for (int r=0; r<inputRupSet.getNumRuptures(); r++) {
			if (rupStatuses.get(r) != RupInRegionStatus.FULLY_INSIDE)
				continue;
			// we're fully inside
			sectionForRups.add(inputRupSet.getSectionsIndicesForRup(r));
			mags[newIndex] = inputRupSet.getMagForRup(r);
			rakes[newIndex] = inputRupSet.getAveRakeForRup(r);
			rupAreas[newIndex] = inputRupSet.getAreaForRup(r);
			if (rupLengths != null)
				rupLengths[newIndex] = inputRupSet.getLengthForRup(r);
			
			double rate = inputSol.getRateForRup(r);
			List<Integer> mappings = rupAssignments.get(r);
			if (mappings != null) {
				for (int mapping : mappings)
					rate += inputSol.getRateForRup(mapping);
			}
			rates[newIndex] = rate;
			
			if (mags[newIndex] >= 6.5)
				numAbove6p5++;
			
			newIndex++;
		}
		Preconditions.checkState(newIndex == numRupsInside, "Indexing mismatch: "+newIndex+" != "+numRupsInside);
		System.out.println(numAbove6p5+" rups with M>=6.5");
		
		FaultSystemRupSet sampledRupSet = new FaultSystemRupSet(inputRupSet.getFaultSectionDataList(),
				inputRupSet.getSlipRateForAllSections(), inputRupSet.getSlipRateStdDevForAllSections(),
				inputRupSet.getAreaForAllSections(), sectionForRups, mags, rakes, rupAreas, rupLengths, info);
		
		FaultSystemSolution sol = new FaultSystemSolution(sampledRupSet, rates);
		
		if (plotDir != null && !excludePartial) {
			try {
				FaultSystemIO.writeSol(sol, new File(plotDir, "downdampled_sol.zip"));
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		return sol;
	}
	
	private static Joiner j = Joiner.on(",");
	
	private static String sectsHash(Collection<Integer> sections) {
		List<Integer> sorted = Lists.newArrayList(sections);
		Collections.sort(sorted);
		return j.join(sorted);
	}
	
	static void plotMFD(Region region, FaultSystemSolution orig, FaultSystemSolution downsampled,
			FaultSystemSolution downsampledExclude) {
		IncrementalMagFreqDist origMFD = orig.calcNucleationMFD_forRegion(region, 5.05, 8.55, 0.1, true);
		origMFD.setName("Original");
		IncrementalMagFreqDist excludeMFD = null;
		if (downsampledExclude != null) {
			excludeMFD = downsampledExclude.calcNucleationMFD_forRegion(region, 5.05, 8.55, 0.1, true);
			excludeMFD.setName("Excluding");
		}
		IncrementalMagFreqDist newMFD = downsampled.calcNucleationMFD_forRegion(region, 5.05, 8.55, 0.1, true);
		newMFD.setName("Downsampled");
		
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		funcs.add(origMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		if (excludeMFD != null) {
			funcs.add(excludeMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GREEN.darker()));
		}
		funcs.add(newMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Downsampled Nucleation MFD Comparison", "Mag", "Incremental Rate (1/yr)");
		spec.setLegendVisible(true);
		GraphWindow gw = new GraphWindow(spec);
		gw.getGraphWidget().setBackgroundColor(Color.WHITE);
		gw.setYLog(true);
		gw.setY_AxisRange(1e-5, 1e-1);
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		if (plotDir != null) {
			try {
				gw.saveAsPNG(new File(plotDir, "mfd_comparison.png").getAbsolutePath());
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		// now plot Mojave MFD
		int parentSect = 301;
		
		origMFD = orig.calcParticipationMFD_forParentSect(parentSect, 5.05, 8.55, origMFD.size());
		origMFD.setName("Original");
		if (downsampledExclude != null) {
			excludeMFD = downsampledExclude.calcParticipationMFD_forParentSect(parentSect, 5.05, 8.55, origMFD.size());
			excludeMFD.setName("Excluding");
		}
		newMFD = downsampled.calcParticipationMFD_forParentSect(parentSect, 5.05, 8.55, origMFD.size());
		newMFD.setName("Downsampled");
		
		funcs = Lists.newArrayList();
		
		funcs.add(origMFD);
		if (excludeMFD != null)
			funcs.add(excludeMFD);
		funcs.add(newMFD);
		
		spec = new PlotSpec(funcs, chars, "Mojave Participation MFD Comparison", "Mag", "Incremental Rate (1/yr)");
		spec.setLegendVisible(true);
		gw = new GraphWindow(spec);
		gw.getGraphWidget().setBackgroundColor(Color.WHITE);
		gw.setYLog(true);
		gw.setY_AxisRange(1e-6, 1e-2);
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		if (plotDir != null) {
			try {
				gw.saveAsPNG(new File(plotDir, "mojave_mfd_comparison.png").getAbsolutePath());
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	static void plotGMPE_Curves(AttenuationRelationship gmpe, List<CybershakeSite> csSites,
			OrderedSiteDataProviderList provs, FaultSystemSolution orig, FaultSystemSolution downsampled,
			FaultSystemSolution downsampledExclude) {
		// set to 3s SA
		gmpe.setParamDefaults();
		gmpe.setIntensityMeasure(SA_Param.NAME);
		SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), 3d);
		HazardCurveCalculator calc = new HazardCurveCalculator();
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(gmpe.getIntensityMeasure());
		
		List<DiscretizedFunc> origCurves = Lists.newArrayList();
		List<DiscretizedFunc> newCurves = Lists.newArrayList();
		List<DiscretizedFunc> excludeCurves = Lists.newArrayList();
		
		SiteTranslator trans = new SiteTranslator();
		List<Site> sites = Lists.newArrayList();
		// create sites with site data
		for (CybershakeSite csSite : csSites) {
			Site site = new Site(csSite.createLocation());
			
			ArrayList<SiteDataValue<?>> datas = provs.getBestAvailableData(site.getLocation());
			
			for (Parameter<?> param : gmpe.getSiteParams()) {
				param = (Parameter<?>) param.clone();
				trans.setParameterValue(param, datas);
				site.addParameter(param);
			}
			sites.add(site);
		}
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(orig);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		erf.getTimeSpan().setDuration(1d);
		erf.updateForecast();
		
		System.out.println("Calculating original curves");
		for (Site site : sites) {
			DiscretizedFunc curve = HazardCurveSetCalculator.getLogFunction(xVals);
			calc.getHazardCurve(curve, site, gmpe, erf);
			origCurves.add(HazardCurveSetCalculator.unLogFunction(xVals, curve));
		}
		
		erf = new FaultSystemSolutionERF(downsampled);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		erf.getTimeSpan().setDuration(1d);
		erf.updateForecast();
		
		System.out.println("Calculating new curves");
		for (Site site : sites) {
			DiscretizedFunc curve = HazardCurveSetCalculator.getLogFunction(xVals);
			calc.getHazardCurve(curve, site, gmpe, erf);
			newCurves.add(HazardCurveSetCalculator.unLogFunction(xVals, curve));
		}
		
		if (downsampledExclude != null) {
			erf = new FaultSystemSolutionERF(downsampledExclude);
			erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
			erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
			erf.getTimeSpan().setDuration(1d);
			erf.updateForecast();
			
			System.out.println("Calculating exclude curves");
			for (Site site : sites) {
				DiscretizedFunc curve = HazardCurveSetCalculator.getLogFunction(xVals);
				calc.getHazardCurve(curve, site, gmpe, erf);
				excludeCurves.add(HazardCurveSetCalculator.unLogFunction(xVals, curve));
			}
		}
		
		for (int i=0; i<sites.size(); i++) {
			List<DiscretizedFunc> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			
			DiscretizedFunc origCurve = origCurves.get(i);
			origCurve.setName("Original");
			funcs.add(origCurve);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			if (downsampledExclude != null) {
				DiscretizedFunc excludeCurve = excludeCurves.get(i);
				excludeCurve.setName("Excluding");
				funcs.add(excludeCurve);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GREEN.darker()));
			}
			DiscretizedFunc newCurve = newCurves.get(i);
			newCurve.setName("Downsampled");
			funcs.add(newCurve);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
			
			String siteName = csSites.get(i).short_name;
			
			PlotSpec spec = new PlotSpec(funcs, chars, siteName+" GMPE Curves", "3sec SA (g)", "Exceed. Prob");
			spec.setLegendVisible(true);
			GraphWindow gw = new GraphWindow(spec);
			gw.getGraphWidget().setBackgroundColor(Color.WHITE);
			gw.setXLog(true);
			gw.setYLog(true);
			gw.setX_AxisRange(1e-2, 3);
			gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
			if (plotDir != null) {
				try {
					gw.saveAsPNG(new File(plotDir, "curve_"+siteName+"_comparison.png").getAbsolutePath());
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	private static File plotDir = null;

	public static void main(String[] args) throws IOException, DocumentException {
		FaultSystemSolution inputSol = FaultSystemIO.loadSol(
//				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
//						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/UCERF3_ERF/"
						+ "cached_dep10.0_depMean_rakeGEOLOGIC.zip"));
		
//		plotDir = new File("/home/kevin/CyberShake/ucerf3/subset_plots");
		plotDir = new File("/home/kevin/CyberShake/ucerf3/subset_plots_both_fm");
		
		if (plotDir != null && !plotDir.exists())
			plotDir.mkdir();
		
		Region region = new CaliforniaRegions.RELM_SOCAL();
		
		FaultSystemSolution downsampled = getDownsampledSol(inputSol, region);
		FaultSystemSolution downsampledExclude = getDownsampledSol(inputSol, region, true);
		plotMFD(region, inputSol, downsampled, downsampledExclude);
		
		OrderedSiteDataProviderList provs = HazardCurvePlotter.createProviders(5); // CVM-S4.26
		
		List<CybershakeSite> curveSites = Lists.newArrayList();
		DBAccess db = Cybershake_OpenSHA_DBApplication.getDB();
		CybershakeSiteInfo2DB sites2db = new CybershakeSiteInfo2DB(db);
		curveSites.add(sites2db.getSiteFromDB("STNI"));
		curveSites.add(sites2db.getSiteFromDB("SBSM"));
		curveSites.add(sites2db.getSiteFromDB("PAS"));
		curveSites.add(sites2db.getSiteFromDB("FIL"));
		db.destroy();
		
		NGAWest_2014_Averaged_AttenRel gmpe = new NGAWest_2014_Averaged_AttenRel(null, false);
		
		plotGMPE_Curves(gmpe, curveSites, provs, inputSol, downsampled, downsampledExclude);
	}

}
