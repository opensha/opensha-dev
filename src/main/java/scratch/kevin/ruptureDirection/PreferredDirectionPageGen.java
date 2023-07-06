package scratch.kevin.ruptureDirection;

import java.awt.Color;
import java.awt.Stroke;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletableFuture;
import java.util.function.Supplier;

import org.jfree.chart.annotations.XYPolygonAnnotation;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.NamedFaults;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.mod.ModAttenRelRef;
import org.opensha.sha.imr.mod.ModAttenuationRelationship;
import org.opensha.sha.imr.mod.impl.BaylessSomerville2013DirectivityModifier;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableSet;
import com.google.common.primitives.Ints;

import scratch.UCERF3.erf.FaultSystemSolutionERF;

public class PreferredDirectionPageGen {

	public static void main(String[] args) throws IOException {
		File markdownDir = new File("/home/kevin/markdown/preferred-rup-direction");
		
		String dirName;
		
		String erfName = "NSHM23 V2";
		dirName = "nshm23";
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
		// unload rupture MFDs (just use mean mags)
		sol.removeModuleInstances(RupMFDsModule.class);
		
		double fractInPreferred = 0.75;
		dirName += "-fractPref"+(float)fractInPreferred;
		double maxBilateralFract = 0.25;
		dirName += "-maxBilateral"+(float)maxBilateralFract;
		
		boolean plotModERFOnlyTest = false;
		
		int numHypos = 20;
		AttenRelRef gmmRef = AttenRelRef.ASK_2014;
		// must always be SA for Bayless
		double[] periods = {1d, 3d, 5d, 10d};
		
		List<PreferredDirection> directions = new ArrayList<>();
		
		Location mojaveSJtarget = new Location(34.2650935644413, -117.511879127267);
		
		List<? extends FaultSection> allSects = sol.getRupSet().getFaultSectionDataList();
		
		int mojaveSouthID = FaultSectionUtils.findParentSectionID(allSects, "mojave", "south");
		Preconditions.checkState(mojaveSouthID > 0);
		int mojaveNorthID = FaultSectionUtils.findParentSectionID(allSects, "mojave", "north");
		Preconditions.checkState(mojaveNorthID > 0);
		directions.add(build(sol, fractInPreferred, maxBilateralFract, "Mojave", mojaveSJtarget,
				mojaveSouthID, mojaveNorthID));
		dirName += "-mojave";
		
		NamedFaults namedFaults = sol.getRupSet().requireModule(NamedFaults.class);
		List<Integer> sjIDs = namedFaults.getParentIDsForFault("San Jacinto");
		directions.add(build(sol, fractInPreferred, maxBilateralFract, "San Jacinto", mojaveSJtarget,
				Ints.toArray(sjIDs)));
		dirName += "-san_jacinto";
		
		Location haywardTarget = new Location(37.039598378082836, -121.49086622123069);
		int haywardSouthID = FaultSectionUtils.findParentSectionID(allSects, "hayward", "south)");
		Preconditions.checkState(haywardSouthID > 0);
		int haywardSouthExtensionID = FaultSectionUtils.findParentSectionID(allSects, "hayward", "south", "extension");
		Preconditions.checkState(haywardSouthExtensionID > 0);
		int haywardNorthID = FaultSectionUtils.findParentSectionID(allSects, "hayward", "north)");
		Preconditions.checkState(haywardNorthID > 0);
		directions.add(build(sol, fractInPreferred, maxBilateralFract, "Hayward", haywardTarget,
				haywardSouthID, haywardSouthExtensionID, haywardNorthID));
		dirName += "-hayward";
		
		List<Location> debugSiteLocs = new ArrayList<>();
		List<String> debugSiteNames = new ArrayList<>();
		
		debugSiteLocs.add(new Location(34.10, -117.30));
		debugSiteNames.add("San Bernardino");
		
		debugSiteLocs.add(new Location(34.05, -118.25));
		debugSiteNames.add("Los Angeles");
		
		debugSiteLocs.add(new Location(34.79583396054309, -118.85289221749382));
		debugSiteNames.add("Gorman");
		
		debugSiteLocs.add(new Location(33.1441920080516, -116.13254218644279));
		debugSiteNames.add("Ocotillo Wells");
		
		debugSiteLocs.add(new Location(37.35, -121.90));
		debugSiteNames.add("San Jose");
		
		debugSiteLocs.add(new Location(37.935708337060895, -122.3475369629483));
		debugSiteNames.add("Richmond");
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.getTimeSpan().setDuration(1d);
		
		PreferredDirectionERF modERF = new PreferredDirectionERF(erf, numHypos, directions);
		modERF.updateForecast();
		
		File outputDir = new File(markdownDir, dirName);
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# Preferred Rupture Direction");
		lines.add("");
		lines.add("This page gives hazard calculations using "+erfName+", and assuming preffered rupture directions on "
				+ "the select sections.");
		lines.add("");
		DecimalFormat pDF = new DecimalFormat("0.##%");
		lines.add("For each rupture involving those faults, ruptures that nucleate on the given fault will rupture in the "
				+ "preferred direction "+pDF.format(fractInPreferred)+" of the time. A rupture is deemed to be in the "
				+ "preferred direction so long as no more than "+pDF.format(maxBilateralFract)+" of the rupture propagates "
				+ "in the opposite direction.");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Preferred-Direction Faults");
		lines.add(topLink); lines.add("");
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		Random rand = new Random(rupSet.getNumRuptures()*directions.size());
		for (PreferredDirection direction : directions) {
			lines.add("### "+direction.getName());
			lines.add(topLink); lines.add("");
			
			String prefix = direction.getName().replaceAll("\\W+", "_");
			writeColorWheelDirectionPlot(rupSet, List.of(direction), null, resourcesDir, prefix);
			
			lines.add("#### "+direction.getName()+" Direction Map");
			lines.add(topLink); lines.add("");
			
			lines.add("![Direction Map]("+resourcesDir.getName()+"/"+prefix+".png)");
			lines.add("");
			
			// first try to find one that uses all sections
			ImmutableSet<Integer> allRups = direction.getAffectedRuptures();
			HashSet<Integer> targetSectIDs = new HashSet<>();
			for (FaultSection sect : direction.getTargetSects())
				targetSectIDs.add(sect.getSectionId());
			int largestAllOnIndex = -1;
			int largestAllOnSize = 0;
			for (int rupIndex : allRups) {
				boolean allOn = true;
				List<Integer> rupSects = rupSet.getSectionsIndicesForRup(rupIndex);
				if (rupSects.size() <= largestAllOnSize)
					continue;
				for (int sectIndex : rupSects) {
					if (!targetSectIDs.contains(sectIndex)) {
						allOn = false;
						break;
					}
				}
				if (allOn) {
					largestAllOnIndex = rupIndex;
					largestAllOnSize = rupSects.size();
				}
			}
			
			lines.add("#### "+direction.getName()+" Example Rupture Preferred Hypocenters");
			lines.add(topLink); lines.add("");
			
			if (largestAllOnIndex > 0) {
				String rupPrefix = prefix+"_rup"+largestAllOnIndex;
				writeExampleRupPlot(rupSet, direction, largestAllOnIndex, resourcesDir, rupPrefix, false);
				
				lines.add("__Largest Rupture Exclusively on "+direction.getName()+"__");
				lines.add("");
				lines.add("![Rupture Example Map]("+resourcesDir.getName()+"/"+rupPrefix+".png)");
				lines.add("");
			}
			
			int numRand = 4;
			List<Integer> rupsList = new ArrayList<>(allRups);
			Collections.sort(rupsList);
			Collections.shuffle(rupsList, rand);
			rupsList = rupsList.subList(0, numRand);
			TableBuilder table = MarkdownUtils.tableBuilder();
			table.initNewLine();
			for (int rupIndex : rupsList) {
				String rupPrefix = prefix+"_rup"+rupIndex;
				writeExampleRupPlot(rupSet, direction, rupIndex, resourcesDir, rupPrefix, false);
				
				table.addColumn("![Rupture Example Map]("+resourcesDir.getName()+"/"+rupPrefix+".png)");
			}
			table.finalizeLine();
			lines.add("__Random Example Ruptures__");
			lines.add("");
			lines.addAll(table.wrap(2, 0).build());
			lines.add("");
		}
		
		ScalarIMR rawGMM = gmmRef.get();
		ScalarIMR directivityGMM = buildDirectivityGMM(gmmRef);
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		DiscretizedFunc[] periodXVals = new DiscretizedFunc[periods.length];
		DiscretizedFunc[] periodLogXVals = new DiscretizedFunc[periods.length];
		String[] periodPrefixes = new String[periods.length];
		String[] periodNames = new String[periods.length];
		String[] periodUnits = new String[periods.length];
		IMT_Info imtInfo = new IMT_Info();
		DecimalFormat oDF = new DecimalFormat("0.##");
		for (int p=0; p<periods.length; p++) {
			if (periods[p] == 0d) {
				periodXVals[p] = imtInfo.getDefaultHazardCurve(PGA_Param.NAME);
				periodPrefixes[p] = "pga";
				periodNames[p] = "PGA";
				periodUnits[p] = "(g)";
			} else {
				Preconditions.checkState(periods[p] > 0);
				periodXVals[p] = imtInfo.getDefaultHazardCurve(SA_Param.NAME);
				periodPrefixes[p] = "sa_"+oDF.format(periods[p])+"s";
				periodNames[p] = oDF.format(periods[p])+"s SA";
				periodUnits[p] = "(g)";
			}
			periodLogXVals[p] = new ArbitrarilyDiscretizedFunc();
			for (Point2D pt : periodXVals[p])
				periodLogXVals[p].set(Math.log(pt.getX()), 0d);
		}
		
		modERF.cacheSourcesParallel();
		
		ArrayDeque<ScalarIMR> rawGMMDeque = new ArrayDeque<>();
		ArrayDeque<ScalarIMR> directivityGMMDeque = new ArrayDeque<>();
		
		GMM_Manager gmmManager = new GMM_Manager(gmmRef);
		
		for (int i=0; i<debugSiteLocs.size(); i++) {
			Site site = new Site(debugSiteLocs.get(i));
			site.setName(debugSiteNames.get(i));
			site.addParameterList(rawGMM.getSiteParams());
			
			lines.add("## "+site.getName()+" Hazard");
			lines.add(topLink); lines.add("");
			
			String sitePrefix = site.getName().replaceAll("\\W+", "_");
			
			String mapPrefix = sitePrefix+"_map";
			writeColorWheelDirectionPlot(rupSet, getWithinDist(site.getLocation(), directions, 200d),
					List.of(site.getLocation()), resourcesDir, mapPrefix);
			
			lines.add("### "+site.getName()+" Map");
			lines.add(topLink); lines.add("");
			lines.add("");
			lines.add("![Site map]("+resourcesDir.getName()+"/"+mapPrefix+".png)");
			lines.add("");
			
			lines.add("### "+site.getName()+" Hazard Curves");
			lines.add(topLink); lines.add("");
			
			System.out.println("Calculating hazard for "+site.getName());
			TableBuilder curveTable = MarkdownUtils.tableBuilder();
			curveTable.initNewLine();
			
			List<CompletableFuture<DiscretizedFunc>> origFutures = new ArrayList<>(periods.length);
			List<CompletableFuture<DiscretizedFunc>> modERFFutures = new ArrayList<>(periods.length);
			List<CompletableFuture<DiscretizedFunc>> directivityFutures = new ArrayList<>(periods.length);
			for (int p=0; p<periods.length; p++) {
				directivityFutures.add(CompletableFuture.supplyAsync(new CurveCalcCallable(site, periods[p],
						periodLogXVals[p], periodXVals[p], modERF, gmmManager, true)));
				modERFFutures.add(CompletableFuture.supplyAsync(new CurveCalcCallable(site, periods[p],
						periodLogXVals[p], periodXVals[p], modERF, gmmManager, false)));
				origFutures.add(CompletableFuture.supplyAsync(new CurveCalcCallable(site, periods[p],
						periodLogXVals[p], periodXVals[p], erf, gmmManager, false)));
			}
			
			TableBuilder twoIn50Table = MarkdownUtils.tableBuilder();
			twoIn50Table.addLine("", "Original", "Directivity & Direction Modified", "Gain");
			double twoInFiftyProb = ReturnPeriods.TWO_IN_50.oneYearProb;
			String twoInFiftyLabel = ReturnPeriods.TWO_IN_50.label;
			
			for (int p=0; p<periods.length; p++) {
				System.out.println("Period: "+periodNames[p]);
				
				if (periods[p] == 0d) {
					rawGMM.setIntensityMeasure(PGA_Param.NAME);
					directivityGMM.setIntensityMeasure(PGA_Param.NAME);
				} else {
					Preconditions.checkState(periods[p] > 0d);
					rawGMM.setIntensityMeasure(SA_Param.NAME);
					SA_Param.setPeriodInSA_Param(rawGMM.getIntensityMeasure(), periods[p]);
					directivityGMM.setIntensityMeasure(SA_Param.NAME);
					SA_Param.setPeriodInSA_Param(directivityGMM.getIntensityMeasure(), periods[p]);
				}
				
				DiscretizedFunc origCurve = origFutures.get(p).join();
				
				DiscretizedFunc modERFOnlyCurve = modERFFutures.get(p).join();
				
				DiscretizedFunc directivityCurve = directivityFutures.get(p).join();
				
				System.out.println("X\tOrig\tModERF\tDirectivity");
				for (int j=0; j<origCurve.size(); j++)
					System.out.println((float)origCurve.getX(j)+"\t"+(float)origCurve.getY(j)
						+"\t"+(float)modERFOnlyCurve.getY(j)+"\t"+(float)directivityCurve.getY(j));
				
				List<DiscretizedFunc> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				origCurve.setName("No Directivity");
				funcs.add(origCurve);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
				
				if (plotModERFOnlyTest) {
					modERFOnlyCurve.setName("Mod-ERF-Only Test");
					funcs.add(modERFOnlyCurve);
					chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.GRAY));
				}
				
				directivityCurve.setName("Directivity & Preffered Direction");
				funcs.add(directivityCurve);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.RED));
				
				double origGM = origCurve.getFirstInterpolatedX_inLogXLogYDomain(twoInFiftyProb);
				double directivityGM = directivityCurve.getFirstInterpolatedX_inLogXLogYDomain(twoInFiftyProb);
				double gain = directivityGM/origGM;
				twoIn50Table.addLine(periodNames[p], (float)origGM+" "+periodUnits[p], (float)directivityGM+" "+periodUnits[p], (float)gain);
				
				PlotSpec spec = new PlotSpec(funcs, chars, site.getName(),
						periodNames[p]+" "+periodUnits[p], "Annual Probability of Exceedance");
				spec.setLegendInset(true);
				
				Range xRange = new Range(1e-3, 1e1);
				Range yRange = new Range(1e-6, 1e0);
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				gp.drawGraphPanel(spec, true, true, xRange, yRange);
				
				String prefix = sitePrefix+"_curves_"+periodPrefixes[p];
				
				PlotUtils.writePlots(resourcesDir, prefix, gp, 800, 750, true, true, false);
				
				curveTable.addColumn("![Hazard Curves]("+resourcesDir.getName()+"/"+prefix+".png)");
			}
			curveTable.finalizeLine();
			
			lines.addAll(curveTable.wrap(2, 0).build());
			lines.add("");
			lines.add("__"+twoInFiftyLabel+" Hazard__");
			lines.add("");
			lines.addAll(twoIn50Table.build());
			lines.add("");
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private static class CurveCalcCallable implements Supplier<DiscretizedFunc> {
		private Site site;
		private double period;
		private DiscretizedFunc logXVals;
		private DiscretizedFunc xVals;
		private AbstractERF erf;
		private GMM_Manager gmmManager;
		private boolean directivity;

		public CurveCalcCallable(Site site, double period, DiscretizedFunc logXVals, DiscretizedFunc xVals,
				AbstractERF erf, GMM_Manager gmmManager, boolean directivity) {
			this.site = site;
			this.period = period;
			this.logXVals = logXVals;
			this.xVals = xVals;
			this.erf = erf;
			this.gmmManager = gmmManager;
			this.directivity = directivity;
		}

		@Override
		public DiscretizedFunc get() {
			ScalarIMR gmm = gmmManager.getGMM(directivity);
			HazardCurveCalculator calc = gmmManager.getCalc();
			
			if (period == 0d) {
				gmm.setIntensityMeasure(PGA_Param.NAME);
			} else {
				gmm.setIntensityMeasure(SA_Param.NAME);
				SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), period);
			}
			
			logXVals = logXVals.deepClone();
			
			calc.getHazardCurve(logXVals, site, gmm, erf);
			
			gmmManager.returnGMM(gmm, directivity);
			gmmManager.returnCalc(calc);
			
			return cloneLinearXVals(logXVals, xVals);
		}
		
	}
	
	private static class GMM_Manager {
		private ArrayDeque<HazardCurveCalculator> calcDeque = new ArrayDeque<>();
		
		private ArrayDeque<ScalarIMR> rawDeque = new ArrayDeque<>();
		private ArrayDeque<ScalarIMR> directivityDeque = new ArrayDeque<>();
		private AttenRelRef gmmRef;
		
		public GMM_Manager(AttenRelRef gmmRef) {
			this.gmmRef = gmmRef;
		}
		
		public synchronized ScalarIMR getGMM(boolean directivity) {
			if (directivity) {
				if (directivityDeque.isEmpty())
					return buildDirectivityGMM(gmmRef);
				return directivityDeque.remove();
			} else {
				if (rawDeque.isEmpty())
					return gmmRef.get();
				return rawDeque.remove();
			}
		}
		
		public synchronized void returnGMM(ScalarIMR gmm, boolean directivity) {
			if (directivity) {
				directivityDeque.push(gmm);;
			} else {
				rawDeque.push(gmm);;
			}
		}
		
		public synchronized HazardCurveCalculator getCalc() {
			if (calcDeque.isEmpty())
				return new HazardCurveCalculator();
			return calcDeque.remove();
		}
		
		public synchronized void returnCalc(HazardCurveCalculator calc) {
			calcDeque.push(calc);
		}
	}
	
	private static DiscretizedFunc cloneLinearXVals(DiscretizedFunc logCurve, DiscretizedFunc xVals) {
		DiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<xVals.size(); i++)
			ret.set(xVals.getX(i), logCurve.getY(i));
		return ret;
	}
	
	private static ScalarIMR buildDirectivityGMM(AttenRelRef ref) {
		ModAttenuationRelationship modGMM = new ModAttenuationRelationship(
				ref, ModAttenRelRef.BAYLESS_SOMERVILLE_2013_DIRECTIVITY);
		modGMM.setParamDefaults();
		return modGMM;
	}
	
	private static PreferredDirection build(FaultSystemSolution sol, double fractInPreferred, double maxBilateralFract,
			String name, Location towardsLoc, int... parentIDs) {
		boolean aseisReducesArea = true;
		
		List<FaultSection> targetSects = new ArrayList<>();
		List<Boolean> targetAlongStrikes = new ArrayList<>();
		
		for (FaultSection sect : sol.getRupSet().getFaultSectionDataList()) {
			boolean match = false;
			for (int parentID : parentIDs) {
				if (sect.getParentSectionId() == parentID) {
					match = true;
					break;
				}
			}
			if (match) {
				FaultTrace trace = sect.getFaultTrace();
				double firstDist = LocationUtils.horzDistance(trace.first(), towardsLoc);
				double lastDist = LocationUtils.horzDistance(trace.last(), towardsLoc);
				boolean alongStrike = firstDist > lastDist;
				targetSects.add(sect);
				targetAlongStrikes.add(alongStrike);
			}
		}
		
		return new PreferredDirection(name, sol, targetSects, targetAlongStrikes,
				fractInPreferred, maxBilateralFract, aseisReducesArea);
	}
	
	private static List<PreferredDirection> getWithinDist(Location siteLoc, List<PreferredDirection> directions, double maxDist) {
		List<PreferredDirection> ret = new ArrayList<>();
		
		for (PreferredDirection direction : directions) {
			boolean match = false;
			for (FaultSection sect : direction.getTargetSects()) {
				for (Location loc : sect.getFaultTrace()) {
					double sectDist = LocationUtils.horzDistanceFast(loc, siteLoc);
					if (sectDist < maxDist) {
						match = true;
						break;
					}
				}
			}
			if (match)
				ret.add(direction);
		}
		
		return ret;
	}
	
	private static void writeColorWheelDirectionPlot(FaultSystemRupSet rupSet, List<PreferredDirection> directions,
			List<Location> sites, File resourcesDir, String prefix) throws IOException {
		MinMaxAveTracker latTrack = new MinMaxAveTracker();
		MinMaxAveTracker lonTrack = new MinMaxAveTracker();
		
		CPT dirCPT = new CPT();
		Color prevColor = null;
		double bright = 0.9d;
		for (int i=0; i<=360; i++) {
			double hue = (double)i/360d;
			Color color = Color.getHSBColor((float)hue, 1f, (float)bright);
			if (prevColor != null)
				dirCPT.add(new CPTVal(i-1, prevColor, i, color));
			prevColor = color;
		}
		
		List<Double> sectAzimuths = new ArrayList<>(rupSet.getNumSections());
		for (int s=0; s<rupSet.getNumSections(); s++)
			sectAzimuths.add(Double.NaN);
		
		for (PreferredDirection direction : directions) {
			List<FaultSection> sects = direction.getTargetSects();
			List<Boolean> alongStrikes = direction.getTargetAlongStrikes();
			
			for (int i=0; i<sects.size(); i++) {
				FaultSection sect = sects.get(i);
				boolean alongStrike = alongStrikes.get(i);
				FaultTrace trace = sect.getFaultTrace();
				latTrack.addValue(trace.first().lat);
				lonTrack.addValue(trace.first().lon);
				latTrack.addValue(trace.last().lat);
				lonTrack.addValue(trace.last().lon);
				double dir = trace.getStrikeDirection();
				if (!alongStrike) {
					// reverse it
					dir += 180;
					if (dir > 360d)
						dir -= 360;
				}
				sectAzimuths.set(sect.getSectionId(), dir);
			}
		}
		
		if (sites != null) {
			for (Location loc : sites) {
				latTrack.addValue(loc.lat);
				lonTrack.addValue(loc.lon);
			}
		}
		
		Region region = bufferedRegion(latTrack, lonTrack);
		
		RupSetMapMaker mapMaker = new RupSetMapMaker(rupSet, region);
		
		if (sites != null && !sites.isEmpty()) {
			mapMaker.plotScatters(sites, Color.BLACK);
			mapMaker.setScatterSymbol(PlotSymbol.FILLED_INV_TRIANGLE, 5f);
		}
		mapMaker.setSkipNaNs(true);
		mapMaker.plotSectScalars(sectAzimuths, dirCPT, "Preferred Propagation Direction (deg)");
		
		PlotSpec spec = mapMaker.buildPlot(" ");
		// build color wheel annotation
		
		// figure out aspect ratio of lat to lon
		Location centerLoc = new Location(0.5*(region.getMinLat()+region.getMaxLat()),
				0.5*(region.getMinLon()+region.getMaxLon()));
		double latLen = LocationUtils.horzDistance(centerLoc, new Location(centerLoc.getLatitude()+1, centerLoc.getLongitude()));
		double lonLen = LocationUtils.horzDistance(centerLoc, new Location(centerLoc.getLatitude(), centerLoc.getLongitude()+1));
		double aspect = latLen/lonLen;
		
		// center location of the wheel
		double centerOffset = 0.12*(region.getMaxLon()-region.getMinLon());
		double wheelCenterX = region.getMinLon() + centerOffset;
		double wheelCenterY = region.getMinLat() + centerOffset/aspect;
//		double widthX = 0.20*(sourceRegion.getMaxLat()-sourceRegion.getMinLat());
//		double aspect = PlotUtils.calcAspectRatio(
//				new Range(sourceRegion.getMinLon(), sourceRegion.getMaxLon()),
//				new Range(sourceRegion.getMinLat(), sourceRegion.getMaxLat()), true);
		
		double innerMult = 0.3;
		double outerMult = 0.4;
//		double sumY = Math.max(1d, hist.calcSumOfY_Vals());
		double halfDelta = 15d;
		for (double centerAz=0d; centerAz<=350d; centerAz+=30d) {
			double startAz = Math.toRadians(centerAz-halfDelta);
			double endAz = Math.toRadians(centerAz+halfDelta);
			
			List<Point2D> points = new ArrayList<>();
			
			double startX = Math.sin(startAz);
			double startY = Math.cos(startAz);
			double endX = Math.sin(endAz);
			double endY = Math.cos(endAz);
			
			points.add(new Point2D.Double(innerMult*startX, innerMult*startY));
			points.add(new Point2D.Double(outerMult*startX, outerMult*startY));
			points.add(new Point2D.Double(outerMult*endX, outerMult*endY));
			points.add(new Point2D.Double(innerMult*endX, innerMult*endY));
			points.add(new Point2D.Double(innerMult*startX, innerMult*startY));
			
			double[] polygon = new double[points.size()*2];
			int cnt = 0;
			for (Point2D pt : points) {
				polygon[cnt++] = pt.getX()+wheelCenterX;
				polygon[cnt++] = pt.getY()/aspect+wheelCenterY;
			}
			Color color = dirCPT.getColor((float)centerAz);
			
			Stroke stroke = PlotLineType.SOLID.buildStroke(1f);
			spec.addPlotAnnotation(new XYPolygonAnnotation(polygon, stroke, Color.DARK_GRAY, color));
		}
		mapMaker.plot(resourcesDir, prefix, spec, 900);
	}
	
	private static void writeExampleRupPlot(FaultSystemRupSet rupSet, PreferredDirection direction,
			int rupIndex, File resourcesDir, String prefix, boolean debug) throws IOException {
		MinMaxAveTracker latTrack = new MinMaxAveTracker();
		MinMaxAveTracker lonTrack = new MinMaxAveTracker();
		
		List<FaultSection> rupSects = rupSet.getFaultSectionDataForRupture(rupIndex);
		List<Location> hypos = direction.drawPreferredDirectionHypos(rupIndex, -1, debug);
		
		for (FaultSection sect : rupSects) {
			FaultTrace trace = sect.getFaultTrace();
			
			latTrack.addValue(trace.first().lat);
			lonTrack.addValue(trace.first().lon);
			latTrack.addValue(trace.last().lat);
			lonTrack.addValue(trace.last().lon);
		}
		
		Region region = bufferedRegion(latTrack, lonTrack);
		RupSetMapMaker mapMaker = new RupSetMapMaker(rupSet, region);
		
		mapMaker.highLightSections(rupSects, new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		mapMaker.plotScatters(hypos, Color.RED.darker());
		mapMaker.setScatterSymbol(PlotSymbol.FILLED_CIRCLE, 3f);
		
		mapMaker.plot(resourcesDir, prefix, direction.getName()+" Preffered Hypocenters for Rupture "+rupIndex);
	}
	
	private static Region bufferedRegion(MinMaxAveTracker latTrack, MinMaxAveTracker lonTrack) {
		double minLat = latTrack.getMin();
		double maxLat = latTrack.getMax();
		double minLon = lonTrack.getMin();
		double maxLon = lonTrack.getMax();
		Location upperLeft = new Location(maxLat, minLon);
		Location bottomRight = new Location(minLat, maxLon);
		// buffer them
		upperLeft = LocationUtils.location(upperLeft, 1.75*Math.PI, 100d);
		bottomRight = LocationUtils.location(bottomRight, 0.75*Math.PI, 100d);
		return new Region(bottomRight, upperLeft);
	}

}
