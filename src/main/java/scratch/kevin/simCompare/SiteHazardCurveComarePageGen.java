package scratch.kevin.simCompare;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Base64;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.ParameterList;
import org.opensha.commons.param.impl.DoubleParameter;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileUtils;
import org.opensha.sha.calc.disaggregation.DisaggregationCalculator;
import org.opensha.sha.calc.disaggregation.DisaggregationCalculatorAPI;
import org.opensha.sha.calc.hazardMap.HazardDataSetLoader;
import org.opensha.sha.calc.params.filters.SourceFilter;
import org.opensha.sha.calc.params.filters.SourceFiltersParam;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import com.google.common.base.Preconditions;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;

import scratch.kevin.simCompare.SimulationDisaggAttenuationRelationshipWrapper.Source;

import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.ReturnPeriodUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;

public abstract class SiteHazardCurveComarePageGen<E> {
	private SimulationRotDProvider<E> simProv;
	SimulationHazardCurveCalc<E> simCalc;
	private String simName;
	
	private static double[] gmpe_truncs = { 3d, 2d, 1d };
//	private static double[] gmpe_fixed_sigmas = { 0.5, 0.3, 0d };
//	private static double[] gmpe_truncs = {  };
	private static double[] gmpe_fixed_sigmas = {  };
	
	private static double[] disagg_fixed_accel_vals = { 0.1, 0.5, 1.0 };
	
	private static ExecutorService exec;
	
	private Map<SimulationHazardCurveCalc<?>, PlotCurveCharacterstics> customPlotChars = new HashMap<>();
	
	private List<SimulationRotDProvider<?>> compSimProvs;
	private List<SimulationHazardCurveCalc<?>> compCurveCals;
	
	private static Map<AttenRelRef, LinkedList<ScalarIMR>> gmpesInstancesCache;
	
	private LinkedList<DisaggCalc> disaggCalcsCache;
	
	private boolean replotCurves = true;
	private boolean replotDisaggs = true;
	
	private double sourceRupContributionSortProb;
	private int sourceRupContributionNum;
	private Table<String, E, Double> sourceRupContributionFracts;
	private boolean sourceRupContributionsMutuallyExclusive;

	public SiteHazardCurveComarePageGen(SimulationRotDProvider<E> simProv, String simName) {
		this(simProv, simName, new ArrayList<>());
	}

	public SiteHazardCurveComarePageGen(SimulationRotDProvider<E> simProv, String simName,
			SimulationRotDProvider<?>... compSimProvs) {
		this(simProv, simName, toList(compSimProvs));
	}

	public SiteHazardCurveComarePageGen(SimulationRotDProvider<E> simProv, String simName,
			List<SimulationRotDProvider<?>> compSimProvs) {
		super();
		this.simProv = simProv;
		simCalc = new SimulationHazardCurveCalc<>(simProv);
		this.simName = simName;
		if (compSimProvs == null)
			compSimProvs = new ArrayList<>();
		this.compSimProvs = compSimProvs;
		compCurveCals = new ArrayList<>();
		for (SimulationRotDProvider<?> compSimProv : compSimProvs)
			compCurveCals.add(new SimulationHazardCurveCalc<>(compSimProv));
		
		gmpesInstancesCache = new HashMap<>();
	}
	
	public void setReplotCurves(boolean replotCurves) {
		this.replotCurves = replotCurves;
	}

	public void setReplotDisaggs(boolean replotDisaggs) {
		this.replotDisaggs = replotDisaggs;
	}
	
	public void setSourceRupContributionFractions(Table<String, E, Double> sourceRupContribFracts,
			double sortProb, int numSourcesToPlot) {
		this.sourceRupContributionFracts = sourceRupContribFracts;
		this.sourceRupContributionSortProb = sortProb;
		this.sourceRupContributionNum = numSourcesToPlot;
		if (sourceRupContribFracts != null) {
			sourceRupContributionsMutuallyExclusive = true;
			for (E rup : sourceRupContribFracts.columnKeySet()) {
				Map<String, Double> col = sourceRupContribFracts.column(rup);
				if (col.size() > 1) {
					int numNonzero = 0;
					for (Double val : col.values())
						if (val > 0)
							numNonzero++;
					if (numNonzero > 1) {
						sourceRupContributionsMutuallyExclusive = false;
					}
				}
			}
			System.out.println("source contributions mutually exclusive ? "+sourceRupContributionsMutuallyExclusive);
		}
	}
	
	public void setCustomPlotColors(SimulationRotDProvider<?> simProv, PlotCurveCharacterstics plotChar) {
		for (int i=0; i<compSimProvs.size(); i++) {
			if (simProv == compSimProvs.get(i)) {
				customPlotChars.put(compCurveCals.get(i), plotChar);
				return;
			}
		}
		throw new IllegalStateException("No comparison calc found for "+simProv.getName());
	}

	private static List<SimulationRotDProvider<?>> toList(SimulationRotDProvider<?>... compSimProvs) {
		if (compSimProvs == null || compSimProvs.length == 0)
			return new ArrayList<>();
		List<SimulationRotDProvider<?>> list = new ArrayList<>();
		for (SimulationRotDProvider<?> simProv : compSimProvs)
			list.add(simProv);
		return list;
	}
	
	public void generateSitePage(Site site, List<? extends RuptureComparison<E>> comps, File outputDir, List<String> headerLines,
			IMT[] imts, AttenRelRef gmpeRef) throws IOException {
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		List<String> lines = new LinkedList<>();
		if (headerLines != null && !headerLines.isEmpty()) {
			lines.addAll(headerLines);
			if (!lines.get(lines.size()-1).isEmpty())
				lines.add("");
		}
		
		lines.add("## Site Information");
		lines.add("");
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.addLine("**Name**", site.getName());
		table.addLine("**Latitude**", (float)site.getLocation().getLatitude()+"");
		table.addLine("**Longitude**", (float)site.getLocation().getLongitude()+"");
		table.addLine("**GMPE Parameters**", "");
		for (Parameter<?> param : site) {
			String name = "**"+param.getName()+"**";
			String units = param.getUnits();
			if (units != null && !units.isEmpty())
				name += " (*"+units+"*)";
			if (param instanceof DoubleParameter)
				table.addLine(name, ((Double)param.getValue()).floatValue()+"");
			else
				table.addLine(name, param.getValue().toString());
		}
		lines.addAll(table.build());
		lines.add("");
		try {
			File mapFile = new File(resourcesDir, "site_location_map.png");
			if (!mapFile.exists())
				FileUtils.downloadURL(new URL(getMiniMap(site.getLocation())), mapFile);
			lines.add("### Site Map");
			lines.add("");
			lines.add("![Site Map]("+resourcesDir.getName()+"/"+mapFile.getName()+")");
			lines.add("");
		} catch (Exception e1) {
			System.out.println("Error fetching map! Skipping");
			e1.printStackTrace();
		}
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		double curveDuration = 1d;
		
		SimulationHazardPlotter<E> curvePlotter = new SimulationHazardPlotter<>(simCalc, compCurveCals, comps, site, curveDuration, gmpeRef);
		curvePlotter.setGMPE_FixedSigmas(gmpe_fixed_sigmas);
		curvePlotter.setGMPE_TruncationLevels(gmpe_truncs);
		curvePlotter.setSourceContributionFractions(sourceRupContributionFracts, sourceRupContributionsMutuallyExclusive);
		for (SimulationHazardCurveCalc<?> key : customPlotChars.keySet())
			curvePlotter.setCustomPlotColors(key, customPlotChars.get(key));
		
		lines.add("## Curve Values Table");
		lines.add(topLink); lines.add("");
		
		table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		table.addColumn("Hazard Level");
		for (IMT imt : imts) {
			table.addColumn(simName+" "+imt.getDisplayName());
			table.addColumn(gmpeRef.getShortName()+" "+imt.getDisplayName());
		}
		table.finalizeLine();
		int[] rps = MultiRupGMPE_ComparePageGen.hazard_curve_rps;
		List<Double> probLevels = new ArrayList<>();
		List<String> probLabels = new ArrayList<>();
		for (int rp : rps) {
			probLabels.add(rp+"yr");
			probLevels.add((curveDuration/(double)rp));
		}
		probLabels.add("RTGM");
		probLevels.add(0d);
		
		for (int i=0; i<probLabels.size(); i++) {
			double level = probLevels.get(i);
			table.initNewLine();
			table.addColumn(probLabels.get(i));
			for (IMT imt : imts) {
				DiscretizedFunc[] curves = {
						curvePlotter.getCalcSimCurve(simCalc, imt),
						curvePlotter.getCalcGMPECurve(imt)
				};
				for (DiscretizedFunc curve : curves) {
					if (level == 0d) {
						// RTGM
						double rtgm = SimulationHazardPlotter.calcRTGM(curve, curveDuration);
						if (rtgm > 0d)
							table.addColumn((float)rtgm+"");
						else
							table.addColumn("N/A");
					} else  {
						if (level < curve.getMinY() || level > curve.getMaxY())
							table.addColumn("N/A");
						else
							table.addColumn((float)curve.getFirstInterpolatedX_inLogXLogYDomain(level)+"");
					}
				}
			}
			table.finalizeLine();
		}
		lines.addAll(table.build());
		lines.add("");
		
		lines.add("## Hazard Spectra");
		lines.add(topLink); lines.add("");
		lines.addAll(curvePlotter.getCurveLegend(true, true, true, false, 0));
		lines.add("");
		
		List<Double> saPeriodList = new ArrayList<>();
		for (IMT imt : imts)
			if (imt.getParamName().equals(SA_Param.NAME))
				saPeriodList.add(imt.getPeriod());
		double[] saPeriods = Doubles.toArray(saPeriodList);
		
		table = MarkdownUtils.tableBuilder();
		String spectraPrefix = site.getName().replaceAll(" ", "_")+"_spectra_"+gmpeRef.getShortName();
		for (int rp : rps) {
			String rpPrefix = spectraPrefix+"_"+rp+"yr";
			double probLevel = 1d/(double)rp;
			if (replotCurves || !new File(resourcesDir, rpPrefix+".png").exists())
				curvePlotter.plotHazardSpectra(resourcesDir, rpPrefix, probLevel, rp+"yr RotD50 (g)", saPeriods);
			
			table.addLine("**"+rp+"yr**", "![Hazard Spectra]("+resourcesDir.getName()+"/"+rpPrefix+".png)");
		}
		lines.addAll(table.build());
		
		if (sourceRupContributionFracts != null) {
			String sourceContribPrefix = spectraPrefix+"_source_contrib";
			
			File[] simFiles = new File[rps.length];
			File[] gmpeFiles = new File[rps.length];
			
			for (int i=0; i<rps.length; i++) {
				int rp = rps[i];
				double probLevel = 1d/(double)rp;
				String yAxisLabel = rp+"yr RotD50 (g)";
				String rpPrefix = sourceContribPrefix+"_"+rp+"yr";
				if (replotCurves || !new File(resourcesDir, rpPrefix+"_gmpe.png").exists()) {
					System.out.println("Plotting simulation source curves");
					simFiles[i] = curvePlotter.plotSourceContributionHazardSpectra(resourcesDir, rpPrefix+"_sim",
							saPeriods, probLevel, yAxisLabel, sourceRupContributionNum, false);
					System.out.println("Calculating/plotting GMPE source curves");
					gmpeFiles[i] = curvePlotter.plotSourceContributionHazardSpectra(resourcesDir, rpPrefix+"_gmpe",
							saPeriods, probLevel, yAxisLabel, sourceRupContributionNum, true);
				} else {
					simFiles[i] = new File(resourcesDir, rpPrefix+"_sim.png");
					gmpeFiles[i] = new File(resourcesDir, rpPrefix+"_gmpe.png");
				}
			}
			
			
			lines.add("### Source Contribution Spectra");
			lines.add(topLink); lines.add("");
			lines.add("");
			lines.add("These plots show the contribution of each fault source to the hazard spectra. The same set of "
					+ "sources are plotted for both simulation values (left) and GMPE values (right). Sources are sorted in the legend "
					+ "(and colored by) their average contrubution in the simulation results at the given return period, "
					+ "and only the top "+sourceRupContributionNum+" sources are plotted.");
			lines.add("");
			if (!sourceRupContributionsMutuallyExclusive) {
				lines.add("*NOTE: Source curves are not mututally exclusive. For the case of multi fault ruptures, "
						+ "a single rupture can be included in the curve for multiple sources*");
				lines.add("");
			}
			
			table = MarkdownUtils.tableBuilder();
			
			table.addLine("**Return Period**", "**Simulation Source Contributions**", "**GMPE Source Contributions**");
			for (int i=0; i<rps.length; i++)
				table.addLine(rps[i]+"yr", "![Hazard Spectra]("+resourcesDir.getName()+"/"+simFiles[i].getName()+")",
						"![Hazard Spectra]("+resourcesDir.getName()+"/"+gmpeFiles[i].getName()+")");
			lines.addAll(table.build());
			lines.add("");
		}
		
		lines.add("## Hazard Curves");
		lines.add(topLink); lines.add("");
		
		lines.addAll(curvePlotter.getCurveLegend(false, true, true, false, 0));
		lines.add("");
		
		GMPESimulationBasedProvider<E> gmpeSimProv = new GMPESimulationBasedProvider<>(
				simProv, comps, gmpeRef.getShortName()+" Simulatios", imts);
		
		boolean hasMultiple = true;
		for (E rup : simProv.getRupturesForSite(site)) {
			hasMultiple = hasMultiple && simProv.getNumSimulations(site, rup) > 1;
			if (!hasMultiple)
				break;
		}
		
		System.out.println("Have multiple sims per rup? "+hasMultiple);
		
		for (IMT imt : imts) {
			System.out.println("Calculating primary hazard curve");
			
			String prefix = site.getName().replaceAll(" ", "_")+"_curves_"+imt.getPrefix()+"_"+gmpeRef.getShortName();
			
			if (replotCurves || !new File(resourcesDir, prefix+".png").exists())
				curvePlotter.plotHazardCurves(resourcesDir, prefix, imt);
			
			lines.add("### "+imt.getDisplayName()+" Hazard Curves");
			lines.add(topLink); lines.add("");
			lines.add("![Hazard Curve]("+resourcesDir.getName()+"/"+prefix+".png)");
			lines.add("");
			
			if (hasMultiple) {
				lines.add("#### "+imt.getDisplayName()+" Percentile Curves");
				lines.add(topLink); lines.add("");
				lines.add("Here, we calculate hazard curves at various percentiles of the rupture variation "
						+ "distribution for each rupture. For example, the p50 curve shows the hazard curve using only "
						+ "the median ground motion for each rupture from the rupture variation distribution.");
				lines.add("");
				lines.addAll(curvePlotter.getCurveLegend(false, false, false, true, 0));
				lines.add("");
				
				String percentilePrefix = prefix+"_percentiles";
				
				if (replotCurves || !new File(resourcesDir, percentilePrefix+".png").exists())
					curvePlotter.plotSimFractileHazardCurves(resourcesDir, percentilePrefix, imt);
				
				lines.add("![Hazard Curve]("+resourcesDir.getName()+"/"+percentilePrefix+".png)");
				lines.add("");
			}
			
			lines.add("#### "+imt.getDisplayName()+" GMPE-Sim Comparison");
			lines.add(topLink); lines.add("");
			int numGMPESims = 100;
			lines.addAll(curvePlotter.getCurveLegend(false, false, false, false, numGMPESims));
			lines.add("");
			
			String gmpeSimPrefix = prefix+"_gmpe_sims";
			
			if (replotCurves || !new File(resourcesDir, gmpeSimPrefix+".png").exists())
				curvePlotter.plotGMPE_SimHazardCurves(resourcesDir, gmpeSimPrefix, imt, gmpeSimProv, numGMPESims);
			
			lines.add("![Hazard Curve]("+resourcesDir.getName()+"/"+gmpeSimPrefix+".png)");
			lines.add("");
			
			if (sourceRupContributionFracts != null) {
				String sourceContribPrefix = prefix+"_source_contrib";
				
				File simFile, gmpeFile;
				
				if (replotCurves || !new File(resourcesDir, sourceContribPrefix+"_gmpe.png").exists()) {
					System.out.println("Plotting simulation source curves");
					simFile = curvePlotter.plotSourceContributionHazardCurves(resourcesDir, sourceContribPrefix+"_sim", imt,
							sourceRupContributionSortProb, sourceRupContributionNum, false);
					System.out.println("Calculating/plotting GMPE source curves");
					gmpeFile = curvePlotter.plotSourceContributionHazardCurves(resourcesDir, sourceContribPrefix+"_gmpe", imt,
							sourceRupContributionSortProb, sourceRupContributionNum, true);
				} else {
					simFile = new File(resourcesDir, sourceContribPrefix+"_sim.png");
					gmpeFile = new File(resourcesDir, sourceContribPrefix+"_gmpe.png");
				}
				
				if (simFile != null && simFile.exists() && gmpeFile != null && gmpeFile.exists()) {
					lines.add("#### "+imt.getDisplayName()+" Source Contributions");
					lines.add(topLink); lines.add("");
					lines.add("");
					String line = "These plots show the contribution of each fault source to the hazard curves. The same set of "
							+ "sources are plotted for both simulation values (left) and GMPE values (right). Sources are sorted in the legend "
							+ "(and colored by) their ";
					if (sourceRupContributionSortProb > 0)
						line += "contrubution in the simulation results at the **POE="+(float)+sourceRupContributionSortProb+"** level";
					else
						line += "risk-targeted ground motion in the simulation results";
					line += ", and only the top "+sourceRupContributionNum+" sources are plotted.";
					lines.add(line);
					lines.add("");
					if (!sourceRupContributionsMutuallyExclusive) {
						lines.add("*NOTE: Source curves are not mututally exclusive. For the case of multi fault ruptures, "
								+ "a single rupture can be included in the curve for multiple sources*");
						lines.add("");
					}
					
					table = MarkdownUtils.tableBuilder();
					
					table.addLine("**Simulation Source Contributions**", "**GMPE Source Contributions**");
					table.addLine("![Hazard Curve]("+resourcesDir.getName()+"/"+simFile.getName()+")",
							"![Hazard Curve]("+resourcesDir.getName()+"/"+gmpeFile.getName()+")");
					lines.addAll(table.build());
					lines.add("");
				}
			}
			
			double[] timeRange = RuptureComparison.getRuptureTimeRange(comps);
			if (timeRange != null && timeRange[1] > timeRange[0]) {
				// we have time based ruptures, generate an animation
				double duration = timeRange[1] - timeRange[0];
				double delta;
				if (duration > 5000000)
					delta = 500000;
				else if (duration > 2500000)
					delta = 250000;
				else if (duration > 1000000)
					delta = 100000;
				else if (duration > 500000)
					delta = 50000;
				else if (duration > 250000)
					delta = 25000;
				else if (duration > 100000)
					delta = 10000;
				else if (duration > 100000)
					delta = 10000;
				else if (duration > 50000)
					delta = 5000;
				else
					delta = 0d;
				if (delta > 0) {
					lines.add("#### "+imt.getDisplayName()+" Simulation Curve Animation");
					lines.add(topLink); lines.add("");
					lines.add("This animation shows the affect of input simulation catalog length on the simulation hazard curve. "
							+ "Each frame adds an additional "+optionalDigitDF.format(delta)+" years of the catalog to the simulation "
									+ "(in both the simulation and GMPE curves). Previous results (for shorter sub-catalogs) are faded out.");
					lines.add("");
					
					File animFile = new File(resourcesDir, prefix+"_animation.gif");
					
					if (replotCurves || !animFile.exists())
						curvePlotter.plotCurveAnimation(animFile, timeRange, delta, imt);
					
					lines.add("![Hazard Curve Animation]("+resourcesDir.getName()+"/"+animFile.getName()+")");
					lines.add("");
				}
			}
		}
		
		lines.add("## Disaggregations");
		lines.add(topLink); lines.add("");
		
		double minMag = Double.POSITIVE_INFINITY;
		boolean hasSimDist = true;
		for (RuptureComparison<E> comp : comps) {
			minMag = Math.min(minMag, comp.getMagnitude());
			hasSimDist = hasSimDist && simProv.getNumSimulations(site, comp.getRupture()) > 1;
		}
		minMag = Math.floor(minMag*10d)/10d;
		
		lines.add("**Note on Epsilons:** Epsilon values are not straightforward for simulations. For a GMPE (in natural log space):");
		lines.add("");
		lines.add("**GMPE Epsilon:** *epsilon = (gmpe_IML - gmpe_mean)/gmpe_sigma*");
		lines.add("");
		if (hasSimDist) {
			lines.add("This simulation has at least 2 simulations per rupture, so we can calculate a 'simulation' "
					+ "mean/sigma for use in epsilon calculations. These disaggregations are labeled 'w/ sim dist for Epsilon':");
			lines.add("");
			lines.add("**Simulation Distribution Epsilon:** *epsilon = (sim_IML - sim_mean)/sim_sigma*");
			lines.add("");
			lines.add("We also include plots which use the GMPE mean and sigma when computing sigma:");
		} else {
			lines.add("This simulation does not have distributions for each rupture, so in order to compute an epsilon, "
					+ "we must use the GMPE mean and sigma values:");
		}
		lines.add("");
		lines.add("**Simulation w/ GMPE Distribution Epsilon:** *epsilon = (sim_IML - gmpe_mean)/gmpe_sigma*");
		lines.add("");
		
		RuptureComparisonERF<E> disaggERF = new RuptureComparisonERF<>(comps);
		disaggERF.updateForecast();
		
		List<Future<File>> disaggFutures = new ArrayList<>();
		for (int p=0; p<imts.length; p++) {
			IMT imt = imts[p];
			if (imt.getParamName().equals(PGV_Param.NAME)) {
				System.err.println("WARNING: skipping PGV disagg until we figure out exceptions");
				continue;
			}
			DiscretizedFunc simCurve = curvePlotter.getCalcSimCurve(simCalc, imt);
			DiscretizedFunc gmpeCurve = curvePlotter.getCalcGMPECurve(imt);
			
			lines.add("### "+imt.getDisplayName()+" Disaggregations");
			lines.add(topLink); lines.add("");
			
			List<String> disaggHeaders = new ArrayList<>();
			List<Boolean> isProbs = new ArrayList<>();
			List<List<String>> disaggLabels = new ArrayList<>();
			List<List<String>> disaggFileLabels = new ArrayList<>();
			List<List<Double>> disaggVals = new ArrayList<>();
			
			// for intersections
			List<Point2D> intersections = getIntersections(simCurve, gmpeCurve, 1e-6, 1e-2);
			if (!intersections.isEmpty()) {
				disaggHeaders.add("Simulation/GMPE Intersections");
				isProbs.add(false);
				List<String> myDisaggLabels = new ArrayList<>();
				List<String> myDisaggFileLabels = new ArrayList<>();
				List<Double> myDisaggVals = new ArrayList<>();
				disaggLabels.add(myDisaggLabels);
				disaggFileLabels.add(myDisaggFileLabels);
				disaggVals.add(myDisaggVals);
				for (Point2D pt : intersections) {
					double prob = pt.getY();
					double iml = pt.getX();

					int rp = (int)(Math.round(1/prob));

					myDisaggLabels.add(rp+" yr<br>"+(float)iml+" g");
					myDisaggFileLabels.add("intersect_"+(float)iml);
					myDisaggVals.add(iml);
				}
			}
			
			// for return periods
			if (MultiRupGMPE_ComparePageGen.hazard_curve_rps.length > 0) {
				disaggHeaders.add("Fixed Return Periods");
				isProbs.add(true);
				List<String> myDisaggLabels = new ArrayList<>();
				List<String> myDisaggFileLabels = new ArrayList<>();
				List<Double> myDisaggVals = new ArrayList<>();
				disaggLabels.add(myDisaggLabels);
				disaggFileLabels.add(myDisaggFileLabels);
				disaggVals.add(myDisaggVals);
				for (int rp : MultiRupGMPE_ComparePageGen.hazard_curve_rps) {
					myDisaggLabels.add(rp+" yr");
					myDisaggFileLabels.add(rp+"yr");
					myDisaggVals.add(1d/(double)rp);
				}
			}
			
			// for fixed locations
			if (disagg_fixed_accel_vals.length > 0 &&
					(imt.getParamName().equals(SA_Param.NAME) || imt.getParamName().equals(PGA_Param.NAME))) {
				disaggHeaders.add("Fixed IMLs");
				isProbs.add(false);
				List<String> myDisaggLabels = new ArrayList<>();
				List<String> myDisaggFileLabels = new ArrayList<>();
				List<Double> myDisaggVals = new ArrayList<>();
				disaggLabels.add(myDisaggLabels);
				disaggFileLabels.add(myDisaggFileLabels);
				disaggVals.add(myDisaggVals);
				for (double iml : disagg_fixed_accel_vals) {
					myDisaggLabels.add((float)iml+" g");
					myDisaggFileLabels.add("iml_"+(float)iml);
					myDisaggVals.add(iml);
				}
			}
			
			Preconditions.checkState(disaggHeaders.size() == disaggLabels.size());
			Preconditions.checkState(disaggHeaders.size() == isProbs.size());
			Preconditions.checkState(disaggHeaders.size() == disaggFileLabels.size());
			Preconditions.checkState(disaggHeaders.size() == disaggVals.size());
			
			for (int i=0; i<disaggHeaders.size(); i++) {
				List<String> myLabels = disaggLabels.get(i);
				List<String> myFileLabels = disaggFileLabels.get(i);
				List<Double> myVals = disaggVals.get(i);
				Preconditions.checkState(myLabels.size() ==  myFileLabels.size());
				Preconditions.checkState(myLabels.size() ==  myVals.size());
				if (myLabels.isEmpty())
					continue;
				boolean isProb = isProbs.get(i);
				
				lines.add("#### "+imt.getDisplayName()+" Disaggregations at "+disaggHeaders.get(i));
				lines.add(topLink); lines.add("");
				
				table = MarkdownUtils.tableBuilder();
				
				table.initNewLine();
				
				table.addColumn("**Disagg Level**");
				
				if (hasSimDist) {
					table.addColumn("**"+simProv.getName()+" w/ sim dist for Epsilon**");
					table.addColumn("**"+simProv.getName()+" w/ GMPE dist for Epsilon**");
				} else {
					table.addColumn("**"+simProv.getName()+" w/ GMPE dist for Epsilon**");
				}
				table.addColumn("**"+gmpeRef.getName()+"**");
				table.finalizeLine();
				
				for (int j=0; j<myLabels.size(); j++) {
					table.initNewLine().addColumn("**"+myLabels.get(j)+"**");
					
					double val = myVals.get(j);
					
					double iml = isProb ? HazardDataSetLoader.getCurveVal(simCurve, false, val) : val;
					double disaggIML = iml;
					if (imt.getParamName().equals(SA_Param.NAME) || imt.getParamName().equals(PGA_Param.NAME))
						// need to disagg on the log of the IML
						disaggIML = Math.log(disaggIML);
					if (!Double.isFinite(iml) || iml > simCurve.getMaxX()) {
						System.out.println("Couldn't get IML for simulation, "+myLabels.get(j)+". Skipping disagg!");
						table.addColumn("N/A");
						table.addColumn("N/A");
					} else {
						if (hasSimDist) {
							String prefix = "disagg_sim_"+imt.getPrefix()+"_"+myFileLabels.get(j);
							ScalarIMR gmpe = new SimulationDisaggAttenuationRelationshipWrapper<>(
									simProv, Source.SIMULATION, Source.SIMULATION, Source.SIMULATION, imts);
							if (replotDisaggs || !new File(resourcesDir, prefix+".png").exists())
								disaggFutures.add(getExec().submit(new DisaggCallable(gmpe, site, disaggERF, disaggIML, minMag, resourcesDir, prefix)));
							table.addColumn("![Disaggregation]("+resourcesDir.getName()+"/"+prefix+".png)");
						}
						
						String prefix = "disagg_sim_gmpe_dist_for_epsilon_"+imt.getPrefix()+"_"+myFileLabels.get(j);
						ScalarIMR gmpe = new SimulationDisaggAttenuationRelationshipWrapper<>(
								simProv, Source.SIMULATION, Source.GMPE, Source.GMPE, imts);
						if (replotDisaggs || !new File(resourcesDir, prefix+".png").exists())
							disaggFutures.add(getExec().submit(new DisaggCallable(gmpe, site, disaggERF, disaggIML, minMag, resourcesDir, prefix)));
						table.addColumn("![Disaggregation]("+resourcesDir.getName()+"/"+prefix+".png)");
					}
					
					iml = isProb ? HazardDataSetLoader.getCurveVal(gmpeCurve, false, val) : val;
					if (!Double.isFinite(iml)|| iml > gmpeCurve.getMaxX()) {
						System.out.println("Couldn't get IML for GMPE, "+myLabels.get(j)+". Skipping disagg!");
						table.addColumn("N/A");
					} else {
						String prefix = "disagg_gmpe_"+imt.getPrefix()+"_"+myFileLabels.get(j);
						ScalarIMR gmpe = new SimulationDisaggAttenuationRelationshipWrapper<>(
								simProv, Source.GMPE, Source.GMPE, Source.GMPE, imts);
						if (replotDisaggs || !new File(resourcesDir, prefix+".png").exists())
							disaggFutures.add(getExec().submit(new DisaggCallable(gmpe, site, disaggERF, disaggIML, minMag, resourcesDir, prefix)));
						table.addColumn("![Disaggregation]("+resourcesDir.getName()+"/"+prefix+".png)");
					}
					
					table.finalizeLine();
				}
				lines.addAll(table.build());
				lines.add("");
			}
		}
		
		System.out.println("Waiting on "+disaggFutures.size()+" disagg futures");
		for (Future<File> f : disaggFutures) {
			try {
				f.get();
			} catch (Exception e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
		System.out.println("DONE with disagg");
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 4));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private List<Point2D> getIntersections(DiscretizedFunc curve1, DiscretizedFunc curve2, double minY, double maxY) {
		List<Point2D> intersections = new ArrayList<>();
		
		Range yRange = new Range(minY, maxY);
		
		for (int i=1; i<curve1.size() && i<curve2.size(); i++) {
			double x1 = curve1.getX(i-1);
			double x2 = curve1.getX(i);
			double x3 = curve2.getX(i-1);
			double x4 = curve2.getX(i);
			Preconditions.checkState(x1 == x3);
			Preconditions.checkState(x2 == x4);
			
			double y1 = curve1.getY(i-1);
			double y2 = curve1.getY(i);
			double y3 = curve2.getY(i-1);
			double y4 = curve2.getY(i);
			
			if (!yRange.contains(y1) || !yRange.contains(y2) || !yRange.contains(y3) || !yRange.contains(y4))
				continue;
			
			if (y1 <= y3 && y2 >= y4 || y1 >= y3 && y2 <= y4) {
				double px = ((x1*y2 - y1*x2)*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4))/((x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4));
				double py = ((x1*y2 - y1*x2)*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4))/((x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4));
				intersections.add(new Point2D.Double(px, py));
			}
		}
		
		return intersections;
	}
	
	public synchronized static ExecutorService getExec() {
		if (exec == null)
			exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		return exec;
	}
	
	protected static ScalarIMR checkOutGMPE(AttenRelRef gmpeRef) {
		synchronized (SiteHazardCurveComarePageGen.class) {
			if (gmpesInstancesCache == null)
				gmpesInstancesCache = new HashMap<>();
		}
		synchronized (gmpesInstancesCache) {
			LinkedList<ScalarIMR> gmpes = gmpesInstancesCache.get(gmpeRef);
			if (gmpes == null) {
				gmpes = new LinkedList<>();
				gmpesInstancesCache.put(gmpeRef, gmpes);
			}
			if (!gmpes.isEmpty())
				return gmpes.pop();
		}
		ScalarIMR gmpe = gmpeRef.instance(null);
		gmpe.setParamDefaults();
		gmpe.setIntensityMeasure(SA_Param.NAME);
		return gmpe;
	}
	
	protected static void checkInGMPE(AttenRelRef gmpeRef, ScalarIMR gmpe) {
		synchronized (gmpesInstancesCache) {
			gmpesInstancesCache.get(gmpeRef).push(gmpe);
		}
	}
	
	public class DisaggCallable implements Callable<File> {
		
		private ScalarIMR gmpe;
		private Site site;
		private AbstractERF erf;
		private double iml;
		private double minMag;
		private File outputDir;
		private String prefix;

		public DisaggCallable(ScalarIMR gmpe, Site site, AbstractERF erf, double iml, double minMag, File outputDir, String prefix) {
			this.gmpe = gmpe;
			this.site = site;
			this.erf = erf;
			this.iml = iml;
			this.minMag = minMag;
			this.outputDir = outputDir;
			this.prefix = prefix;
		}

		@Override
		public File call() {
			DisaggCalc disaggCalc = checkOutDisaggCalc(minMag);
			
			File file;
			try {
				file = disaggCalc.disagg(gmpe, site, erf, iml, outputDir, prefix);
			} catch (IOException e) {
				System.err.println("Disagg failed for iml="+iml+", prefix="+prefix);
				throw ExceptionUtils.asRuntimeException(e);
			}
			
			checkInDisaggCalc(disaggCalc);
			
			return file;
		}
		
	}
	
	protected DisaggCalc checkOutDisaggCalc(double minMag) {
		synchronized (SiteHazardCurveComarePageGen.class) {
			if (disaggCalcsCache == null)
				disaggCalcsCache = new LinkedList<>();
		}
		synchronized (disaggCalcsCache) {
			if (!disaggCalcsCache.isEmpty())
				return disaggCalcsCache.pop();
		}
		return new DisaggCalc(minMag);
	}
	
	protected void checkInDisaggCalc(DisaggCalc disaggCalc) {
		synchronized (disaggCalcsCache) {
			disaggCalcsCache.push(disaggCalc);
		}
	}
	
	private static long prevDisaggMillis;
	private synchronized static void checkSleepDisagg() {
		long curTime = System.currentTimeMillis();
		long diff = curTime - prevDisaggMillis;
		if (diff < 100) {
			try {
				long sleepTime = 100 - diff;
				System.out.println("Sleeping for "+sleepTime+" milliseconds");
				Thread.sleep(sleepTime);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		prevDisaggMillis = curTime;
	}
	
	private class DisaggCalc {
		private final ParameterList disaggParams;
		private final List<SourceFilter> sourceFilters;
		
		// disagg plot settings
		private final double minMag;
		private final int numMags;
		private final double deltaMag;

		private final int numSourcesForDisag = 100;

		private final boolean showSourceDistances = false;

		private final double maxZAxis = Double.NaN;
		
		private DisaggregationCalculatorAPI disaggCalc;
		
		public DisaggCalc(double minMag) {
			this.minMag = minMag;
			this.deltaMag = 0.2;
			this.numMags = (int)((8.6d - minMag)/deltaMag + 0.5);
			
			disaggParams = DisaggregationCalculator.getDefaultParams();
			sourceFilters = SourceFiltersParam.getDefault().getEnabledFilters();
			
			disaggCalc = new DisaggregationCalculator();
			disaggCalc.setMagRange(minMag, numMags, deltaMag);
			disaggCalc.setNumSourcesToShow(numSourcesForDisag);
			disaggCalc.setShowDistances(showSourceDistances);
		}
		
		public File disagg(ScalarIMR gmpe, Site site, AbstractERF erf, double iml, File outputDir, String prefix) throws IOException {
			boolean success = disaggCalc.disaggregate(iml, site, gmpe, erf, sourceFilters, disaggParams);
			if (!success)
				throw new RuntimeException("Disagg calc failed for iml="+iml+", prefix='"+prefix+"'(see errors above, if any).");
			disaggCalc.setMaxZAxisForPlot(maxZAxis);
			System.out.println("Done Disaggregating");
			String metadata = "temp metadata";

			System.out.println("Fetching plot...");
			String address;
			checkSleepDisagg();
			address = disaggCalc.getDisaggregationPlotUsingServlet(metadata);

//			String meanModeText = disaggCalc.getMeanAndModeInfo();
//			String binDataText = disaggCalc.getBinData();
//			String sourceDataText = disaggCalc.getDisaggregationSourceInfo();

			File outputFile = new File(outputDir, prefix);

//			String metadataText = "Custom disagg";

			File pngFile = new File(outputFile.getAbsolutePath()+".png");
//			File pdfFile = new File(outputFile.getAbsolutePath()+".pdf");
//			DisaggregationPlotViewerWindow.saveAsPDF(
//					address+DisaggregationCalculator.DISAGGREGATION_PLOT_PDF_NAME,
//					pdfFile.getAbsolutePath(), meanModeText, metadataText, binDataText, sourceDataText);
			FileUtils.downloadURL(address+DisaggregationCalculator.DISAGGREGATION_PLOT_PNG_NAME,
					pngFile);
//			DisaggregationPlotViewerWindow.saveAsTXT(outputFile.getAbsolutePath()+".txt", meanModeText, metadataText,
//					binDataText, sourceDataText);
			return pngFile;
		}
	}
	
	static final DecimalFormat optionalDigitDF = new DecimalFormat("0.##");
	
	private static final String GMAPS_API = new String(Base64.getDecoder().decode(
			"QUl6YVN5QTE5Y0lWOEJxQTdPdXskld58190u50hCSC1OYzdTcDE0YTJLbktKX19j".replaceAll("skld58190u50", "")));
	private static final String MAPQUEST_API = new String(Base64.getDecoder().decode(
			"NTNBMDVBskld58190u50eThZT2t5eTdJYkV4TFVRWHNtMVskld58190u50VWTVVxdmo=".replaceAll("skld58190u50", "")));
	
	public static String getMiniMap(Location loc) {
		String locStr = loc.getLatitude()+","+loc.getLongitude();
//		return "https://maps.googleapis.com/maps/api/staticmap?center="+locStr+"&zoom=9"
//				+ "&size=400x300&maptype=roadmap&markers="+locStr+"&key="+GMAPS_API;
		return "https://www.mapquestapi.com/staticmap/v5/map?key="+MAPQUEST_API
				+"&size=400,300&zoom=9&center="+locStr+"&locations="+locStr+"|marker-end&imagetype=png";
	}
}