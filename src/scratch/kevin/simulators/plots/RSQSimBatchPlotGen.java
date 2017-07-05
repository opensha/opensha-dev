package scratch.kevin.simulators.plots;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.dom4j.DocumentException;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.util.ClassUtils;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.MagRangeRuptureIdentifier;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.parsers.RSQSimFileReader;

import com.google.common.base.Preconditions;
import com.google.common.base.Splitter;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.griddedSeismicity.GridSourceProvider;
import scratch.UCERF3.utils.FaultSystemIO;

/**
 * TODO:
 * * Multi fault participation
 * * map based recurrence and ratios
 * @author kevin
 *
 */
public class RSQSimBatchPlotGen {
	
	private static final String DEFAULT_NAME = "RSQSim Catalog";
	
	private static final double u3MinMag = 5.05;
	private static final double u3MaxMag = 9.05;
	private static final double u3Delta = 0.1;
	private static final int u3NumMag = (int)((u3MaxMag - u3MinMag)/u3Delta + 0.5);
	
	private static final String ri_mags_default = "6,7,7.5,8";
	
	public enum PlotConfig {
		TOTAL_MFD("mfd", "total-mfd", false, "Plot total regional MFD", null) {
			@Override
			protected List<? extends AbstractPlot> buildPlots(String catalogName, File outputDir,
					Double minMag, FaultSystemSolution u3Sol, String arg, ElementBundles elemBundle) {
				MFDPlot regionMFD = new MFDPlot(minMag);
				if (u3Sol != null) {
					IncrementalMagFreqDist u3MFD = calcU3TotalMFD(u3Sol, false);
					regionMFD.setComparableMFD(u3MFD, "UCERF3 On-Fault");
				}
				regionMFD.initialize(catalogName, outputDir, "mfd_total");
				return Lists.newArrayList(regionMFD);
			}
		},
		PARENT_SECT_MFD("pmfd", "parent-sect-mfds", false, "Plot parent section participation MFDs", null) {
			@Override
			protected List<? extends AbstractPlot> buildPlots(String catalogName, File outputDir,
					Double minMag, FaultSystemSolution u3Sol, String arg, ElementBundles elemBundle) {
				outputDir = new File(outputDir, "parent_sect_partic");
				Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
				Map<Integer, List<SimulatorElement>> faults = elemBundle.getFaultBundledElems();
				if (faults.isEmpty()) {
					System.out.println("No fault info, skipping sub-sect-mfd plots");
					return null;
				}
				List<MFDPlot> plots = Lists.newArrayList();
				for (Integer faultID : faults.keySet()) {
					String name = elemBundle.getFaultName(faultID);
					MFDPlot plot = new MFDPlot(minMag, faults.get(faultID));
					plot.setPlotTitle(name+" MFD");
					if (u3Sol != null) {
						int parentSectID = elemBundle.getParentSectIDforFault(faultID);
						IncrementalMagFreqDist u3MFD = u3Sol.calcParticipationMFD_forParentSect(parentSectID, u3MinMag, u3MaxMag, u3NumMag);
						plot.setComparableMFD(u3MFD, "UCERF3 Parent Sect");
					}
					plot.initialize(catalogName, outputDir, getFileSafe(name));
					plots.add(plot);
				}
				return plots;
			}
		},
		TOTAL_RI("tri", "total-ri", true, "Total recurrence interval distribution plot", ri_mags_default) {
			@Override
			protected List<? extends AbstractPlot> buildPlots(String catalogName, File outputDir, Double minMag,
					FaultSystemSolution u3Sol, String arg, ElementBundles elemBundle) {
				double[] minMags = commaDoubleSplit(arg);
				RecurrenceIntervalPlot riPlot = new RecurrenceIntervalPlot(minMags);
				if (u3Sol != null) {
					double[] compVals = new double[minMags.length];
					IncrementalMagFreqDist u3MFD = calcU3TotalMFD(u3Sol, false);
					EvenlyDiscretizedFunc u3Cumulative = u3MFD.getCumRateDistWithOffset();
					for (int i=0; i<minMags.length; i++)
						compVals[i] = 1d/u3Cumulative.getInterpolatedY_inLogYDomain(minMags[i]);
					riPlot.setComparison(compVals, "UCERF3 On-Fault");
				}
				riPlot.initialize(catalogName, outputDir, "interevent_total");
				return Lists.newArrayList(riPlot);
			}
		},
		RARENT_SECT_RI("pri", "parent-sect-ris", true, "Parent section recurrence interval distribution plots", ri_mags_default) {
			@Override
			protected List<? extends AbstractPlot> buildPlots(String catalogName, File outputDir, Double minMag,
					FaultSystemSolution u3Sol, String arg, ElementBundles elemBundle) {
				outputDir = new File(outputDir, "parent_sect_interevent");
				Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
				
				double[] minMags = commaDoubleSplit(arg);
				
				Map<Integer, List<SimulatorElement>> faults = elemBundle.getFaultBundledElems();
				if (faults.isEmpty()) {
					System.out.println("No fault info, skipping sub-sect-mfd plots");
					return null;
				}
				
				List<RecurrenceIntervalPlot> riPlots = Lists.newArrayList();
				for (Integer faultID : faults.keySet()) {
					String name = elemBundle.getFaultName(faultID);
					RecurrenceIntervalPlot riPlot = new RecurrenceIntervalPlot(faults.get(faultID), minMags);
					riPlot.setPlotTitle(name);
					if (u3Sol != null) {
						int parentSectID = elemBundle.getParentSectIDforFault(faultID);
						IncrementalMagFreqDist u3MFD = u3Sol.calcParticipationMFD_forParentSect(parentSectID, u3MinMag, u3MaxMag, u3NumMag);
						EvenlyDiscretizedFunc u3Cumulative = u3MFD.getCumRateDistWithOffset();
						double[] compVals = new double[minMags.length];
						for (int i=0; i<minMags.length; i++)
							compVals[i] = 1d/u3Cumulative.getInterpolatedY_inLogYDomain(minMags[i]);
						riPlot.setComparison(compVals, "UCERF3 Parent Sect");
					}
					riPlot.initialize(catalogName, outputDir, getFileSafe(name));
					riPlots.add(riPlot);
				}
				return riPlots;
			}
		},
		SUB_SECT_RI_SUMMARY("ssri", "sub-sect-ris", true, "Sub section RI CSV file and scatter (if U3 values supplied)", ri_mags_default) {
			@Override
			protected List<? extends AbstractPlot> buildPlots(String catalogName, File outputDir, Double minMag,
					FaultSystemSolution u3Sol, String arg, ElementBundles elemBundle) {
				double[] minMags = commaDoubleSplit(arg);
				SubSectRecurrenceSummary riPlot = new SubSectRecurrenceSummary(elemBundle.getElements(), minMags);
				if (u3Sol != null) {
					List<EvenlyDiscretizedFunc> sectCumMFDs = Lists.newArrayList();
					for (int i=0; i<u3Sol.getRupSet().getNumSections(); i++)
						sectCumMFDs.add(u3Sol.calcParticipationMFD_forSect(i, u3MinMag, u3MaxMag, u3NumMag).getCumRateDistWithOffset());
					riPlot.setComparison(sectCumMFDs, "UCERF3 On-Fault");
				}
				riPlot.initialize(catalogName, outputDir, "interevent_sub_sects");
				return Lists.newArrayList(riPlot);
			}
		},
		MAG_AREA_SCALING("mas", "mag-area-scaling", false, "Magnitude area scaling plot.", null) {
			@Override
			protected List<? extends AbstractPlot> buildPlots(String catalogName, File outputDir, Double minMag,
					FaultSystemSolution u3Sol, String arg, ElementBundles elemBundle) {
				MagAreaScalingPlot plot = new MagAreaScalingPlot();
				plot.initialize(catalogName, outputDir, "mag_area_scaling");
				return Lists.newArrayList(plot);
			}
		};
		
		private Option op;
		private String defaultArg;
		private PlotConfig(String shortOpt, String longOpt, boolean hasArg, String description, String defaultArg) {
			op = new Option(shortOpt, longOpt, hasArg, description);
			op.setRequired(false);
			this.defaultArg = defaultArg;
		}
		
		protected abstract List<? extends AbstractPlot> buildPlots(
				String catalogName, File outputDir, Double minMag, FaultSystemSolution u3Sol, String arg, ElementBundles elemBundle);
		
		public List<? extends AbstractPlot> buildIfapplicable(CommandLine cmd, boolean forcePlot,
				String catalogName, File outputDir, Double minMag, FaultSystemSolution u3Sol, ElementBundles elemBundle) {
			boolean hasOption = cmd.hasOption(op.getLongOpt());
			if (!hasOption && !forcePlot)
				return null;
			String arg = defaultArg;
			if (hasOption && op.hasArg())
				arg = cmd.getOptionValue(op.getLongOpt());
			return buildPlots(catalogName, outputDir, minMag, u3Sol, arg, elemBundle);
		}
	}
	/*
	 * TODO:
	 * * average slip vs mag (can get slip from Mw if needed)
	 */
	
	public static Options createOptions() {
		Options ops = new Options();

		// global options/requirements
		Option geomFile = new Option("g", "geometry-file", true, "Geometry file");
		geomFile.setRequired(true);
		ops.addOption(geomFile);

		Option catalogFile = new Option("c", "catalog-file", true, "Path to one of the .*List files, or a directory containing them");
		catalogFile.setRequired(true);
		ops.addOption(catalogFile);
		
		Option outputDir = new Option("o", "output-dir", true, "Output directory where plots should be written");
		outputDir.setRequired(true);
		ops.addOption(outputDir);
		
		Option catalogTitle = new Option("n", "name", true, "Short name for this catalog for use in plots (Default: "+DEFAULT_NAME+")");
		catalogTitle.setRequired(false);
		ops.addOption(catalogTitle);
		
		Option minMag = new Option("m", "min-mag", true, "Minimum magnitude to load/plot (default is load all)");
		minMag.setRequired(false);
		ops.addOption(minMag);
		
		Option u3File = new Option("u", "ucerf-sol", true, "UCERF3 Fault System Solution file for comparisons");
		u3File.setRequired(false);
		ops.addOption(u3File);
		
		Option plotAll = new Option("a", "plot-all", false, "Flag to generate all plots");
		plotAll.setRequired(false);
		ops.addOption(plotAll);
		
		Option skipYears = new Option("s", "skip-years", true,
				"Time in years to skip from the beginning of the catalog, relative to the first observed event time");
		skipYears.setRequired(false);
		ops.addOption(skipYears);
		
		// add individual plots
		for (PlotConfig plot : PlotConfig.values())
			ops.addOption(plot.op);

		return ops;
	}
	
	protected static CommandLine parse(Options options, String args[], Class<?> clazz) {
		try {
			CommandLineParser parser = new GnuParser();
			
			CommandLine cmd = parser.parse(options, args);
			return cmd;
		} catch (Exception e) {
			System.out.println(e.getMessage());
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(
					ClassUtils.getClassNameWithoutPackage(clazz),
					options, true );
			System.exit(2);
			return null; // not accessible
		}
	}

	public static void main(String[] args) throws IOException, DocumentException {
		if (args.length == 1 && args[0].equals("--hardcoded")) {
			System.out.println("Hardcoded test!");
			File dir = new File("/home/kevin/Simulators/UCERF3_JG_supraSeisGeo2");
			File geomFile = new File(dir, "UCERF3.D3.1.1km.tri.2.flt");
			File outputDir = new File("/tmp/rsqsim_plots");
			File solFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
					+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
			String argStr = "--geometry-file "+geomFile.getAbsolutePath()+" --catalog-file "+dir.getAbsolutePath()
					+" --output-dir "+outputDir.getAbsolutePath()+" --name TestCatalog"
					+" --ucerf-sol "+solFile.getAbsolutePath();
			argStr += " --plot-all --min-mag 4";
//			argStr += " --mag-area-scaling --min-mag 4";
			argStr += " --skip-years 10000";
			args = Splitter.on(" ").splitToList(argStr).toArray(new String[0]);
		}
		
		Options options = createOptions();
		
		CommandLine cmd = parse(options, args, RSQSimBatchPlotGen.class);
		
		File geomFile = new File(cmd.getOptionValue("geometry-file"));
		Preconditions.checkState(geomFile.exists(), "Geometry file doesn't exist: %s", geomFile.getAbsolutePath());
		System.out.println("Loading geometry");
		List<SimulatorElement> elements = RSQSimFileReader.readGeometryFile(geomFile, 11, 'S');
		
		Double minMag = null;
		List<RuptureIdentifier> loadIdens = null;
		if (cmd.hasOption("min-mag")) {
			minMag = Double.parseDouble(cmd.getOptionValue("min-mag"));
			loadIdens = Lists.newArrayList();
			loadIdens.add(new MagRangeRuptureIdentifier(minMag, Double.POSITIVE_INFINITY));
		}
		
		double skipYears = 0d;
		if (cmd.hasOption("skip-years"))
			skipYears = Double.parseDouble(cmd.getOptionValue("skip-years"));
		
		File catalogFile = new File(cmd.getOptionValue("catalog-file"));
		Preconditions.checkState(catalogFile.exists(), "Catalog file/dir doesn't exist: %s", catalogFile.getAbsolutePath());
		System.out.println("Creating event file iterator");
		Iterable<RSQSimEvent> eventsIterable = RSQSimFileReader.getEventsIterable(catalogFile, elements, loadIdens);
		
		File outputDir = new File(cmd.getOptionValue("output-dir"));
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir(),
				"Output dir doesn't exist or couldn't be created: %s", outputDir.getAbsolutePath());
		
		String catalogName = DEFAULT_NAME;
		if (cmd.hasOption("name"))
			catalogName = cmd.getOptionValue("name");
		
		FaultSystemSolution u3Sol = null;
		FaultSystemRupSet rupSet = null;
		if (cmd.hasOption("ucerf-sol")) {
			u3Sol = FaultSystemIO.loadSol(new File(cmd.getOptionValue("ucerf-sol")));
			rupSet = u3Sol.getRupSet();
		}
		
		System.out.println("Building section/fault element bundles");
		ElementBundles elemBundle = new ElementBundles(elements, rupSet);
		
		List<List<? extends AbstractPlot>> plots = Lists.newArrayList();
		List<PlotConfig> configs = Lists.newArrayList();
		// plots that are applicable to all elements
		List<AbstractPlot> allElementPlots = Lists.newArrayList();
		// element specific plots
		Map<SimulatorElement, List<AbstractPlot>> elementSpecificPlots = Maps.newHashMap();
		for (SimulatorElement elem : elements)
			elementSpecificPlots.put(elem, new ArrayList<AbstractPlot>());
		
		System.out.println("Initializing plots");
		
		boolean plotAll = cmd.hasOption("plot-all");
		
		for (PlotConfig conf : PlotConfig.values()) {
			List<? extends AbstractPlot> myPlots = conf.buildIfapplicable(cmd, plotAll, catalogName, outputDir, minMag, u3Sol, elemBundle);
			if (myPlots != null && !myPlots.isEmpty()) {
				plots.add(myPlots);
				configs.add(conf);
				for (AbstractPlot plot : myPlots) {
					Collection<SimulatorElement> applicableElems = plot.getApplicableElements();
					if (applicableElems == null)
						allElementPlots.add(plot);
					else
						for (SimulatorElement elem : applicableElems)
							elementSpecificPlots.get(elem).add(plot);
				}
			}
		}
		
		System.out.println("Loading/processing events");
		
		DecimalFormat countDF = new DecimalFormat("#");
		countDF.setGroupingUsed(true);
		countDF.setGroupingSize(3);
		
		int count = 0;
		int mod = 1000;
		int nextMod = mod*10;
		double firstEventTimeYears = Double.NaN;
		double skipPeriodEnd = 0d;
		int numSkipped = 0;
		for (SimulatorEvent e : eventsIterable) {
			if (skipYears > 0) {
				if (Double.isNaN(firstEventTimeYears)) {
					firstEventTimeYears = e.getTimeInYears();
					skipPeriodEnd = firstEventTimeYears + skipYears;
					System.out.println("Skipping "+(float)skipYears+" years. "
							+ "Catalog starts at "+(float)firstEventTimeYears+", will start at year "+(float)skipPeriodEnd);
				}
				if (e.getTimeInYears() < skipPeriodEnd) {
					numSkipped++;
					continue;
				}
			}
			if (count == 0 && numSkipped > 0)
				System.out.println("Skipped "+numSkipped+" events at start of catalog");
			if (count == 0)
				System.out.println("First event processed at "+(float)e.getTimeInYears()+" years");
			
			// first process plots specific to all elements
			for (AbstractPlot plot : allElementPlots)
				plot.processEvent(e);
			// now process plots specific to the elements of this event
			HashSet<AbstractPlot> elemPlots = new HashSet<AbstractPlot>();
			for (EventRecord rec : e)
				for (SimulatorElement elem : rec.getElements())
					for (AbstractPlot plot : elementSpecificPlots.get(elem))
						elemPlots.add(plot);
			for (AbstractPlot plot : elemPlots)
				plot.processEvent(e);
			count++;
			if (count >= nextMod) {
				mod = nextMod;
				nextMod = mod*10;
			}
			if (count % mod == 0)
				System.out.println("Processed "+countDF.format(count)+" events");
		}
		System.out.println("Done processing "+count+" events");
		
		System.out.println("Finalizing plots");
		for (int i=0; i<plots.size(); i++) {
			PlotConfig conf = configs.get(i);
			System.out.println("Writing plot(s): "+conf);
			for (AbstractPlot plot : plots.get(i)) {
				try {
					plot.finalize();
				} catch (Exception e1) {
					System.err.println("Error processing plot, skipping");
					e1.printStackTrace();
				}
			}
		}
		
		System.out.println("DONE");
	}
	
	private static class ElementBundles {
		
		private List<SimulatorElement> elements;
		
		private int minSectIndex;
		private Map<Integer, List<SimulatorElement>> sectBundledElems;
		private Map<Integer, String> sectNamesMap;
		private Map<Integer, List<SimulatorElement>> faultBundledElems;
		private Map<Integer, String> faultNamesMap;
		private Map<Integer, Integer> faultIDtoParentSectMap;
		
		public ElementBundles(List<SimulatorElement> elements, FaultSystemRupSet rupSet) {
			this.elements = elements;
			
			minSectIndex = Integer.MAX_VALUE;
			for (SimulatorElement elem : elements) {
				Integer sectID = elem.getSectionID();
				if (sectID < minSectIndex)
					minSectIndex = sectID;
			}
			
			sectBundledElems = Maps.newHashMap();
			sectNamesMap = Maps.newConcurrentMap();
			faultBundledElems = Maps.newHashMap();
			faultNamesMap = Maps.newHashMap();
			faultIDtoParentSectMap = Maps.newHashMap();
			for (SimulatorElement elem : elements) {
				Integer sectID = elem.getSectionID();
				Integer faultID = elem.getFaultID();
				if (sectID < minSectIndex)
					minSectIndex = sectID;
				if (sectID >= 0) {
					List<SimulatorElement> elemsForSect = sectBundledElems.get(sectID);
					if (elemsForSect == null) {
						elemsForSect = Lists.newArrayList();
						sectBundledElems.put(sectID, elemsForSect);
						sectNamesMap.put(sectID, elem.getSectionName());
					}
					elemsForSect.add(elem);
					if (faultID >= 0 && rupSet != null) {
						int parentSectID = rupSet.getFaultSectionData(sectID - minSectIndex).getParentSectionId();
						if (faultIDtoParentSectMap.containsKey(faultID))
							Preconditions.checkState(faultIDtoParentSectMap.get(faultID) == parentSectID);
						else
							faultIDtoParentSectMap.put(faultID, parentSectID);
					}
				}
				if (faultID >= 0) {
					List<SimulatorElement> elemsForFault = faultBundledElems.get(faultID);
					if (elemsForFault == null) {
						elemsForFault = Lists.newArrayList();
						faultBundledElems.put(faultID, elemsForFault);
						String faultName = elem.getSectionName();
						if (faultName.contains("Subsection"))
							faultName = faultName.substring(0, faultName.indexOf("Subsection")).trim();
						while (faultName.endsWith(","))
							faultName = faultName.substring(0, faultName.length()-1);
						faultNamesMap.put(faultID, faultName);
					}
					elemsForFault.add(elem);
				}
			}
		}
		
		public List<SimulatorElement> getElements() {
			return elements;
		}
		
		public Map<Integer, List<SimulatorElement>> getSectBundledElems() {
			return sectBundledElems;
		}
		
		public String getSectionName(int sectID) {
			return sectNamesMap.get(sectID);
		}
		
		public int getSubsectionIndex(int sectID) {
			return sectID - minSectIndex;
		}
		
		public Map<Integer, List<SimulatorElement>> getFaultBundledElems() {
			return faultBundledElems;
		}
		
		public String getFaultName(int faultID) {
			return faultNamesMap.get(faultID);
		}
		
		public int getParentSectIDforFault(int faultID) {
			Integer parentID = faultIDtoParentSectMap.get(faultID);
			if (parentID == null)
				return -1;
			return parentID;
		}
	}
	
	private static String getFileSafe(String str) {
		return str.replaceAll("\\W+", "_");
	}
	
	private static HashSet<Integer> getElemIndexesSet(Collection<SimulatorElement> elements) {
		HashSet<Integer> indexes = new HashSet<Integer>();
		for (SimulatorElement elem : elements)
			indexes.add(elem.getID());
		return indexes;
	}
	
	
	private static double[] commaDoubleSplit(String str) {
		String[] strs = str.split(",");
		double[] ret = new double[strs.length];
		for (int i=0; i<strs.length; i++)
			ret[i] = Double.parseDouble(strs[i]);
		return ret;
	}

	private static IncrementalMagFreqDist u3OnFaultMFD, u3TotalMFD;
	private static synchronized IncrementalMagFreqDist calcU3TotalMFD(FaultSystemSolution u3Sol, boolean offFault) {
		if (offFault) {
			if (u3TotalMFD == null)
				u3TotalMFD = calcU3TotalMFD(u3Sol, u3MinMag, u3MaxMag, u3Delta, true);
			return u3TotalMFD;
		} else {
			if (u3OnFaultMFD == null)
				u3OnFaultMFD = calcU3TotalMFD(u3Sol, u3MinMag, u3MaxMag, u3Delta, false);
			return u3OnFaultMFD;
		}
	}

	private static IncrementalMagFreqDist calcU3TotalMFD(FaultSystemSolution u3Sol, double u3MinMag, double u3MaxMag,
			double u3Delta, boolean offFault) {
		IncrementalMagFreqDist u3MFD = u3Sol.calcTotalNucleationMFD(u3MinMag, u3MaxMag, u3Delta);
		if (u3Sol.getGridSourceProvider() != null) {
			GridSourceProvider gridProv = u3Sol.getGridSourceProvider();
			for (int i=0; i<gridProv.size(); i++) {
				IncrementalMagFreqDist nodeMFD;
				if (offFault)
					nodeMFD = gridProv.getNodeMFD(i);
				else
					nodeMFD = gridProv.getNodeSubSeisMFD(i);
				if (nodeMFD != null) {
					for (int j=0; j<nodeMFD.size(); j++) {
						double myX = nodeMFD.getX(j);
						double myY = nodeMFD.getY(j);
						if (myX < u3MFD.getMinX() || myY == 0)
							continue;
						int index = u3MFD.getXIndex(myX);
						Preconditions.checkState(index >= 0, "Bad x mapping. myX=%s", myX);
						u3MFD.add(index, nodeMFD.getY(j));
					}
				}
			}
		}
		return u3MFD;
	}
	
	

}
