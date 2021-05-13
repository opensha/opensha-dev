package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.exceptions.ParameterException;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.ClassUtils;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.StdDevTypeParam;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceJBParameter;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceRupParameter;
import org.opensha.sha.imr.param.SiteParams.DepthTo1pt0kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;

import com.google.common.base.Preconditions;
import com.google.common.base.Splitter;
import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.griddedSeismicity.AbstractGridSourceProvider;
import scratch.UCERF3.griddedSeismicity.GridSourceProvider;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.kevin.ucerf3.etas.MPJ_ETAS_Simulator;

public class UCERF3_IM_Calculator {
	
	private FaultSystemSolution sol;
	private Site site;
	private ScalarIMR gmpe;
	private boolean doGridded;
	private File outputFile;
	private List<ETAS_EqkRupture> etasCatalog;
	
	public UCERF3_IM_Calculator(CommandLine cmd) throws IOException, DocumentException {
		File fssFile = new File(cmd.getOptionValue("sol-file"));
		FaultSystemSolution sol = FaultSystemIO.loadSol(fssFile);

		double lat = Double.parseDouble(cmd.getOptionValue("latitude"));
		double lon = Double.parseDouble(cmd.getOptionValue("longitude"));
		double vs30 = Double.parseDouble(cmd.getOptionValue("vs30"));
		
		Site site = new Site(new Location(lat, lon));
		Vs30_Param vs30Param = new Vs30_Param();
		vs30Param.setValue(vs30);
		site.addParameter(vs30Param);
		
		if (cmd.hasOption("z10")) {
			double z10 = Double.parseDouble(cmd.getOptionValue("z10"));
			DepthTo1pt0kmPerSecParam z10Param = new DepthTo1pt0kmPerSecParam(z10, true);
			z10Param.setValue(z10);
			site.addParameter(z10Param);
		}
		
		if (cmd.hasOption("z25")) {
			double z25 = Double.parseDouble(cmd.getOptionValue("z25"));
			DepthTo2pt5kmPerSecParam z25Param = new DepthTo2pt5kmPerSecParam(z25, true);
			z25Param.setValue(z25);
			site.addParameter(z25Param);
		}
		
		AttenRelRef gmpe = AttenRelRef.valueOf(cmd.getOptionValue("gmpe"));
		
		String imt = cmd.getOptionValue("imt");
		double period = Double.NaN;
		if (imt.equals(SA_Param.NAME)) {
			Preconditions.checkArgument(cmd.hasOption("period"));
			period = Double.parseDouble(cmd.getOptionValue("period"));
		}
		
		boolean doGridded = cmd.hasOption("do-gridded");
		
		File outputFile = new File(cmd.getOptionValue("output-file"));
		
		List<ETAS_EqkRupture> etasCatalog = null;
		if (cmd.hasOption("etas-catalog"))
			etasCatalog = ETAS_CatalogIO.loadCatalog(new File(cmd.getOptionValue("etas-catalog")));
		
		init(sol, site, gmpe, doGridded, outputFile, etasCatalog, imt, period);
	}
	
	private void init(FaultSystemSolution sol, Site site, AttenRelRef gmpe, boolean doGridded,
			File outputFile, List<ETAS_EqkRupture> etasCatalog, String imt, double period) {
		this.sol = sol;
		this.site = site;
		this.doGridded = doGridded;
		this.outputFile = outputFile;
		this.etasCatalog = etasCatalog;
		
		this.gmpe = gmpe.instance(null);
		this.gmpe.setParamDefaults();
		this.gmpe.setIntensityMeasure(imt);
		if (imt.equals(SA_Param.NAME))
			SA_Param.setPeriodInSA_Param(this.gmpe.getIntensityMeasure(), period);
		
		for (Parameter<?> param : this.gmpe.getSiteParams()) {
			if (!site.containsParameter(param))
				site.addParameter(param);
		}
	}
	
	private static final double ERF_DURATION = 1d;
	
	public void calculate() throws IOException {
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		if (doGridded)
			erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
		else
			erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		erf.getTimeSpan().setDuration(ERF_DURATION);
		erf.updateForecast();
		
		CSVFile<String> csv = new CSVFile<>(true);
		
		List<String> header = Lists.newArrayList("FSS Index", "Grid Node Index", "Mag", "Rake", "U3 Annual Rate",
				gmpe.getShortName()+" Log Mean", gmpe.getShortName()+" Total Std Dev");
		
		boolean hasInterStdDev;
		boolean hasIntraStdDev;
		try {
			StdDevTypeParam stdDevtypeParam = (StdDevTypeParam)gmpe.getParameter(StdDevTypeParam.NAME);
			hasInterStdDev = stdDevtypeParam.isAllowed(StdDevTypeParam.STD_DEV_TYPE_INTER);
			hasIntraStdDev = stdDevtypeParam.isAllowed(StdDevTypeParam.STD_DEV_TYPE_INTRA);
		} catch (ParameterException e) {
			hasInterStdDev = false;
			hasIntraStdDev = false;
		}
		if (hasInterStdDev)
			header.add("Iter-Event Std Dev");
		if (hasIntraStdDev)
			header.add("Itra-Event Std Dev");
		boolean hasRrup;
		boolean hasRJB;
		try {
			gmpe.getParameter(DistanceRupParameter.NAME);
			hasRrup = true;
		} catch (ParameterException e) {
			hasRrup = false;
		}
		try {
			gmpe.getParameter(DistanceJBParameter.NAME);
			hasRJB = true;
		} catch (ParameterException e) {
			hasRJB = false;
		}
		if (hasRrup)
			header.add("Distance Rup (km)");
		if (hasRJB)
			header.add("Distance J-B (km)");
		if (etasCatalog != null)
			header.add(0, "ETAS ID");
		csv.addLine(header);
		
		int numSourcesFSS = erf.getNumFaultSystemSources();
		
		if (etasCatalog == null) {
			// do the whole ERF
			for (int sourceID=0; sourceID<erf.getNumSources(); sourceID++) {
				ProbEqkSource source = erf.getSource(sourceID);
				int fssIndex = -1;
				int gridNodeIndex = -1;
				if (sourceID < numSourcesFSS)
					fssIndex = erf.getFltSysRupIndexForSource(sourceID);
				else
					gridNodeIndex = sourceID - numSourcesFSS;
				for (int rupID=0; rupID<source.getNumRuptures(); rupID++) {
					ProbEqkRupture rup = source.getRupture(rupID);
					csv.addLine(processRupture(rup, site, gmpe, null, fssIndex, gridNodeIndex, hasInterStdDev, hasIntraStdDev, hasRrup, hasRJB));
				}
			}
		} else {
			for (ETAS_EqkRupture etasRup : etasCatalog) {
				int fssIndex = etasRup.getFSSIndex();
				int gridNodeIndex = etasRup.getGridNodeIndex();
				Preconditions.checkState(fssIndex >= 0 || gridNodeIndex >= 0);
				if (fssIndex >= 0) {
					int sourceID = erf.getSrcIndexForFltSysRup(fssIndex);
					Preconditions.checkState(sourceID >= 0 && sourceID < numSourcesFSS);
					for (ProbEqkRupture rup : erf.getSource(sourceID))
						csv.addLine(processRupture(rup, site, gmpe, etasRup,
								fssIndex, gridNodeIndex, hasInterStdDev, hasIntraStdDev, hasRrup, hasRJB));
				} else if (doGridded) {
					if (etasRup.getMag() < AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF)
						continue;
					int sourceID = gridNodeIndex + numSourcesFSS;
					Preconditions.checkState(sourceID >= numSourcesFSS && sourceID < erf.getNumSources(),
							"Bad SourceID=%s for FSSIndex=%s, GridNodeIndex=%s, NumFSS=%s", sourceID, fssIndex, gridNodeIndex, numSourcesFSS);
					Location etasLoc = etasRup.getHypocenterLocation();
					for (ProbEqkRupture rup : erf.getSource(sourceID)) {
						Location rupLoc = rup.getRuptureSurface().getFirstLocOnUpperEdge();
						double distDegrees = Math.sqrt(Math.pow(rupLoc.getLatitude() - etasLoc.getLatitude(), 2)
								+ Math.pow(rupLoc.getLongitude() - etasLoc.getLongitude(), 2));
						Preconditions.checkState(distDegrees < 0.2);
						double magDiff = Math.abs(rup.getMag() - etasRup.getMag());
						if ((float)magDiff > 0.06)
							continue;
						csv.addLine(processRupture(rup, site, gmpe, etasRup,
								fssIndex, gridNodeIndex, hasInterStdDev, hasIntraStdDev, hasRrup, hasRJB));
					}
				}
			}
		}
		
		csv.writeToFile(outputFile);
	}
	
	private static double max_val = Double.NEGATIVE_INFINITY;
	
	private List<String> processRupture(ProbEqkRupture rup, Site site, ScalarIMR gmpe, ETAS_EqkRupture etasRup, int fssIndex, int gridNodeIndex,
			boolean hasInterStdDev, boolean hasIntraStdDev, boolean hasRrup, boolean hasRJB) {
		List<String> line = new ArrayList<>();
		if (etasRup != null)
			line.add(etasRup.getID()+"");
		line.add(fssIndex+"");
		line.add(gridNodeIndex+"");
		line.add((float)rup.getMag()+"");
		line.add((float)rup.getAveRake()+"");
		line.add((float)rup.getMeanAnnualRate(ERF_DURATION)+"");
		gmpe.setAll(rup, site, gmpe.getIntensityMeasure());
		line.add((float)gmpe.getMean()+"");
		max_val = Math.max(max_val, gmpe.getMean());
		line.add((float)gmpe.getStdDev()+"");
		if (hasInterStdDev) {
			gmpe.getParameter(StdDevTypeParam.NAME).setValue(StdDevTypeParam.STD_DEV_TYPE_INTER);
			line.add((float)gmpe.getStdDev()+"");
		}
		if (hasIntraStdDev) {
			gmpe.getParameter(StdDevTypeParam.NAME).setValue(StdDevTypeParam.STD_DEV_TYPE_INTRA);
			line.add((float)gmpe.getStdDev()+"");
		}
		if (hasInterStdDev || hasIntraStdDev)
			gmpe.getParameter(StdDevTypeParam.NAME).setValue(StdDevTypeParam.STD_DEV_TYPE_TOTAL);
		if (hasRrup)
			line.add(((Double)gmpe.getParameter(DistanceRupParameter.NAME).getValue()).floatValue()+"");
		if (hasRJB)
			line.add(((Double)gmpe.getParameter(DistanceJBParameter.NAME).getValue()).floatValue()+"");
		
		return line;
	}
	
	public static Options createOptions() {
		Options ops = new Options();
		
		Option solFile = new Option("s", "sol-file", true, "UCERF3 Fault System Solution file");
		solFile.setRequired(true);
		ops.addOption(solFile);
		
		Option lat = new Option("lat", "latitude", true, "Site latitude");
		lat.setRequired(true);
		ops.addOption(lat);
		
		Option lon = new Option("lon", "longitude", true, "Site longitude");
		lon.setRequired(true);
		ops.addOption(lon);
		
		Option vs30 = new Option("vs", "vs30", true, "Site Vs30 (m/s)");
		vs30.setRequired(true);
		ops.addOption(vs30);
		
		Option z10 = new Option("z10", "depth-1", true, "Site Z1.0 (m)");
		z10.setRequired(false);
		ops.addOption(z10);
		
		Option z25 = new Option("z25", "depth-25", true, "Site Z2.5 (km)");
		z25.setRequired(false);
		ops.addOption(z25);
		
		Option gmpe = new Option("g", "gmpe", true, "GMPE Name (e.g. ASK_2014)");
		gmpe.setRequired(true);
		ops.addOption(gmpe);
		
		Option imt = new Option("i", "imt", true, "GMPE Intensity Measure Type (e.g. PGA, SA, PGV)");
		imt.setRequired(true);
		ops.addOption(imt);
		
		Option period = new Option("p", "period", true, "Spectral period, required if IMT is SA");
		period.setRequired(false);
		ops.addOption(period);
		
		Option gridded = new Option("grid", "do-gridded", false, "Flag to enable gridded seismicity");
		gridded.setRequired(false);
		ops.addOption(gridded);
		
		Option outputFile = new Option("o", "output-file", true, "Path to output CSV file");
		outputFile.setRequired(true);
		ops.addOption(outputFile);
		
		Option etasCatalog = new Option("etas", "etas-catalog", true, "Optional path to ETAS catalog. If supplied, only ruptures in this catalog will be output");
		etasCatalog.setRequired(false);
		ops.addOption(etasCatalog);
		
		return ops;
	}
	
	public static void printHelp(Options options, String appName) {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp( appName, options, true );
		System.exit(2);
	}
	
	public static void printUsage(Options options, String appName) {
		HelpFormatter formatter = new HelpFormatter();
		PrintWriter pw = new PrintWriter(System.out);
		formatter.printUsage(pw, 80, appName, options);
		pw.flush();
		System.exit(2);
	}

	public static void main(String[] args) throws IOException, DocumentException {
		Options options = createOptions();
		
		if (args.length == 1 && args[0].equals("--hardcoded")) {
			String argStr = "--sol-file /home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
					+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip";
			argStr += " --latitude 34.038165 --longitude -118.266111";
			// Vs30 in m/s, Z1.0 in m, Z2.5 in km (I know, it's weird that they're different units)
			// Z1.0 and Z2.5 are optional
			argStr += " --vs30 225.6";
//			argStr += " --vs30 760 --z10 100 --z25 1.5";
			argStr += " --gmpe ASK_2014";
			argStr += " --imt PGA";
			argStr += " --period 0";
			// optional argument to specify an etas catalog
//			argStr += " --etas-catalog /path/to/catalog"
			// flag to enable gridded seismicity
			argStr += " --do-gridded";
			// output file (.csv)
			argStr += " --output-file /tmp/u3_ims.csv";
			args = Lists.newArrayList(Splitter.on(" ").split(argStr)).toArray(new String[0]);
		} else if (args.length == 0) {
			printUsage(options, ClassUtils.getClassNameWithoutPackage(UCERF3_IM_Calculator.class));
		}
		
		CommandLineParser parser = new GnuParser();

		CommandLine cmd = null;
		try {
			cmd = parser.parse(options, args);
		} catch (ParseException e) {
			e.printStackTrace();
			printUsage(options, ClassUtils.getClassNameWithoutPackage(UCERF3_IM_Calculator.class));
		}
		
		UCERF3_IM_Calculator calc = new UCERF3_IM_Calculator(cmd);
		calc.calculate();
//		System.out.println("Max value: "+max_val);
	}

}
