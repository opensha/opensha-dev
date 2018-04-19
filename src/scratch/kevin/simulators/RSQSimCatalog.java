package scratch.kevin.simulators;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.charset.Charset;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.dom4j.Document;
import org.dom4j.DocumentException;
import org.dom4j.Element;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.metadata.XMLSaveable;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.commons.util.FileUtils;
import org.opensha.commons.util.XMLUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.CatalogLengthLoadIden;
import org.opensha.sha.simulators.iden.EventIDsRupIden;
import org.opensha.sha.simulators.iden.LogicalAndRupIden;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;
import org.opensha.sha.simulators.iden.MagRangeRuptureIdentifier;
import org.opensha.sha.simulators.iden.RegionIden;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.iden.SkipYearsLoadIden;
import org.opensha.sha.simulators.parsers.RSQSimFileReader;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimStateTransitionFileReader;
import org.opensha.sha.simulators.srf.RSQSimTransValidIden;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Table;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.IDPairing;
import scratch.UCERF3.utils.UCERF3_DataUtils;
import scratch.kevin.simulators.plots.AbstractPlot;
import scratch.kevin.simulators.plots.MFDPlot;
import scratch.kevin.simulators.plots.MagAreaScalingPlot;
import scratch.kevin.simulators.plots.NormalizedFaultRecurrenceIntervalPlot;
import scratch.kevin.simulators.plots.RecurrenceIntervalPlot;
import scratch.kevin.simulators.plots.RuptureVelocityPlot;
import scratch.kevin.simulators.plots.SectionRecurrenceComparePlot;
import scratch.kevin.simulators.plots.SectionRecurrenceComparePlot.SectType;
import scratch.kevin.simulators.plots.StationarityPlot;
import scratch.kevin.util.MarkdownUtils;
import scratch.kevin.util.MarkdownUtils.TableBuilder;

public class RSQSimCatalog implements XMLSaveable {
	
	public enum Catalogs {
		BRUCE_2142("bruce/rundir2142", "Bruce 2142", "Bruce Shaw", cal(2017, 6, 16),
				"Old projection; slip weakening; stress loaded; no creep correction",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2194("bruce/rundir2194", "Bruce 2194", "Bruce Shaw", cal(2017, 7, 5),
				"Catalog with decent large event scaling and distribution of sizes while not using"
				+ " any of the enhanced frictional weakening terms.", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		JG_UCERF3_millionElement("JG_UCERF3_millionElement", "U3 1mil Element Test", "Jacqui Gilchrist", cal(2017, 9, 27),
				"Test 1 million element catalog on UCERF3 fault system, ~0.25 km^2 trianglar elements",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2194_LONG("rundir2194_long", "Bruce 2194 Long", "Bruce Shaw (extended by Jacqui Gilchrist)", cal(2017, 8, 31),
				"Catalog with decent large event scaling and distribution of sizes while not using"
				+ " any of the enhanced frictional weakening terms.", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2273("bruce/rundir2273", "Bruce 2273", "Bruce Shaw", cal(2017, 10, 13),
				"Stress loading, more refined geometry, does not contain projection fix (some location discrepancies "
				+ "are present relative to UCERF3 faults).", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2310("bruce/rundir2310", "Bruce 2310", "Bruce Shaw", cal(2017, 10, 16),
				"Backslip loading, more refined geometry, projection fix (but all faults surface breaking)",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2320("bruce/rundir2320", "Bruce 2320", "Bruce Shaw", cal(2017, 10, 17),
				"Backslip loading, less refined geometry, projection fix (but all\n" + 
				"faults surface breaking), same as rundir2310 but less resolved",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2336("bruce/rundir2336", "Bruce 2336", "Bruce Shaw", cal(2017, 10, 20),
				"Larger slip velocity (1.5 m/s), backslipFromStress loading",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2337("bruce/rundir2337", "Bruce 2337", "Bruce Shaw", cal(2017, 10, 20),
				"Larger slip velocity (2.0 m/s), backslipFromStress loading",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2326("bruce/rundir2326", "Bruce 2326", "Bruce Shaw", cal(2017, 10, 23),
				"reference_1: a=.001 b=.008  Veq=1.0  sigmaN=100. backslipFromStress",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2342("bruce/rundir2342", "Bruce 2342", "Bruce Shaw", cal(2017, 10, 23),
				"larger Veq=1.2          relative to reference_1",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2343("bruce/rundir2343", "Bruce 2343", "Bruce Shaw", cal(2017, 10, 23),
				"smaller Veq=0.8        relative to reference_1",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2349("bruce/rundir2349", "Bruce 2349", "Bruce Shaw", cal(2017, 10, 23),
				"smaller sigmaN=80   relative to reference_1",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		JG_2194_K2("rundir2194_K2", "JG 2194 K2", "Jacqui Gilchrist", cal(2017, 10, 16),
				"Keith's fault geometry, normal backslip with U3 geologic long-term slip rates,"
				+ " and the same parameter values as Bruce's 2194",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		JG_tuneBase1m("tuneBase1m", "JG Tune Base 1M", "Jacqui Gilchrist", cal(2017, 11, 2),
				"U3 fault geometry with 1km^2 triangles, normal backslip loading with U3 geologic slip rates,"
				+ "calibrated to U3 supraseismogenic recurrence intervals, and default a/b",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		JG_modLoad_testB("modLoad_testB", "JG Mod Load Test B", "Jacqui Gilchrist", cal(2017, 11, 14),
				"Bruce's modified loading with higher values of frictional parameters",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		JG_baseCatalog2("baseCatalog2", "JG Base Catalog 2", "Jacqui Gilchrist", cal(2017, 11, 16),
				"Untuned version of tuneBase1m. Same fault model and frictional parameters, without any stress adjustments",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		JG_tuneBaseSW_1e5("tuneBaseSW_1e5", "JG Tune Base SW", "Jacqui Gilchrist", cal(2017, 11, 20),
				"Tuned, additional slip weakening parameters using Keith's fault geometry. muSlipAmp = 0.2, muSlipInvDist_1 = 2.0, cohesion = 6.",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		JG_baseCatalogSW_10("baseCatalogSW_10", "JG Base SW", "Jacqui Gilchrist", cal(2017, 11, 20),
				"Untuned, additional slip weakening parameters using Keith's fault geometry. muSlipAmp = 0.2, muSlipInvDist_1 = 2.0, cohesion = 6.",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		JG_tunedBase1m_ddotEQmod("tunedBase1m_ddotEQmod", "JG Tune Base Mod Vel", "Jacqui Gilchrist", cal(2017, 11, 26),
				"New version of tuneBase1m, with patch-specific slip velocities.",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2381("bruce/rundir2381", "Bruce 2381", "Bruce Shaw", cal(2017, 12, 22),
				"fracCreep=0.5", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2388("bruce/rundir2388", "Bruce 2388", "Bruce Shaw", cal(2017, 12, 22),
				"fracCreep=0.25", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2410("bruce/rundir2410", "Bruce 2410", "Bruce Shaw", cal(2017, 12, 22),
				"factorNormal=2.0;  50ky reference time; 212ky spinup", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2411("bruce/rundir2411", "Bruce 2411", "Bruce Shaw", cal(2017, 12, 22),
				"factorNormal=2.0;  50ky reference time; 85ky spinup", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2412("bruce/rundir2412", "Bruce 2412", "Bruce Shaw", cal(2017, 12, 22),
				"factorNormal=2.0;  50ky reference time; 85ky spinup; srt(slipRate)", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2413("bruce/rundir2413", "Bruce 2413", "Bruce Shaw", cal(2017, 12, 22),
				"factorNormal=2.0;  100ky reference time; 85ky spinup; sqrt(slipRate)", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2457("bruce/rundir2457", "Bruce 2457", "Bruce Shaw", cal(2018, 1, 14),
				"new loading;  fCreep=0.25", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2495("bruce/rundir2495", "Bruce 2495", "Bruce Shaw", cal(2018, 1, 29),
				"flat loaded.  fracCreep=0.5.  maxDepth=14.4", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2579("bruce/rundir2579", "Bruce 2579", "Bruce Shaw", cal(2018, 2, 07),
				"straight loaded;  fracCreep=0.5;  H=18 (2,12,4);  stressMult=1.2;  neighbors",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2585("bruce/rundir2585", "Bruce 2585", "Bruce Shaw", cal(2018, 2, 10),
				"Longer run; else same as r2579. straight loaded;  fracCreep=0.5; H=18 (2,12,4); stressMult=1.2; neighbors",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2592("bruce/rundir2592", "Bruce 2592", "Bruce Shaw", cal(2018, 2, 11),
				"straight loaded;  fracCreep=0.5;  H=16 (2,11,3); stressMult=1.2; neighbors",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2620("bruce/rundir2620", "Bruce 2620", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   b=.006", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2621("bruce/rundir2621", "Bruce 2621", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   b=.007", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2622("bruce/rundir2622", "Bruce 2622", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   b=.009", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2623("bruce/rundir2623", "Bruce 2623", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   b=.010", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2624("bruce/rundir2624", "Bruce 2624", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   a=.0008", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2625("bruce/rundir2625", "Bruce 2625", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   a=.0009", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2626("bruce/rundir2626", "Bruce 2626", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   a=.0011", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2627("bruce/rundir2627", "Bruce 2627", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   a=.0012", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2628("bruce/rundir2628", "Bruce 2628", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   N=80", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2629("bruce/rundir2629", "Bruce 2629", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   N=90", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2630("bruce/rundir2630", "Bruce 2630", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   N=110", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2631("bruce/rundir2631", "Bruce 2631", "Bruce Shaw", cal(2018, 3, 8),
				"sensitivity test, diff r2585   N=120", FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2585_1MYR("rundir2585_1myr", "Bruce 2585 1myr", "Bruce Shaw/Jacqui Gilchrist", cal(2018, 3, 12),
				"Extended version of Bruce's 2585 to 1 million years",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2616("bruce/rundir2616", "Bruce 2616", "Bruce Shaw", cal(2018, 3, 19),
				"similar to r2585, but bigger seismogenic depth: H=18 (2,13,3)",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2636("bruce/rundir2636", "Bruce 2636", "Bruce Shaw", cal(2018, 3, 28),
				"sensitivity test, diff r2585 a=.0013",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2637("bruce/rundir2637", "Bruce 2637", "Bruce Shaw", cal(2018, 3, 28),
				"sensitivity test, diff r2585  N=130",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2665("bruce/rundir2665", "Bruce 2665", "Bruce Shaw", cal(2018, 4, 16),
				"dx/2, LatCut=37, rateCut=0.2mm/yr, interpolated nearest",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC),
		BRUCE_2666("bruce/rundir2666", "Bruce 2666", "Bruce Shaw", cal(2018, 4, 17),
				"dx/4, LatCut=37, rateCut=2.0mm/yr, interpolated nearest",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC);
		
		private String dirName;
		private RSQSimCatalog catalog;
		
		private Catalogs(String dirName, String name, String author, GregorianCalendar date, String metadata,
				FaultModels fm, DeformationModels dm) {
			this.dirName = dirName;
			catalog = new RSQSimCatalog(name, author, date, metadata, fm, dm);
		}
		
		public RSQSimCatalog instance(File baseDir) {
			File dir = new File(baseDir, dirName);
			Preconditions.checkState(dir.exists(), "Catalog dir doesn't exist: %s", dir.getAbsolutePath());
			catalog.dir = dir;
			return catalog;
		}
	}
	
	private File dir;
	private String name;
	private String author;
	private GregorianCalendar date;
	private String metadata;
	private FaultModels fm;
	private DeformationModels dm;
	
	private double constSlipVel = Double.NaN;
	private Map<Integer, Double> slipVels = null;
	private double aveArea = Double.NaN;
	private int numEvents = -1;
	private double durationYears = Double.NaN;
	private Map<String, String> params;
	
	private List<SimulatorElement> elements;
	private RSQSimStateTransitionFileReader transReader;
	private List<FaultSectionPrefData> subSects;
	private Map<Integer, Double> subSectAreas;
	private Map<IDPairing, Double> subSectDistsCache;
	
	private static Table<FaultModels, DeformationModels, FaultSystemSolution> compSolsTable = HashBasedTable.create();
	
	public static final String XML_METADATA_NAME = "RSQSimCatalog";
	
	private RSQSimCatalog(String name, String author, GregorianCalendar date, String metadata,
			FaultModels fm, DeformationModels dm) {
		this(null, name, author, date, metadata, fm, dm);
	}

	public RSQSimCatalog(File dir, String name, String author, GregorianCalendar date, String metadata,
			FaultModels fm, DeformationModels dm) {
		this.dir = dir;
		this.name = name;
		this.author = author;
		this.date = date;
		this.metadata = metadata;
		this.fm = fm;
		this.dm = dm;
	}
	
	public File getCatalogDir() {
		return dir;
	}

	public String getName() {
		return name;
	}

	public String getAuthor() {
		return author;
	}

	public GregorianCalendar getDate() {
		return date;
	}

	public String getMetadata() {
		return metadata;
	}
	
	public FaultModels getFaultModel() {
		return fm;
	}

	public DeformationModels getDeformationModel() {
		return dm;
	}
	
	public synchronized Map<Integer, Double> getSlipVelocities() throws IOException {
		if (slipVels == null) {
			if (Doubles.isFinite(constSlipVel) && constSlipVel > 0) {
				List<SimulatorElement> elems = getElements();
				slipVels = new HashMap<>();
				for (int i=0; i<elems.size(); i++)
					slipVels.put(elems.get(i).getID(), constSlipVel);
			} else {
				Map<String, String> params = getParams();
				String ddotEQFname = params.get("ddotEQFname");
				if (ddotEQFname != null && !ddotEQFname.trim().isEmpty()) {
					File ddotEQFile = new File(getCatalogDir(), ddotEQFname);
					Preconditions.checkState(ddotEQFile.exists(),
							"ddotEQFname = %s doesn't exist in %s", ddotEQFname, getCatalogDir().getAbsolutePath());
					double[] velArray = loadDoubleInputFile(ddotEQFile);
					List<SimulatorElement> elems = getElements();
					Preconditions.checkState(velArray.length == elems.size(), "expected %s patch velocities, have %s", elems.size(), velArray.length);
					slipVels = new HashMap<>();
					for (int i=0; i<elems.size(); i++)
						slipVels.put(elems.get(i).getID(), velArray[i]);
					constSlipVel = Double.NaN;
				} else {
					String ddotEQ = params.get("ddotEQ_1");
					Preconditions.checkNotNull(ddotEQ, "ddotEQ_1 not in params file");
					constSlipVel = Double.parseDouble(ddotEQ);
					List<SimulatorElement> elems = getElements();
					slipVels = new HashMap<>();
					for (int i=0; i<elems.size(); i++)
						slipVels.put(elems.get(i).getID(), constSlipVel);
				}
			}
		}
		return slipVels;
	}
	
	private static double[] loadDoubleInputFile(File file) throws IOException {
		List<Double> vals = new ArrayList<>();
		for (String line : Files.readLines(file, Charset.defaultCharset())) {
			line = line.trim();
			if (line.isEmpty() || line.startsWith("#"))
				continue;
			vals.add(Double.parseDouble(line));
		}
		return Doubles.toArray(vals);
	}
	
	private static DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd");
	private static DecimalFormat areaDF = new DecimalFormat("0.00");
	private static DecimalFormat groupedIntDF = new DecimalFormat("#");
	static {
		groupedIntDF.setGroupingUsed(true);
		groupedIntDF.setGroupingSize(3);
	}
	
	public List<String> getMarkdownMetadataTable() {
		TableBuilder builder = MarkdownUtils.tableBuilder();
		builder.addLine("**Catalog**", getName());
		builder.addLine("**Author**", getAuthor()+", "+dateFormat.format(getDate().getTime()));
		builder.addLine("**Description**", getMetadata());
		builder.addLine("**Fault/Def Model**", fm+", "+dm);
		try {
			Map<Integer, Double> slipVels = getSlipVelocities();
			String velStr;
			if (Double.isFinite(constSlipVel)) {
				velStr = (float)constSlipVel+" m/s";
			} else {
				MinMaxAveTracker velTrack = new MinMaxAveTracker();
				for (Double vel : slipVels.values())
					velTrack.addValue(vel);
				velStr = "Variable, range=["+(float)velTrack.getMin()+" "+(float)velTrack.getMax()+"], mean="+(float)velTrack.getAverage();
			}
			builder.addLine("**Slip Velocity**", velStr);
		} catch (IOException e) {
			e.printStackTrace();
		}
		try {
			double  aveArea = getAveArea();
			builder.addLine("**Average Element Area**", areaDF.format(aveArea)+" km^2");
		} catch (IOException e) {
			e.printStackTrace();
		}
		try {
			int numEvents = getNumEvents();
			double durationYears = getDurationYears();
			builder.addLine("**Length**", groupedIntDF.format(numEvents)+" events in "
					+groupedIntDF.format(durationYears)+" years");
		} catch (IOException e) {
			e.printStackTrace();
		}
		try {
			Map<String, String> params = getParams();
			if (params != null) {
				double a = Double.parseDouble(params.get("A_1"));
				double b = Double.parseDouble(params.get("B_1"));
				String ddotEQ = params.get("ddotEQ_1");
				if (ddotEQ.contains("."))
					ddotEQ = Float.parseFloat(ddotEQ)+"";
				builder.addLine("**Frictional Params**", "a="+(float)a+", b="+(float)b+", (b-a)="+(float)(b-a)+", ddotEQ="+ddotEQ);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return builder.build();
	}
	
	public synchronized double getAveArea() throws IOException {
		if (Double.isNaN(aveArea)) {
			List<SimulatorElement> elements = getElements();
			aveArea = 0d;
			for (SimulatorElement e : elements)
				aveArea += e.getArea()*1e-6;
			aveArea /= elements.size();
		}
		return aveArea;
	}
	
	public synchronized int getNumEvents() throws IOException {
		if (numEvents < 0)
			numEvents = RSQSimFileReader.getNumEvents(getCatalogDir());
		return numEvents;
	}
		
	public synchronized double getDurationYears() throws IOException {
		if (Double.isNaN(durationYears))
			durationYears = RSQSimFileReader.getDurationYears(getCatalogDir());
		return durationYears;
	}
	
	private synchronized Map<String, String> getParams() throws IOException {
		if (params == null) {
			File paramFile = findParamFile();
			if (paramFile != null) {
				System.out.println("Loading params from "+paramFile.getAbsolutePath());
				params = new HashMap<>();
				for (String line : Files.readLines(paramFile, Charset.defaultCharset())) {
					line = line.trim();
					if (line.contains("=")) {
						int ind = line.indexOf("=");
						String key = line.substring(0, ind).trim();
						String val = line.substring(ind+1).trim();
						params.put(key, val);
					}
				}
			}
		}
		return params;
	}
	
	private File findParamFile() throws IOException {
		File dir = getCatalogDir();
		File bruceInFile = new File(dir, "multiparam.in");
		if (bruceInFile.exists())
			return bruceInFile;
		for (File file : getCatalogDir().listFiles()) {
			String name = file.getName();
			if (!name.endsWith(".in"))
				continue;
			String lower = name.toLowerCase();
//			if (lower.contains("deepen") || lower.contains("dotmod"))
//				continue;
			if (isParamFile(file))
				return file;
		}
		return null;
	}
	
	private boolean isParamFile(File file) throws IOException {
		int max = 1000;
		int count = 0;
		for (String line : Files.readLines(file, Charset.defaultCharset())) {
			line = line.trim();
			if (line.startsWith("A_1"))
				return true;
			if (count++ > max)
				return false;
		}
		return false;
	}
	
	public void writeMarkdownSummary(File dir, boolean plots, boolean replot) throws IOException {
		List<String> lines = new LinkedList<>();
		String topLink = "*[(top)](#"+MarkdownUtils.getAnchorName(getName())+")*";
		lines.add("# "+getName());
		lines.add("## Metadata");
		lines.addAll(getMarkdownMetadataTable());
		lines.add("");
		int tocIndex = lines.size();
		
		List<String> eventLinks = new ArrayList<>();
		List<String> eventNames = new ArrayList<>();
		
		List<String> gmpeLinks = new ArrayList<>();
		List<String> gmpeNames = new ArrayList<>();
		
		List<String> gmpeGriddedLinks = new ArrayList<>();
		List<String> gmpeGriddedNames = new ArrayList<>();
		
		List<String> gmpeRGLinks = new ArrayList<>();
		List<String> gmpeRGNames = new ArrayList<>();
		
		List<String> hazardLinks = new ArrayList<>();
		List<String> hazardNames = new ArrayList<>();
		
		List<String> hazardClusterLinks = new ArrayList<>();
		List<String> hazardClusterNames = new ArrayList<>();
		
		String rotDDLink = null;
		String multiFaultLink = null;
		String extremeEventLink = null;
		
		File[] dirList = dir.listFiles();
		Arrays.sort(dirList, new FileNameComparator());
		for (File subDir : dirList) {
			if (!subDir.isDirectory())
				continue;
			File mdFile = new File(subDir, "README.md");
			if (!mdFile.exists())
				continue;
			String name = subDir.getName();
			if (name.startsWith("event_")) {
				eventNames.add(MarkdownUtils.getTitle(mdFile));
				eventLinks.add(name);
			} else if (name.startsWith("gmpe_bbp_comparisons_")) {
				String subName = name.substring("gmpe_bbp_comparisons_".length());
				if (name.endsWith("_GriddedSites")) {
					subName = subName.substring(0, subName.indexOf("_GriddedSites"));
					gmpeGriddedNames.add(subName);
					gmpeGriddedLinks.add(name);
				} else {
					gmpeNames.add(subName);
					gmpeLinks.add(name);
				}
			} else if (name.startsWith("gmpe_bbp_rg_comparisons_")) {
				gmpeRGNames.add(name.substring("gmpe_bbp_rg_comparisons_".length()));
				gmpeRGLinks.add(name);
			} else if (name.equals("catalog_rotd_ratio_comparisons")) {
				Preconditions.checkState(rotDDLink == null, "Duplicate RotDD dirs! %s and %s", name, rotDDLink);
				rotDDLink = name;
			} else if (name.startsWith("hazard_")) {
				String hazName;
				if (name.contains("_pga")) {
					hazName = "PGA";
				} else if (name.contains("t_dependence")) {
					hazName = "T Dependence";
				} else {
					Preconditions.checkState(name.contains("_sa"));
					String periodStr = name.substring(name.indexOf("_sa_")+4);
					periodStr = periodStr.substring(0, periodStr.indexOf("s"));
					double period = Double.parseDouble(periodStr);
					hazName = (float)period+"s SA";
				}
				
				if (name.contains("_sigma")) {
					String sigmaStr = name.substring(name.indexOf("_sigma")+6);
					if (sigmaStr.contains("_"))
						sigmaStr = sigmaStr.substring(0, sigmaStr.indexOf("_"));
					hazName += ", Fixed σ="+sigmaStr;
				}
				
				if (name.contains("_gmpe")) {
					String gmpeStr = name.substring(name.indexOf("_gmpe")+5);
					if (gmpeStr.contains("_"))
						gmpeStr = gmpeStr.substring(0, gmpeStr.indexOf("_"));
					hazName += ", "+gmpeStr;
				}
				
				if (name.contains("_sectArea")) {
					String areaStr = name.substring(name.indexOf("_sectArea")+9);
					if (areaStr.contains("_"))
						areaStr = areaStr.substring(0, areaStr.indexOf("_"));
					hazName += ", SectAreaFract="+areaStr;
				}
				
				if (name.contains("_cluster")) {
					hazardClusterLinks.add(name);
					hazardClusterNames.add(hazName);
				} else {
					hazardLinks.add(name);
					hazardNames.add(hazName);
				}
			} else if (name.equals("multi_fault")) {
				Preconditions.checkState(multiFaultLink == null, "Duplicate Multi Fault dirs! %s and %s", name, multiFaultLink);
				multiFaultLink = name;
			} else if (name.equals("extreme_events")) {
				Preconditions.checkState(extremeEventLink == null, "Duplicate Extreme Event dirs! %s and %s", name, multiFaultLink);
				extremeEventLink = name;
			}
		}
		
		if (!eventNames.isEmpty()) {
			lines.add("");
			lines.add("## Single Event Comparisons");
			lines.add(topLink);
			lines.add("");
			for (int i=0; i<eventNames.size(); i++)
				lines.add("* ["+eventNames.get(i)+"]("+eventLinks.get(i)+"/)");
		}
		if (!gmpeNames.isEmpty() || !gmpeGriddedNames.isEmpty()) {
			lines.add("");
			lines.add("## Full Catalog GMPE Comparisons");
			lines.add(topLink);
			lines.add("");
//			System.out.print("Have "+gmpeNames.size()+" regulars and "+gmpeGriddedNames.size()+" gridded");
			boolean both = !gmpeNames.isEmpty() && !gmpeGriddedNames.isEmpty();
			if (both) {
				lines.add("### Points Of Interest");
				lines.add("");
				for (int i=0; i<gmpeNames.size(); i++)
					lines.add("* ["+gmpeNames.get(i)+"]("+gmpeLinks.get(i)+"/)");
				lines.add("");
				lines.add("### Gridded Sites");
				lines.add("");
				for (int i=0; i<gmpeGriddedNames.size(); i++)
					lines.add("* ["+gmpeGriddedNames.get(i)+"]("+gmpeGriddedLinks.get(i)+"/)");
			} else {
				for (int i=0; i<gmpeNames.size(); i++)
					lines.add("* ["+gmpeNames.get(i)+"]("+gmpeLinks.get(i)+"/)");
				for (int i=0; i<gmpeGriddedNames.size(); i++)
					lines.add("* ["+gmpeGriddedNames.get(i)+"]("+gmpeGriddedLinks.get(i)+"/)");
			}
		}
		if (!gmpeRGNames.isEmpty()) {
			lines.add("");
			lines.add("## Full Catalog GMPE Comparisons with BBP Rupture Generator");
			lines.add(topLink);
			lines.add("");
			for (int i=0; i<gmpeRGNames.size(); i++)
				lines.add("* ["+gmpeRGNames.get(i)+"]("+gmpeRGLinks.get(i)+"/)");
		}
		if (rotDDLink != null) {
			lines.add("");
			lines.add("## Full Catalog RotD100/RotD50 Ratios");
			lines.add(topLink);
			lines.add("");
			lines.add("[Full Catalog RotD100/RotD50 Ratios Plotted Here]("+rotDDLink+"/)");
		}
		if (!hazardLinks.isEmpty()) {
			lines.add("");
			lines.add("## Hazard Comparisons");
			lines.add(topLink);
			lines.add("");
			for (int i=0; i<hazardNames.size(); i++)
				lines.add("* ["+hazardNames.get(i)+"]("+hazardLinks.get(i)+"/)");
		}
		if (!hazardClusterLinks.isEmpty()) {
			lines.add("");
			lines.add("## Hazard Clustering Comparisons");
			lines.add(topLink);
			lines.add("");
			for (int i=0; i<hazardClusterNames.size(); i++)
				lines.add("* ["+hazardClusterNames.get(i)+"]("+hazardClusterLinks.get(i)+"/)");
		}
		if (multiFaultLink != null) {
			lines.add("");
			lines.add("## Multi-Fault Rupture Comparisons");
			lines.add(topLink);
			lines.add("");
			lines.add("[Multi-Fault Rupture Comparisons here]("+multiFaultLink+"/)");
		}
		if (extremeEventLink != null) {
			lines.add("");
			lines.add("## Extreme Event Examples");
			lines.add(topLink);
			lines.add("");
			lines.add("[Extreme Event Examples Here]("+extremeEventLink+"/)");
		}
		
		if (plots) {
			File resourcesDir = new File(dir, "resources");
			Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
			lines.add("");
			int skipYears;
			double duration = getDurationYears();
			if (duration > 20000)
				skipYears = 5000;
			else if (duration > 10000)
				skipYears = 3000;
			else if (duration > 1000)
				skipYears = 1000;
			else
				skipYears = 0;
			lines.addAll(writeStandardDiagnosticPlots(resourcesDir, skipYears, 6d, replot, topLink));
		}
		
		File inputFile = findParamFile();
		if (params != null) {
			lines.add("");
			lines.add("## Input File");
			lines.add(topLink);
			lines.add("");
			lines.add("```");
			for (String line : Files.readLines(inputFile, Charset.defaultCharset()))
				lines.add(line);
			lines.add("```");
		}
		
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		
		MarkdownUtils.writeReadmeAndHTML(lines, dir);
		
		// write metadata
		XMLUtils.writeObjectToXMLAsRoot(this, new File(dir, "catalog.xml"));
	}

	public File getGeomFile() throws FileNotFoundException {
		for (File file : dir.listFiles()) {
			String name = file.getName().toLowerCase();
			if (name.endsWith(".flt"))
				return file;
			if (name.startsWith("zfault") && name.endsWith(".in") && !name.contains("deepen_"))
				return file;
		}
		throw new FileNotFoundException("No geometry file found in "+dir.getAbsolutePath());
	}
	
	public synchronized List<SimulatorElement> getElements() throws IOException {
		if (elements == null) {
			File geomFile = getGeomFile();
			elements = RSQSimFileReader.readGeometryFile(geomFile, 11, 'N');
		}
		return elements;
	}
	
	public Loader loader() throws IOException {
		return new Loader(getElements(), getCatalogDir());
	}
	
	private static File getTransFile(File dir) throws FileNotFoundException {
		for (File file : dir.listFiles()) {
			String name = file.getName().toLowerCase();
			if (name.startsWith("trans.") && name.endsWith(".out"))
				return file;
		}
		throw new FileNotFoundException("No transitions file found in "+dir.getAbsolutePath());
	}
	
	public synchronized RSQSimStateTransitionFileReader getTransitions() throws IOException {
		if (transReader == null) {
			File transFile = getTransFile(getCatalogDir());
			transReader = new RSQSimStateTransitionFileReader(transFile, getElements());
		}
		return transReader;
	}
	
	public synchronized RSQSimEventSlipTimeFunc getSlipTimeFunc(RSQSimEvent event) throws IOException {
		return new RSQSimEventSlipTimeFunc(getTransitions().getTransitions(event), getSlipVelocities());
	}

	private static GregorianCalendar cal(int year, int month, int day) {
		return new GregorianCalendar(year, month-1, day);
	}
	
	public synchronized List<FaultSectionPrefData> getU3SubSects() {
		if (subSects == null)
			subSects = RSQSimUtils.getUCERF3SubSectsForComparison(getFaultModel(), getDeformationModel());
		return subSects;
	}
	
	private File getSolCacheDir() {
		File scratchDir = UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR;
		if (scratchDir.exists()) {
			// eclipse project
			File dir = new File(scratchDir, "ucerf3_fm_dm_sols");
			if (!dir.exists())
				Preconditions.checkState(dir.mkdir());
			return dir;
		} else {
			// use home dir
			String path = System.getProperty("user.home");
			File homeDir = new File(path);
			Preconditions.checkState(homeDir.exists(), "user.home dir doesn't exist: "+path);
			File openSHADir = new File(homeDir, ".opensha");
			if (!openSHADir.exists())
				Preconditions.checkState(openSHADir.mkdir(),
						"Couldn't create OpenSHA store location: "+openSHADir.getAbsolutePath());
			File uc3Dir = new File(openSHADir, "ucerf3_fm_dm_sols");
			if (!uc3Dir.exists())
				Preconditions.checkState(uc3Dir.mkdir(),
						"Couldn't create UCERF3 ERF store location: "+uc3Dir.getAbsolutePath());
			return uc3Dir;
		}
	}
	
	public FaultSystemSolution getU3CompareSol() throws IOException {
		synchronized (compSolsTable) {
			FaultSystemSolution sol = compSolsTable.get(fm, dm);
			
			if (sol == null) {
				File solDir = getSolCacheDir();
				File solFile = new File(solDir, fm.encodeChoiceString()+"_"+dm.encodeChoiceString()+"_MEAN_BRANCH_AVG_SOL.zip");
				if (!solFile.exists()) {
					// download it
					String addr = "http://opensha.usc.edu/ftp/kmilner/ucerf3/2013_05_10-ucerf3p3-production-10runs_fm_dm_sub_plots/"
							+ fm.encodeChoiceString()+"_"+dm.encodeChoiceString()+"/"+solFile.getName();
					FileUtils.downloadURL(addr, solFile);
				}
				
				try {
					sol = FaultSystemIO.loadSol(solFile);
				} catch (DocumentException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
				compSolsTable.put(fm, dm, sol);
			}
			
			return sol;
		}
	}
	
	private synchronized Map<Integer, Double> getSubSectAreas() throws IOException {
		if (subSectAreas == null)
			subSectAreas = RSQSimUtils.calcSubSectAreas(getElements(), getU3SubSects());
		return subSectAreas;
	}
	
	public synchronized Map<IDPairing, Double> getSubSectDistsCache() {
		if (subSectDistsCache == null)
			subSectDistsCache = new HashMap<>();
		return subSectDistsCache;
	}
	
	public EqkRupture getGMPE_Rupture(RSQSimEvent event, double minFractForInclusion) {
		List<SimulatorElement> elements;
		Map<Integer, Double> subSectAreas = null;
		try {
			elements = getElements();
			if (minFractForInclusion > 0d)
				subSectAreas = getSubSectAreas();
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		return RSQSimUtils.buildSubSectBasedRupture(event, getU3SubSects(), elements,
				minFractForInclusion, subSectAreas, getSubSectDistsCache());
	}
	
	public List<FaultSectionPrefData> getSubSectsForRupture(RSQSimEvent event, double minFractForInclusion) {
		List<SimulatorElement> elements;
		Map<Integer, Double> subSectAreas = null;
		try {
			elements = getElements();
			if (minFractForInclusion > 0d)
				subSectAreas = getSubSectAreas();
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		int minSectIndex = RSQSimUtils.getSubSectIndexOffset(elements, getU3SubSects());
		
		List<List<FaultSectionPrefData>> bundled =  RSQSimUtils.getSectionsForRupture(event, minSectIndex,
				getU3SubSects(), getSubSectDistsCache(), minFractForInclusion, subSectAreas);
		List<FaultSectionPrefData> allSects = new ArrayList<>();
		for (List<FaultSectionPrefData> sects : bundled)
			allSects.addAll(sects);
		return allSects;
	}
	
	public class Loader {
		private List<SimulatorElement> elements;
		private File catalogDir;
		
		private List<RuptureIdentifier> loadIdens;
		
		private Loader(List<SimulatorElement> elements, File catalogDir) {
			super();
			this.elements = elements;
			this.catalogDir = catalogDir;
			
			loadIdens = new ArrayList<>();
		}
		
		public Loader minMag(double minMag) {
			loadIdens.add(new MagRangeRuptureIdentifier(minMag, Double.POSITIVE_INFINITY));
			return this;
		}
		
		public Loader maxMag(double maxMag) {
			loadIdens.add(new MagRangeRuptureIdentifier(Double.NEGATIVE_INFINITY, maxMag));
			return this;
		}
		
		public Loader skipYears(double years) {
			if (years > 0)
				loadIdens.add(new SkipYearsLoadIden(years));
			return this;
		}
		
		public Loader maxDuration(double years) {
			loadIdens.add(new CatalogLengthLoadIden(years));
			return this;
		}
		
		public Loader matches(RuptureIdentifier iden) {
			loadIdens.add(iden);
			return this;
		}
		
		public Loader withinCutoffDist(double maxDist, Collection<Location> locs) {
			return withinCutoffDist(maxDist, locs.toArray(new Location[0]));
		}
		
		public Loader withinCutoffDist(double maxDist, Location... locs) {
			List<RegionIden> regIdens = new ArrayList<>();
			for (Location loc : locs)
				regIdens.add(new RegionIden(new Region(loc, maxDist)));
			return matches(new LogicalOrRupIden(regIdens));
		}
		
		public Loader hasTransitions() throws IOException {
			loadIdens.add(new RSQSimTransValidIden(getTransitions(), getSlipVelocities()));
			return this;
		}
		
		public RSQSimEvent byID(int eventID) throws IOException {
			List<RSQSimEvent> events = this.byIDs(eventID);
			Preconditions.checkState(events.size() == 1, "Event "+eventID+" not found");
			return events.get(0);
		}
		
		public List<RSQSimEvent> byIDs(int... eventIDs) throws IOException {
			loadIdens.add(new EventIDsRupIden(eventIDs));
			return this.load();
		}
		
		public List<RSQSimEvent> load() throws IOException {
			LogicalAndRupIden loadIden = new LogicalAndRupIden(loadIdens);
			List<RuptureIdentifier> rupIdens = new ArrayList<>();
			rupIdens.add(loadIden);
			return RSQSimFileReader.readEventsFile(catalogDir, elements, rupIdens);
		}
		
		public Iterable<RSQSimEvent> iterable() throws IOException {
			LogicalAndRupIden loadIden = new LogicalAndRupIden(loadIdens);
			List<RuptureIdentifier> rupIdens = new ArrayList<>();
			rupIdens.add(loadIden);
			return RSQSimFileReader.getEventsIterable(catalogDir, elements, rupIdens);
		}
	}
	
	public FaultSystemSolution buildSolution(Loader loader, double minMag, double minFractForInclusion) throws IOException {
		return buildSolution(loader.load(), minMag, minFractForInclusion);
	}
	
	public FaultSystemSolution buildSolution(List<RSQSimEvent> events, double minMag, double minFractForInclusion) throws IOException {
		return RSQSimUtils.buildFaultSystemSolution(getU3SubSects(), getElements(), events, minMag, minFractForInclusion);
	}
	
	public List<String> writeStandardDiagnosticPlots(File outputDir, int skipYears, double minMag, boolean replot, String topLink)
			throws IOException {
		List<String> lines = new ArrayList<>();
		lines.add("## Plots");
		List<AbstractPlot> plots = new ArrayList<>();
		
		TableBuilder table;
		
		if (replot || !new File(outputDir, "mfd.png").exists()) {
			MFDPlot mfdPlot = new MFDPlot(minMag);
			mfdPlot.initialize(getName(), outputDir, "mfd");
			plots.add(mfdPlot);
		}
		lines.add("### Magnitude-Frequency Plot");
		lines.add(topLink);
		lines.add("");
		lines.add("![MFD]("+outputDir.getName()+"/mfd.png)");
		
		if (replot || !new File(outputDir, "mag_area_hist2D.png").exists()) {
			MagAreaScalingPlot magAreaPlot = new MagAreaScalingPlot();
			magAreaPlot.initialize(getName(), outputDir, "mag_area");
			plots.add(magAreaPlot);
		}
		lines.add("### Magnitude-Area Plots");
		lines.add(topLink);
		lines.add("");
		table = MarkdownUtils.tableBuilder();
		table.addLine("Scatter", "2-D Hist");
		table.initNewLine();
		table.addColumn("![MA Scatter]("+outputDir.getName()+"/mag_area.png)");
		table.addColumn("![MA Hist]("+outputDir.getName()+"/mag_area_hist2D.png)");
		table.finalizeLine();
		lines.addAll(table.build());
		
		if (replot || !new File(outputDir, "rupture_velocity_scatter.png").exists()) {
			RuptureVelocityPlot rupVelPlot = new RuptureVelocityPlot(getElements(), minMag);
			rupVelPlot.initialize(getName(), outputDir, "rupture_velocity");
			plots.add(rupVelPlot);			
		}
		lines.add("### Rupture Velocity Plots");
		lines.add(topLink);
		lines.add("");
		table = MarkdownUtils.tableBuilder();
		table.initNewLine().addColumn("**Scatter**");
		table.addColumn("![Rupture Velocity Scatter]("+outputDir.getName()+"/rupture_velocity_scatter.png)");
		table.finalizeLine().initNewLine().addColumn("**Distance/Velocity**");
		table.addColumn("![Rupture Velocity vs Dist]("+outputDir.getName()+"/rupture_velocity_vs_dist.png)");
		table.finalizeLine();
		lines.addAll(table.build());
		
		double[] riMinMags = {6d, 6.5, 7d, 7.5};
		while (minMag > riMinMags[0])
			riMinMags = Arrays.copyOfRange(riMinMags, 1, riMinMags.length);
		if (replot || !new File(outputDir, "interevent_times_m7.5.png").exists()) {
			RecurrenceIntervalPlot riPlot = new RecurrenceIntervalPlot(riMinMags);
			riPlot.initialize(getName(), outputDir, "interevent_times");
			plots.add(riPlot);
		}
		lines.add("### Global Interevent-Time Distributions");
		lines.add(topLink);
		lines.add("");
		table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		for (double riMinMag : riMinMags)
			if (riMinMag == Math.round(riMinMag))
				table.addColumn("**M≥"+(int)riMinMag+"**");
			else
				table.addColumn("**M≥"+(float)riMinMag+"**");
		table.finalizeLine().initNewLine();
		for (double riMinMag : riMinMags)
			if (riMinMag == Math.round(riMinMag))
				table.addColumn("![Interevent Times]("+outputDir.getName()+"/interevent_times_m"+(int)riMinMag+".png)");
			else
				table.addColumn("![Interevent Times]("+outputDir.getName()+"/interevent_times_m"+(float)riMinMag+".png)");
		table.finalizeLine();
		lines.addAll(table.build());
		
		double minFractForInclusion = 0.2;
		if (replot || !new File(outputDir, "norm_ri_elem_m7.5.png").exists()) {
			List<NormalizedFaultRecurrenceIntervalPlot> myPlots = new ArrayList<>();
			myPlots.add(new NormalizedFaultRecurrenceIntervalPlot(getElements(), riMinMags));
			myPlots.add(new NormalizedFaultRecurrenceIntervalPlot(getElements(), SectType.SUBSECTION,
					getU3SubSects(), minFractForInclusion, riMinMags));
			myPlots.add(new NormalizedFaultRecurrenceIntervalPlot(getElements(), SectType.PARENT,
					getU3SubSects(), minFractForInclusion, riMinMags));
			for (NormalizedFaultRecurrenceIntervalPlot plot : myPlots)
				plot.initialize(getName(), outputDir, "norm_ri_"+plot.getSectType().getPrefix());
			plots.addAll(myPlots);
		}
		lines.add("### Normalized Fault Interevent-Time Distributions");
		lines.add(topLink);
		lines.add("");
		table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		table.addColumn("");
		for (double riMinMag : riMinMags)
			if (riMinMag == Math.round(riMinMag))
				table.addColumn("**M≥"+(int)riMinMag+"**");
			else
				table.addColumn("**M≥"+(float)riMinMag+"**");
		table.finalizeLine();
		for (SectType type : new SectType[] {SectType.ELEMENT, SectType.SUBSECTION, SectType.PARENT }) {
			table.initNewLine();
			table.addColumn("**"+type.getSimType()+"s**");
			for (double riMinMag : riMinMags)
				if (riMinMag == Math.round(riMinMag))
					table.addColumn("![Norm RIs]("+outputDir.getName()+"/norm_ri_"+type.getPrefix()+"_m"+(int)riMinMag+".png)");
				else
					table.addColumn("![Norm RIs]("+outputDir.getName()+"/norm_ri_"+type.getPrefix()+"_m"+(float)riMinMag+".png)");
			table.finalizeLine();
		}
		lines.addAll(table.build());
		
		if (replot || !new File(outputDir, "stationarity.png").exists()) {
			StationarityPlot stationarityPlot = new StationarityPlot(minMag, 7d);
			stationarityPlot.initialize(getName(), outputDir, "stationarity");
			plots.add(stationarityPlot);
		}
		lines.add("### Stationarity Plot");
		lines.add(topLink);
		lines.add("");
		lines.add("![Stationarity]("+outputDir.getName()+"/stationarity.png)");
		
		String testMagStr;
		if (riMinMags[0] == Math.floor(riMinMags[0]))
			testMagStr = (int)riMinMags[0]+"";
		else
			testMagStr = (float)riMinMags[0]+"";
		if (replot || !new File(outputDir, "interevent_elements_m"+testMagStr+"_scatter.png").exists()) {
			SectionRecurrenceComparePlot elemCompare = new SectionRecurrenceComparePlot(getElements(), getU3CompareSol(), "UCERF3",
					SectionRecurrenceComparePlot.SectType.ELEMENT, 0, riMinMags);
			elemCompare.initialize(getName(), outputDir, "interevent_elements");
			plots.add(elemCompare);
			
			SectionRecurrenceComparePlot subSectCompare = new SectionRecurrenceComparePlot(getElements(), getU3CompareSol(), "UCERF3",
					SectionRecurrenceComparePlot.SectType.SUBSECTION, minFractForInclusion, riMinMags);
			subSectCompare.initialize(getName(), outputDir, "interevent_sub_sects");
			plots.add(subSectCompare);
		}
		lines.add("### Element/Subsection Interevent Time Comparisons");
		lines.add("");
		lines.add("#### Element Interevent Time Comparisons");
		lines.add(topLink);
		lines.add("");
		table = MarkdownUtils.tableBuilder();
		table.addLine("Min Mag", "Scatter", "2-D Hist");
		for (double riMinMag : riMinMags) {
			String prefix = "interevent_elements_m";
			if (riMinMag == Math.floor(riMinMag))
				prefix += (int)riMinMag;
			else
				prefix += (float)riMinMag;
			table.initNewLine();
			table.addColumn("**M≥"+(float)riMinMag+"**");
			table.addColumn("![Element Scatter]("+outputDir.getName()+"/"+prefix+"_scatter.png)");
			table.addColumn("![Element 2-D Hist]("+outputDir.getName()+"/"+prefix+"_hist2D.png)");
			table.finalizeLine();
		}
		lines.addAll(table.build());
		lines.add("");
		lines.add("#### Subsection Interevent Time Comparisons");
		lines.add(topLink);
		lines.add("");
		if (minFractForInclusion > 0) {
			lines.add("*Subsections participate in a rupture if at least "+(float)(minFractForInclusion*100d)+" % of its area ruptures*");
			lines.add("");
		}
		table = MarkdownUtils.tableBuilder();
		table.addLine("Min Mag", "Scatter", "2-D Hist");
		for (double riMinMag : riMinMags) {
			String prefix = "interevent_sub_sects_m";
			if (riMinMag == Math.floor(riMinMag))
				prefix += (int)riMinMag;
			else
				prefix += (float)riMinMag;
			table.initNewLine();
			table.addColumn("**M≥"+(float)riMinMag+"**");
			table.addColumn("![Subsection Scatter]("+outputDir.getName()+"/"+prefix+"_scatter.png)");
			table.addColumn("![Subsection 2-D Hist]("+outputDir.getName()+"/"+prefix+"_hist2D.png)");
			table.finalizeLine();
		}
		lines.addAll(table.build());
		
		if (plots.isEmpty())
			return lines;
		
		for (AbstractPlot p : plots)
			p.setPlotSize(650, 600);
		
		Loader l = loader().minMag(minMag).skipYears(skipYears);
		
		Iterable<RSQSimEvent> iterable = l.iterable();
		
		for (RSQSimEvent e : iterable)
			for (AbstractPlot p : plots)
				p.processEvent(e);
		
		for (AbstractPlot p : plots)
			p.finalizePlot();
		
		return lines;
	}

	@Override
	public Element toXMLMetadata(Element root) {
		Element el = root.addElement(XML_METADATA_NAME);
		
		el.addAttribute("name", name);
		el.addAttribute("author", author);
		el.addAttribute("dateMillis", date.getTimeInMillis()+"");
		el.addAttribute("metadata", metadata);
		if (fm != null)
			el.addAttribute("fm", fm.name());
		if (dm != null)
			el.addAttribute("dm", dm.name());
		try {
			Map<Integer, Double> slipVels = getSlipVelocities();
			el.addAttribute("slipVel", constSlipVel+"");
			if (!Double.isFinite(constSlipVel)) {
				// write individual
				Element slipVelsEl = el.addElement("SlipVelocities");
				for (int patchID : slipVels.keySet()) {
					Element patchEl = slipVelsEl.addElement("Patch");
					patchEl.addAttribute("id", patchID+"");
					patchEl.addAttribute("velocity", slipVels.get(patchID)+"");
				}
			}
			el.addAttribute("aveArea", getAveArea()+"");
			el.addAttribute("numEvents", getNumEvents()+"");
			el.addAttribute("durationYears", getDurationYears()+"");
		} catch (Exception e) {}
		
		return root;
	}
	
	static RSQSimCatalog fromXMLMetadata(Element el) {
		String name = el.attributeValue("name");
		String author = el.attributeValue("author");
		long dateMillis = Long.parseLong(el.attributeValue("dateMillis"));
		GregorianCalendar cal = new GregorianCalendar();
		cal.setTimeInMillis(dateMillis);
		String metadata = el.attributeValue("metadata");
		FaultModels fm = null;
		if (el.attribute("fm") != null)
			fm = FaultModels.valueOf(el.attributeValue("fm"));
		DeformationModels dm = null;
		if (el.attribute("dm") != null)
			dm = DeformationModels.valueOf(el.attributeValue("dm"));
		double constSlipVel = Double.parseDouble(el.attributeValue("slipVel"));
		Map<Integer, Double> slipVels = null;
		if (!Double.isFinite(constSlipVel)) {
			Element slipVelEl = el.element("SlipVelocities");
			if (slipVelEl != null) {
				slipVels = new HashMap<>();
				for (Element patchEl : XMLUtils.getSubElementsList(slipVelEl)) {
					int patchID = Integer.parseInt(patchEl.attributeValue("id"));
					double patchVel = Double.parseDouble(patchEl.attributeValue("velocity"));
					slipVels.put(patchID, patchVel);
				}
			}
		}
		double aveArea = Double.NaN;
		if (el.attribute("aveArea") != null)
			aveArea = Double.parseDouble(el.attributeValue("aveArea"));
		int numEvents = -1;
		if (el.attribute("numEvents") != null)
			numEvents = Integer.parseInt(el.attributeValue("numEvents"));
		double durationYears = Double.NaN;
		if (el.attribute("durationYears") != null)
			durationYears = Double.parseDouble(el.attributeValue("durationYears"));
		
		RSQSimCatalog cat = new RSQSimCatalog(name, author, cal, metadata, fm, dm);
		cat.aveArea = aveArea;
		cat.numEvents = numEvents;
		cat.durationYears = durationYears;
		cat.constSlipVel = constSlipVel;
		cat.slipVels = slipVels;
		return cat;
	}
	
	public static void writeCatalogsIndex(File dir) throws IOException, DocumentException {
		// sort by date, newest first
		List<Long> times = new ArrayList<>();
		List<RSQSimCatalog> catalogs = new ArrayList<>();
		for (File subDir : dir.listFiles()) {
			if (!subDir.isDirectory())
				continue;
			File xmlFile = new File(subDir, "catalog.xml");
			if (!xmlFile.exists())
				continue;
			Document doc = XMLUtils.loadDocument(xmlFile);
			Element root = doc.getRootElement();
			Element el = root.element(XML_METADATA_NAME);
			RSQSimCatalog cat = fromXMLMetadata(el);
			times.add(cat.getDate().getTimeInMillis());
			cat.dir = subDir;
			catalogs.add(cat);
		}
		catalogs = ComparablePairing.getSortedData(times, catalogs);
		Collections.reverse(catalogs);
		
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.addLine("Date", "Name", "Duration", "Element Area", "Description");
		for (RSQSimCatalog cat : catalogs) {
			table.initNewLine();
			
			table.addColumn(dateFormat.format(cat.getDate().getTime()));
			table.addColumn("["+cat.getName()+"]("+cat.dir.getName()+"#"+MarkdownUtils.getAnchorName(cat.getName())+")");
			try {
				table.addColumn(groupedIntDF.format(cat.getDurationYears())+" yrs");
			} catch (IOException e) {
				table.addColumn("");
			}
			try {
				table.addColumn(areaDF.format(cat.getAveArea())+" km");
			} catch (IOException e) {
				table.addColumn("");
			}
			table.addColumn(cat.getMetadata());
			
			table.finalizeLine();
		}
		
		List<String> lines = new LinkedList<>();
		lines.add("# RSQSim Catalog Analysis");
		lines.add("");
		lines.addAll(table.build());
		
		MarkdownUtils.writeReadmeAndHTML(lines, dir);
	}
	
	static class CatEnumDateComparator implements Comparator<Catalogs> {

		@Override
		public int compare(Catalogs o1, Catalogs o2) {
			// reverse sorted, newest first
			return o2.catalog.getDate().compareTo(o1.catalog.getDate());
		}
		
	}
	
	public static void main(String args[]) throws IOException, DocumentException {
		File gitDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		
		boolean overwriteIndividual = true;
		boolean replot = false;
		
		File baseDir = new File("/data/kevin/simulators/catalogs");
		
		Catalogs[] cats = Catalogs.values();
		Arrays.sort(cats, new CatEnumDateComparator());
		GregorianCalendar minDate = cal(2000, 1, 1);
//		GregorianCalendar minDate = cal(2018, 2, 1);
		
		for (Catalogs cat : cats) {
//		for (Catalogs cat : new Catalogs[] { Catalogs.BRUCE_2585 }) {
			if (cat.catalog.getDate().before(minDate))
				continue;
			RSQSimCatalog catalog = cat.instance(baseDir);
			System.out.print(catalog.getName()+" ? ");
			File catGitDir = new File(gitDir, catalog.getCatalogDir().getName());
			Preconditions.checkState(catGitDir.exists() || catGitDir.mkdir());
			File xmlFile = new File(catGitDir, "catalog.xml");
			if (xmlFile.exists())
				System.out.println("exists");
			else
				System.out.println("missing");
			if (overwriteIndividual || !xmlFile.exists()) {
				System.out.println("\twriting summary");
				catalog.writeMarkdownSummary(catGitDir, true, replot);
			}
		}
		
		writeCatalogsIndex(gitDir);
	}
}
