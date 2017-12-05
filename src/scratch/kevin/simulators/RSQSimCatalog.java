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
import org.opensha.commons.metadata.XMLSaveable;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.XMLUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.EventIDsRupIden;
import org.opensha.sha.simulators.iden.LogicalAndRupIden;
import org.opensha.sha.simulators.iden.MagRangeRuptureIdentifier;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.iden.SkipYearsLoadIden;
import org.opensha.sha.simulators.parsers.RSQSimFileReader;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimStateTransitionFileReader;
import org.opensha.sha.simulators.srf.RSQSimTransValidIden;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.IDPairing;
import scratch.kevin.MarkdownUtils;
import scratch.kevin.MarkdownUtils.TableBuilder;
import scratch.kevin.simulators.plots.AbstractPlot;
import scratch.kevin.simulators.plots.MFDPlot;
import scratch.kevin.simulators.plots.MagAreaScalingPlot;
import scratch.kevin.simulators.plots.RecurrenceIntervalPlot;
import scratch.kevin.simulators.plots.RuptureVelocityPlot;
import scratch.kevin.simulators.plots.StationarityPlot;

public class RSQSimCatalog implements XMLSaveable {
	
	public enum Catalogs {
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
		lines.add("");;
		int tocIndex = lines.size();
		
		List<String> eventLinks = new ArrayList<>();
		List<String> eventNames = new ArrayList<>();
		
		List<String> gmpeLinks = new ArrayList<>();
		List<String> gmpeNames = new ArrayList<>();
		
		List<String> gmpeRGLinks = new ArrayList<>();
		List<String> gmpeRGNames = new ArrayList<>();
		
		String rotDDLink = null;
		
		for (File subDir : dir.listFiles()) {
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
				gmpeNames.add(name.substring("gmpe_bbp_comparisons_".length()));
				gmpeLinks.add(name);
			} else if (name.startsWith("gmpe_bbp_rg_comparisons_")) {
				gmpeRGNames.add(name.substring("gmpe_bbp_rg_comparisons_".length()));
				gmpeRGLinks.add(name);
			} else if (name.equals("catalog_rotd_ratio_comparisons")) {
				Preconditions.checkState(rotDDLink == null, "Duplicate RotDD dirs! %s and %s", name, rotDDLink);
				rotDDLink = name;
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
		if (!gmpeNames.isEmpty()) {
			lines.add("");
			lines.add("## Full Catalog GMPE Comparisons");
			lines.add(topLink);
			lines.add("");
			for (int i=0; i<gmpeNames.size(); i++)
				lines.add("* ["+gmpeNames.get(i)+"]("+gmpeLinks.get(i)+"/)");
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
		
		if (plots) {
			File resourcesDir = new File(dir, "resources");
			Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
			lines.add("");
			lines.addAll(writeStandardDiagnosticPlots(resourcesDir, 5000, 6d, replot, topLink));
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

	private static File getGeomFile(File dir) throws FileNotFoundException {
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
			File geomFile = getGeomFile(dir);
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
		throw new FileNotFoundException("No geometry file found in "+dir.getAbsolutePath());
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
	
	private synchronized Map<Integer, Double> getSubSectAreas() throws IOException {
		if (subSectAreas == null)
			subSectAreas = RSQSimUtils.calcSubSectAreas(getElements());
		return subSectAreas;
	}
	
	private synchronized Map<IDPairing, Double> getSubSectDistsCache() {
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
			loadIdens.add(new SkipYearsLoadIden(years));
			return this;
		}
		
		public Loader matches(RuptureIdentifier iden) {
			loadIdens.add(iden);
			return this;
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
		table.addColumn("![MFD Scatter]("+outputDir.getName()+"/mag_area.png)");
		table.addColumn("![MFD Hist]("+outputDir.getName()+"/mag_area_hist2D.png)");
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
		lines.add("### Interevent-Time Distributions");
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
		
		if (replot || !new File(outputDir, "stationarity.png").exists()) {
			StationarityPlot stationarityPlot = new StationarityPlot(minMag, 7d);
			stationarityPlot.initialize(getName(), outputDir, "stationarity");
			plots.add(stationarityPlot);
		}
		lines.add("### Stationarity Plot");
		lines.add(topLink);
		lines.add("");
		lines.add("![Stationarity]("+outputDir.getName()+"/stationarity.png)");
		
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
	
	private static class CatEnumDateComparator implements Comparator<Catalogs> {

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
		
		for (Catalogs cat : cats) {
//		for (Catalogs cat : new Catalogs[] { Catalogs.BRUCE_2310 }) {
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
