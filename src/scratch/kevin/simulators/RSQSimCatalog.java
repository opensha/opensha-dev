package scratch.kevin.simulators;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.opensha.commons.util.ExceptionUtils;
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
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Preconditions;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.kevin.MarkdownUtils;
import scratch.kevin.MarkdownUtils.TableBuilder;

public class RSQSimCatalog {
	
	public enum Catalogs {
		JG_UCERF3_millionElement("JG_UCERF3_millionElement", "U3 1mil Element Test", "Jacqui Gilchrist", cal(2017, 9, 27),
				"Test 1 million element catalog on UCERF3 fault system, ~0.25 km^2 trianglar elements",
				FaultModels.FM3_1, DeformationModels.GEOLOGIC, 1d),
		BRUCE_2194_LONG("rundir2194_long", "Bruce 2194 Long", "Bruce Shaw (extended by Jacqui Gilchrist)", cal(2017, 8, 31),
				"Catalog with decent large event scaling and distribution of sizes while not using"
				+ " any of the enhanced frictional weakening terms.", FaultModels.FM3_1, DeformationModels.GEOLOGIC, 1d);
		
		private String dirName;
		private RSQSimCatalog catalog;
		
		private Catalogs(String dirName, String name, String author, GregorianCalendar date, String metadata,
				FaultModels fm, DeformationModels dm, double slipVel) {
			this.dirName = dirName;
			catalog = new RSQSimCatalog(name, author, date, metadata, fm, dm, slipVel);
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
	private double slipVel;
	
	private List<SimulatorElement> elements;
	private RSQSimStateTransitionFileReader transReader;
	private List<FaultSectionPrefData> subSects;
	private Map<Integer, Double> subSectAreas;
	
	private RSQSimCatalog(String name, String author, GregorianCalendar date, String metadata,
			FaultModels fm, DeformationModels dm, double slipVel) {
		this(null, name, author, date, metadata, fm, dm, slipVel);
	}

	public RSQSimCatalog(File dir, String name, String author, GregorianCalendar date, String metadata,
			FaultModels fm, DeformationModels dm, double slipVel) {
		this.dir = dir;
		this.name = name;
		this.author = author;
		this.date = date;
		this.metadata = metadata;
		this.fm = fm;
		this.dm = dm;
		this.slipVel = slipVel;
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
	
	public double getSlipVelocity() {
		return slipVel;
	}
	
	private static DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd");
	private static DecimalFormat areaDF = new DecimalFormat("0.00");
	
	public List<String> getMarkdownMetadataTable() {
		TableBuilder builder = MarkdownUtils.tableBuilder();
		builder.addLine("Catalog", getName());
		builder.addLine("Author", getAuthor()+", "+dateFormat.format(getDate().getTime()));
		builder.addLine("Description", getMetadata());
		builder.addLine("Fault/Def Model", fm+", "+dm);
		builder.addLine("Slip Velocity", (float)slipVel+" m/s");
		try {
			double aveArea = 0d;
			for (SimulatorElement e : getElements())
				aveArea += e.getArea()*1e-6;
			aveArea /= getElements().size();
			builder.addLine("Average Element Area", areaDF.format(aveArea)+" km^2");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return builder.build();
	}
	
	public void writeMarkdownSummary(File dir) throws IOException {
		List<String> lines = new LinkedList<>();
		lines.add("# "+getName());
		lines.addAll(getMarkdownMetadataTable());
		
		List<String> eventLinks = new ArrayList<>();
		List<String> eventNames = new ArrayList<>();
		
		List<String> gmpeLinks = new ArrayList<>();
		List<String> gmpeNames = new ArrayList<>();
		
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
			}
		}
		
		if (!eventNames.isEmpty()) {
			lines.add("");
			lines.add("## Event Comparisons");
			for (int i=0; i<eventNames.size(); i++)
				lines.add("* ["+eventNames.get(i)+"]("+eventLinks.get(i)+"/)");
		}
		if (!gmpeNames.isEmpty()) {
			lines.add("");
			lines.add("## GMPE Comparisons");
			for (int i=0; i<gmpeNames.size(); i++)
				lines.add("* ["+gmpeNames.get(i)+"]("+gmpeLinks.get(i)+"/)");
		}
		
		MarkdownUtils.writeReadmeAndHTML(lines, dir);
	}

	private static File getGeomFile(File dir) throws FileNotFoundException {
		for (File file : dir.listFiles()) {
			String name = file.getName().toLowerCase();
			if (name.endsWith(".flt"))
				return file;
			if (name.startsWith("zfault") && name.endsWith(".in") && !name.contains("Deepen_"))
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
		return new RSQSimEventSlipTimeFunc(getTransitions().getTransitions(event), slipVel);
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
	
	public EqkRupture getGMPE_Rupture(RSQSimEvent event, double minFractForInclusion) {
		List<RSQSimEvent> events = new ArrayList<>();
		events.add(event);
		List<SimulatorElement> elements;
		Map<Integer, Double> subSectAreas = null;
		try {
			elements = getElements();
			if (minFractForInclusion > 0d)
				subSectAreas = getSubSectAreas();
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		return RSQSimUtils.buildSubSectBasedRupture(event, getU3SubSects(), elements, minFractForInclusion, subSectAreas);
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

}
