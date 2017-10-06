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

import org.opensha.commons.util.ExceptionUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.EventIDsRupIden;
import org.opensha.sha.simulators.iden.LogicalAndRupIden;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;
import org.opensha.sha.simulators.iden.MagRangeRuptureIdentifier;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.parsers.RSQSimFileReader;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimStateTransitionFileReader;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

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
		builder.addLine("Slim Velocity", (float)slipVel+" m/s");
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
		
		MarkdownUtils.writeReadmeAndHTML(lines, dir);
	}

	private static File getGeomFile(File dir) throws FileNotFoundException {
		for (File file : dir.listFiles()) {
			String name = file.getName().toLowerCase();
			if (name.endsWith(".flt"))
				return file;
			if (name.startsWith("zfault") && name.endsWith(".in"))
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
	
	public RSQSimEvent loadEventByID(int eventID) throws IOException {
		List<RSQSimEvent> events = loadEvents(new EventIDsRupIden(new int[] {eventID}));
		Preconditions.checkState(events.size() == 1);
		return events.get(0);
	}
	
	public List<RSQSimEvent> loadEventsByID(int... eventIDs) throws IOException {
		return loadEvents(new EventIDsRupIden(eventIDs));
	}
	
	public List<RSQSimEvent> loadEventsByMag(double minMag) throws IOException {
		return loadEvents(new MagRangeRuptureIdentifier(minMag, Double.POSITIVE_INFINITY));
	}
	
	public List<RSQSimEvent> loadEvents(RuptureIdentifier loadIden) throws IOException {
		return RSQSimFileReader.readEventsFile(dir, getElements(), Lists.newArrayList(loadIden));
	}
	
	public List<RSQSimEvent> loadEventsLogicalAnd(RuptureIdentifier... loadIdens) throws IOException {
		return loadEvents(new LogicalAndRupIden(loadIdens));
	}
	
	public List<RSQSimEvent> loadEventsLogicalOr(RuptureIdentifier... loadIdens) throws IOException {
		return loadEvents(new LogicalOrRupIden(loadIdens));
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
	
	public RSQSimEventSlipTimeFunc getSlipTimeFunc(RSQSimEvent event) throws IOException {
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
	
	public EqkRupture getGMPE_Rupture(RSQSimEvent event, double minFractForInclusion) {
		List<RSQSimEvent> events = new ArrayList<>();
		events.add(event);
		List<SimulatorElement> elements;
		try {
			elements = getElements();
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		FaultSystemRupSet rupSet = RSQSimUtils.buildFaultSystemSolution(
				getU3SubSects(), elements, events, 0d, minFractForInclusion).getRupSet();
		Preconditions.checkState(rupSet.getNumRuptures() == 1);
		
		RuptureSurface surf = rupSet.getSurfaceForRupupture(0, 1d, false);
		
		return new EqkRupture(event.getMagnitude(), rupSet.getAveRakeForRup(0), surf, null);
	}

}
