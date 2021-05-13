package scratch.kevin.simulators;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.ElementMagRangeDescription;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.parsers.EQSIMv06FileReader;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class SimAnalysisCatLoader {
	
	private List<RuptureIdentifier> rupIdens;
	private List<Color> rupIdenColors;
	private List<? extends SimulatorEvent> events;
	
	private General_EQSIM_Tools tools;
	
	public SimAnalysisCatLoader(boolean longCat, int[] include_elems, double idenMinMag, double idenMaxMag) throws IOException {
		if (include_elems != null && include_elems.length > 0) {
			rupIdens = Lists.newArrayList();
			rupIdenColors = Lists.newArrayList();
			
			loadElemMagIdens(include_elems, rupIdens, rupIdenColors, idenMinMag, idenMaxMag);
		}
		
		init(longCat, false, idenMinMag == 7d);
	}
	
	public SimAnalysisCatLoader(boolean longCat, List<RuptureIdentifier> rupIdens, boolean m7File) throws IOException {
		this.rupIdens = rupIdens;
		init(longCat, false, m7File);
	}
	
	// geom only
	private SimAnalysisCatLoader() throws IOException {
		init(false, true, false);
	}
	
	private static List<File> dirsToCheck = Lists.newArrayList();
	static {
		dirsToCheck.add(new File("/home/kevin/Simulators"));
		dirsToCheck.add(new File("/home/scec-02/kmilner/simulators"));
	}
	
	private void init(boolean longCat, boolean noEvents, boolean m7File) throws IOException {
		File dir = null;
		for (File testDir : dirsToCheck) {
			if (testDir.exists()) {
				dir = testDir;
				break;
			}
		}
		Preconditions.checkNotNull(dir, "Simulator data files not found!");
		File geomFile = new File(dir, "ALLCAL2_1-7-11_Geometry.dat");
		System.out.println("Loading geometry...");
		tools = new General_EQSIM_Tools(geomFile);
		if (noEvents)
			return;
		String eventFileName;
		if (longCat)
			eventFileName = "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.long.barall";
		else
			eventFileName = "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.barall";
		if (m7File)
			eventFileName += ".m7";
		File eventFile = new File(dir, eventFileName);
		
		System.out.println("Loading events...");
		events = EQSIMv06FileReader.readEventsFile(eventFile, tools.getElementsList(), rupIdens);
	}
	
	public static List<SimulatorElement> loadGeomOnly() throws IOException {
		return new SimAnalysisCatLoader().tools.getElementsList();
	}
	
	static String getMagStr(double mag) {
		if (mag == Math.floor(mag))
			return (int)mag +"";
		return (float)mag+"";
	}
	
	public static void loadElemMagIdens(int[] include_elems, List<RuptureIdentifier> rupIdens, List<Color> colors,
			double minMag, double maxMag) {
		String magStr;
		if (maxMag >= 10)
			magStr = getMagStr(minMag)+"+";
		else
			magStr = getMagStr(minMag)+"=>"+getMagStr(maxMag);
		for (int elemID : include_elems) {
			String name;
			Color color;
			switch (elemID) {
			case ElementMagRangeDescription.SAF_CHOLAME_ELEMENT_ID:
				name = "SAF Cholame";
				color = Color.RED;
				break;
			case ElementMagRangeDescription.SAF_CARRIZO_ELEMENT_ID:
				name = "SAF Carrizo";
				color = Color.BLUE;
				break;
			case ElementMagRangeDescription.GARLOCK_WEST_ELEMENT_ID:
				name = "Garlock";
				color = Color.GREEN;
				break;
			case ElementMagRangeDescription.SAF_MOJAVE_ELEMENT_ID:
				name = "SAF Mojave";
				color = Color.BLACK;
				break;
			case ElementMagRangeDescription.SAF_COACHELLA_ELEMENT_ID:
				name = "SAF Coachella";
				color = Color.RED;
				break;
			case ElementMagRangeDescription.SAN_JACINTO__ELEMENT_ID:
				name = "San Jacinto";
				color = Color.CYAN;
				break;
			case ElementMagRangeDescription.ELSINORE_ELEMENT_ID:
				name = "Elsinore";
				color = Color.GRAY;
				break;
			case ElementMagRangeDescription.PUENTE_HILLS_ELEMENT_ID:
				name = "Puente Hills";
				color = Color.MAGENTA;
				break;
			case ElementMagRangeDescription.NEWP_INGL_ONSHORE_ELEMENT_ID:
				name = "Newport-Inglewood (onshore)";
				color = Color.MAGENTA;
				break;
			case ElementMagRangeDescription.CALAVERAS:
				name = "Calaveras";
				color = Color.MAGENTA;
				break;
			case ElementMagRangeDescription.HAYWARD:
				name = "Hayward";
				color = Color.MAGENTA;
				break;
			case ElementMagRangeDescription.SAF_SANTA_CRUZ:
				name = "SAF Santa Cruz";
				color = Color.MAGENTA;
				break;
			case ElementMagRangeDescription.SAF_MID_PENINSULA:
				name = "SAF Mid Peninsula";
				color = Color.MAGENTA;
				break;

			default:
				throw new IllegalStateException("Unknown elem: "+elemID);
			}
			name += " "+magStr;
			rupIdens.add(new ElementMagRangeDescription(name,
					elemID, minMag, maxMag));
			if (colors != null)
				colors.add(color);
		}
	}

	public List<RuptureIdentifier> getRupIdens() {
		return rupIdens;
	}

	public List<Color> getRupIdenColors() {
		return rupIdenColors;
	}

	public List<? extends SimulatorEvent> getEvents() {
		return events;
	}
	
	public List<SimulatorElement> getElements() {
		return tools.getElementsList();
	}
	
	public static void main(String[] args) throws IOException {
		SimAnalysisCatLoader load = new SimAnalysisCatLoader(true, null, false);
		
		System.out.println("Num elements: "+load.getElements().size());
		System.out.println("Num events: "+load.getEvents().size());
		double minMag = Double.POSITIVE_INFINITY;
		double maxMag = 0d;
		List<Double> allMags = Lists.newArrayList();
		for (SimulatorEvent e : load.getEvents()) {
			double mag = e.getMagnitude();
			if (!Double.isNaN(mag)) {
				minMag = Math.min(minMag, mag);
				maxMag = Math.max(maxMag, mag);
				allMags.add(mag);
			}
		}
		HistogramFunction hist = HistogramFunction.getEncompassingHistogram(minMag, maxMag, 0.1);
		System.out.println("Min mag: "+minMag);
		for (double mag : allMags) {
			hist.add(mag, 1d);
		}
		new GraphWindow(hist, "Mag Histogram");
	}

}
