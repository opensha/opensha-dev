package scratch.kevin.simulators;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashSet;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.MagRangeRuptureIdentifier;
import org.opensha.sha.simulators.iden.RegionIden;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.parsers.EQSIMv06FileReader;
import org.opensha.sha.simulators.parsers.RSQSimFileReader;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.griddedSeismicity.AbstractGridSourceProvider;
import scratch.UCERF3.griddedSeismicity.GridSourceProvider;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.aftershockStatistics.AftershockStatsCalc;

import org.opensha.sha.simulators.EQSIM_Event;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;

public class MFDCalc {
	
	public static IncrementalMagFreqDist calcMFD(
			List<? extends SimulatorEvent> events, HashSet<Integer> elementsInRegion,
			double duration, double minMag, int num, double delta) {
		
		IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(minMag, num, delta);
		double myMin = minMag-0.5*delta;
		double myMax = mfd.getMaxX()+0.5*delta;
		for (SimulatorEvent e : events) {
			double mag = e.getMagnitude();
			if (mag < myMin || mag > myMax)
				continue;
			int ind = mfd.getClosestXIndex(mag);
			double eventRate;
			if (elementsInRegion == null)
				eventRate = 1;
			else
				eventRate = getFractInsideRegion(e, elementsInRegion);
			mfd.set(ind, mfd.getY(ind)+eventRate);
		}
		if (duration > 0)
			for (int i=0; i<mfd.size(); i++)
				mfd.set(i, mfd.getY(i) / duration);
		
		return mfd;
	}
	
	public static double getFractInsideRegion(SimulatorEvent e, HashSet<Integer> elementsInRegion) {
		double tot = e.getNumElements();
		double inside = 0d;
		for (int elemID : e.getAllElementIDs())
			if (elementsInRegion.contains(elemID))
				inside++;
		return inside / tot;
	}
	
	public static HashSet<Integer> getElementsInsideRegion(
			List<SimulatorElement> elements, Region region) {
		HashSet<Integer> elementsInRegion = new HashSet<Integer>();
		for (SimulatorElement elem : elements) {
			double lat = 0; 
			double lon = 0;
			int num = 0;
			// just averaging to get middle, should be fine for this use
			for (Location loc : elem.getVertices()) {
				lat += loc.getLatitude();
				lon += loc.getLongitude();
				num++;
			}
			lat /= (double)num;
			lon /= (double)num;
			if (region.contains(new Location(lat, lon)))
				elementsInRegion.add(elem.getID());
		}
		return elementsInRegion;
	}
	
	public static double calcMinBelow(double eventMin, double delta) {
		double min;
		for (min=0; min<=eventMin; min+=delta);
		return min - delta;
	}
	
	public static int calcNum(double min, double eventMax, double delta) {
		int num = 0;
		for (double max=min; max<=eventMax; max+=delta)
			num++;
		return num;
	}
	
	public static void writeMFDPlots(List<SimulatorElement> elements, List<? extends SimulatorEvent> events, File outputDir,
			Region... regions) throws IOException {
		if (regions.length == 0)
			regions = new Region[] { null };
		
		double duration = events.get(events.size()-1).getTimeInYears() - events.get(0).getTimeInYears();
		double minMag = 10d;
		for (SimulatorEvent e : events)
			minMag = Math.min(minMag, e.getMagnitude());
		// round minMag
		minMag = (int)(minMag*100d + 0.5)/100d;
		double delta = 0.1;
		int num = (int)((8.5d - minMag)/delta + 0.5);
		
		for (Region reg : regions) {
			HashSet<Integer> elementsInRegion;
			if (reg == null)
				elementsInRegion = null;
			else
				elementsInRegion = getElementsInsideRegion(elements, reg);
			IncrementalMagFreqDist mfd = calcMFD(events, elementsInRegion, duration, minMag+0.5*delta, num, delta);
			
			double bVal = estimateB(events, minMag);
			System.out.println("Estimated b-value: "+bVal);
			
			double[] bVals = { bVal, 1d, 1.5d, 2d };
			
			List<DiscretizedFunc> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			funcs.add(mfd);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			mfd.setName("Catalog");
			
			for (double b : bVals) {
				GutenbergRichterMagFreqDist grMFD = getGR(b, mfd);
				funcs.add(grMFD);
				if (b == bVal)
					chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLUE));
				else
					chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.GRAY));
				grMFD.setName("G-R b="+bValDF.format(b));
			}
			
			String title = (int)(duration+0.5)+" yr MFD";
			if (reg != null)
				title += ", "+reg.getName();
			String xAxisLabel = "Magnitude";
			String yAxisLabel = "Incremental Rate (1/yr)";
			PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
			spec.setLegendVisible(true);
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setTickLabelFontSize(18);
			gp.setAxisLabelFontSize(20);
			gp.setPlotLabelFontSize(21);
			gp.setBackgroundColor(Color.WHITE);
			
			gp.setUserBounds(minMag, mfd.getMaxX()+0.5*delta, 1e-7, 1e1);
			gp.drawGraphPanel(spec, false, true);
			gp.getChartPanel().setSize(1000, 800);
			
			String prefix = "mfd";
			if (reg == null)
				prefix += "_all";
			else
				prefix += "_"+reg.getName().replaceAll(" ", "_");
			File file = new File(outputDir, prefix);
			
			gp.saveAsPDF(file.getAbsolutePath()+".pdf");
			gp.saveAsPNG(file.getAbsolutePath()+".png");
			gp.saveAsTXT(file.getAbsolutePath()+".txt");
		}
	}
	
	private static final DecimalFormat bValDF = new DecimalFormat("0.0#");
	
	public static void writeMFDPlots(List<SimulatorElement> elements, File outputDir, String prefix, Region region,
			double minMag, FaultSystemSolution fssForComparison, File... eventFiles) throws IOException {
		double delta = 0.1;
		int num = (int)((8.5d - minMag)/delta + 0.5);
		double maxX = minMag + (num-1)*delta;
		
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		List<Color> colors = Lists.newArrayList(Color.RED, Color.GREEN.darker(), Color.BLUE.darker(),
				Color.ORANGE.darker(), Color.MAGENTA, Color.BLACK);
		
		List<RuptureIdentifier> rupIdens = Lists.newArrayList();
		rupIdens.add(new MagRangeRuptureIdentifier(minMag, 10d));
		HashSet<Integer> elementsInRegion;
		if (region != null) {
			rupIdens.add(new RegionIden(region));
			elementsInRegion = getElementsInsideRegion(elements, region);
		} else {
			elementsInRegion = null;
		}
		
		if (fssForComparison != null) {
			IncrementalMagFreqDist compSupraSeis = fssForComparison.calcTotalNucleationMFD(
					AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF, 8.55, 0.1);
			compSupraSeis.setName("UCERF3 Supra-Seis");
			SummedMagFreqDist compSubSeis = null;
			SummedMagFreqDist compTotal = null;
			if (fssForComparison.getGridSourceProvider() != null) {
				GridSourceProvider prov = fssForComparison.getGridSourceProvider();
				GriddedRegion reg = prov.getGriddedRegion();
				compSubSeis = new SummedMagFreqDist(compSupraSeis.getMinX(), compSupraSeis.size(), compSupraSeis.getDelta());
				compTotal = new SummedMagFreqDist(compSupraSeis.getMinX(), compSupraSeis.size(), compSupraSeis.getDelta());
				for (int index=0; index<prov.size(); index++) {
					if (region == null || region.contains(reg.getLocation(index))) {
						IncrementalMagFreqDist associatedMFD = prov.getNodeSubSeisMFD(index);
						if (associatedMFD != null)
							compSubSeis.addResampledMagFreqDist(associatedMFD, true);
						compTotal.addResampledMagFreqDist(prov.getNodeMFD(index), true);
					}
				}
				// now add in supra
				compSubSeis.addResampledMagFreqDist(compSupraSeis, true);
				compSubSeis.setName("UCERF3 Fault");
				compTotal.addResampledMagFreqDist(compSupraSeis, true);
				compTotal.setName("UCERF3 Total");
				
				funcs.add(compSubSeis);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.DARK_GRAY));
				funcs.add(compSupraSeis);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.DARK_GRAY));
				funcs.add(compTotal);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.DARK_GRAY));
			} else {
				funcs.add(compSupraSeis);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.DARK_GRAY));
			}
		}
		
		double meanCumRate = 0d;
		
		for (int i=0; i<eventFiles.length; i++) {
			File eventFile = eventFiles[i];
			
			System.out.println(eventFile.getName());
			
			List<? extends SimulatorEvent> events;
			if (eventFile.isDirectory() || eventFile.getName().endsWith(".bin"))
				events = RSQSimFileReader.readEventsFile(eventFile, elements, rupIdens);
			else
				events = EQSIMv06FileReader.readEventsFile(eventFile, elements, rupIdens);
			
			double duration = events.get(events.size()-1).getTimeInYears() - events.get(0).getTimeInYears();
			
			IncrementalMagFreqDist mfd = calcMFD(events, elementsInRegion, duration, minMag+0.5*delta, num, delta);
			meanCumRate += mfd.getCumRateDistWithOffset().getY(0)/(double)eventFiles.length;
			
			Color c = colors.get(i % colors.size());
			
			double bVal = estimateB(events, minMag);
			System.out.println("Estimated b-value: "+bVal);
			String name = eventFile.getName();
			if (name.contains("."))
				name = name.substring(0, name.lastIndexOf("."));
			mfd.setName(name+" (b="+bValDF.format(bVal)+")");
			
			funcs.add(mfd);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, c));
		}
		
		double[] bVals = { 1d, 1.5d, 2d };
		
		for (double b : bVals) {
			GutenbergRichterMagFreqDist grMFD = getGR(b, meanCumRate, minMag, maxX, num);
			funcs.add(grMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.GRAY));
			grMFD.setName("G-R b="+bValDF.format(b));
		}
		
		String title;
		if (region == null)
			title = "Total MFDs";
		else
			title = region.getName()+" MFDs";
		String xAxisLabel = "Magnitude";
		String yAxisLabel = "Incremental Rate (1/yr)";
		PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		gp.setBackgroundColor(Color.WHITE);
		
		gp.setUserBounds(minMag, maxX+0.5*delta, 1e-7, 1e1);
		gp.drawGraphPanel(spec, false, true);
		gp.getChartPanel().setSize(1000, 800);
		
		File file = new File(outputDir, prefix);
		
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsTXT(file.getAbsolutePath()+".txt");
	}
	
	private static GutenbergRichterMagFreqDist getGR(double b, IncrementalMagFreqDist mfd) {
		double totCmlRate = mfd.getCumRateDistWithOffset().getY(0);
		return getGR(b, totCmlRate, mfd.getMinX(), mfd.getMaxX(), mfd.size());
	}
	
	private static GutenbergRichterMagFreqDist getGR(double b, double totCmlRate, double minX, double maxX, int size) {
		GutenbergRichterMagFreqDist grMFD = new GutenbergRichterMagFreqDist(b, 1d, minX, maxX, size);
		grMFD.scaleToCumRate(0, totCmlRate);
		return grMFD;
	}

	public static double estimateB(List<? extends SimulatorEvent> events, double minMag) {
		double magMean = 0d;
		for (SimulatorEvent e : events)
			if (e.getMagnitude() > minMag)
				magMean += e.getMagnitude();
		magMean /= events.size();
		double magComplete = minMag;
		double magPrecision = 0d;
		
		// estimate b
		return AftershockStatsCalc.getMaxLikelihood_b_value(magMean, magComplete, magPrecision);
	}

	/**
	 * @param args
	 * @throws IOException 
	 * @throws DocumentException 
	 */
	public static void main(String[] args) throws IOException, DocumentException {
//		File internDir = new File("/home/kevin/Simulators/UCERF3_interns");
//		File outputDir = internDir;
//		File geomFile = new File(new File(internDir, "sigmahigh"), "UCERF3.D3.1.1km.tri.2.flt");
//		List<SimulatorElement> elements = RSQSimFileReader.readGeometryFile(geomFile, 11, 'S');
//		
//		File[] eventDirs = { new File(internDir, "base"), new File(internDir, "sigmahigh"),
//				new File(internDir, "sigmalow"), new File(internDir, "state")};
		
		File outputDir = new File("/home/kevin/Simulators/bruce/rundir1435");
		File geomFile = new File(outputDir, "zfault_Deepen.in");
		List<SimulatorElement> elements = RSQSimFileReader.readGeometryFile(geomFile, 11, 'S');
		File[] eventDirs = { outputDir };
		
		double minMag = 5d;
		
		FaultSystemSolution fssForComparison = FaultSystemIO.loadSol(new File(""
				+ "/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "FM3_1_GEOL_MEAN_BRANCH_AVG_SOL.zip"));
		
		writeMFDPlots(elements, outputDir, "total_mfds", null, minMag, fssForComparison, eventDirs);
		
//		File dir = new File("/home/kevin/Simulators/UCERF3_interns/UCERF3sigmahigh");
//		Region[] regions =  { new CaliforniaRegions.RELM_SOCAL(), new CaliforniaRegions.RELM_TESTING() };
//		File eventDir = dir;
//		double minMag = 5.5d;
//		List<RSQSimEvent> events = RSQSimFileReader.readEventsFile(eventDir, elements,
//				Lists.newArrayList(new MagRangeRuptureIdentifier(minMag, 10d)));
//		
//		writeMFDPlots(elements, events, new File("/tmp"), regions);
	}

}
