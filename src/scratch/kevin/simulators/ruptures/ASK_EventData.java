package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.awt.Font;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.imr.attenRelImpl.ngaw2.FaultStyle;

import com.google.common.base.Preconditions;
import com.google.common.collect.Range;

public class ASK_EventData {
	
	public final int recID;
	public final int eqID;
	public final double mag;
	public final double rake;
	public final double zTOR;
	public final double rRup;
	public final double rJB;
	public final double rX;
	public final double rY;
	public final double vs30;
	public final double rockPGA;
	public final double lnY;
	
	private ASK_EventData(int recID, int eqID, double mag, double rake, double zTOR, double rRup, double rJB, double rX,
			double rY, double vs30, double rockPGA, double lnY) {
		super();
		this.recID = recID;
		this.eqID = eqID;
		this.mag = mag;
		this.rake = rake;
		this.zTOR = zTOR;
		this.rRup = rRup;
		this.rJB = rJB;
		this.rX = rX;
		this.rY = rY;
		this.vs30 = vs30;
		this.rockPGA = rockPGA;
		this.lnY = lnY;
	}
	
	private static final Map<Float, String> periodFileMap;
	static {
		periodFileMap = new HashMap<>();
		periodFileMap.put(1f, "resources/resid_T1.000.out2.txt");
		periodFileMap.put(2f, "resources/resid_T2.000.out2.txt");
		periodFileMap.put(3f, "resources/resid_T3.000.out2.txt");
	}
	
	public static Map<Integer, List<ASK_EventData>> load(double period) throws IOException {
		Preconditions.checkState(periodFileMap.containsKey((float)period), "No data file for period: %s", (float)period);
		String fileName = periodFileMap.get((float)period);
		BufferedReader br = new BufferedReader(new InputStreamReader(ASK_EventData.class.getResourceAsStream(fileName)));
		
		br.readLine(); // skip header
		
		Map<Integer, List<ASK_EventData>> map = new HashMap<>();
		
		String line;
		while ((line = br.readLine()) != null) {
			line = line.trim();
			String[] split = line.split("\\s+");
			Preconditions.checkState(split.length == 28, "expected 29 items, got %s:\n%s", split.length, line);
			int index = 0;
			index++; // i
			index++; // j
			int recID = Integer.parseInt(split[index++]); // recID
			int eqID = Integer.parseInt(split[index++]); // eqid
			index++; // region
			index++; // SOF
			index++; // magBin
			double mag = Double.parseDouble(split[index++]); // mag
			double rake = Double.parseDouble(split[index++]); // rake
			double zTOR = Double.parseDouble(split[index++]); // ZTOR
			double rRup = Double.parseDouble(split[index++]); // Rrup
			double rJB = Double.parseDouble(split[index++]); // RJB
			double rX = Double.parseDouble(split[index++]); // RX
			double rY = Double.parseDouble(split[index++]); // RY
			index++; // FLen
			double vs30 = Double.parseDouble(split[index++]); // VS30
			double rockPGA = Double.parseDouble(split[index++]); // rockPGA
			double lnY = Double.parseDouble(split[index++]); // lnY
			// intra-Resid
			// Median
			// eventTerm
			// srcSiteA
			// z1
			// dip
			// Z1_1D
			// Z1_Ref
			// iregion(i)
			// epsilon_within
			ASK_EventData data = new ASK_EventData(recID, eqID, mag, rake, zTOR, rRup, rJB, rX, rY, vs30, rockPGA, lnY);
			List<ASK_EventData> eventList = map.get(eqID);
			if (eventList == null) {
				eventList = new ArrayList<>();
				map.put(eqID, eventList);
			}
			eventList.add(data);
		}
		
		return map;
	}
	
	public static Map<Integer, List<ASK_EventData>> getMatches(Map<Integer, List<ASK_EventData>> rawData, Range<Double> magRange,
			Range<Double> distRange, FaultStyle faultStyle, double rakeTolerance) {
		Map<Integer, List<ASK_EventData>> ret = new HashMap<>();
		for (Integer eventID : rawData.keySet()) {
			List<ASK_EventData> datas = rawData.get(eventID);
			if (datas.isEmpty())
				continue;
			double mag = datas.get(0).mag;
			if (magRange != null && !magRange.contains(mag))
				continue;
			if (faultStyle != null) {
				double rake = datas.get(0).rake;
				boolean match;
				switch (faultStyle) {
				case STRIKE_SLIP:
					match = (rake >= -180 && rake <= -180+rakeTolerance) || (rake >= -rakeTolerance && rake <= rakeTolerance)
							|| (rake >= 180-rakeTolerance && rake<=180d);
					break;
				case REVERSE:
					match = rake >= 90-rakeTolerance && rake <= 90+rakeTolerance;
					break;
				case NORMAL:
					match = rake >= -90-rakeTolerance && rake <= -90+rakeTolerance;
					break;

				default:
					throw new IllegalStateException("Unsupported fault type: "+faultStyle);
				}
				if (!match) {
//					System.out.println("Event with rake "+(float)rake+" is not a match for "+faultStyle);
					continue;
				}
			}
			if (distRange == null) {
				ret.put(eventID, datas);
			} else {
				List<ASK_EventData> matches = new ArrayList<>();
				for (ASK_EventData data : datas)
					if (distRange.contains(data.rRup))
						matches.add(data);
				if (!matches.isEmpty())
					ret.put(eventID, matches);
			}
		}
		
		return ret;
	}
	
	private static MinMaxAveTracker getStats(Map<Integer, List<ASK_EventData>> datas) {
		MinMaxAveTracker track = new MinMaxAveTracker();
		for (List<ASK_EventData> list : datas.values())
			track.addValue(list.size());
		return track;
	}
	
	private static List<String> csvLine(String label, double[] periods, Map<Double, Map<Integer, List<ASK_EventData>>> periodData,
			Range<Double> magRange, Range<Double> distRange, FaultStyle faultStyle, double rakeTolerance) {
		List<String> line = new ArrayList<>();
		line.add(label);
		for (double period : periods) {
			Map<Integer, List<ASK_EventData>> data = periodData.get(period);
			
			Map<Integer, List<ASK_EventData>> allMatches = getMatches(data, magRange, distRange, faultStyle, rakeTolerance);
			MinMaxAveTracker allTrack = getStats(allMatches);
			
			line.add(allMatches.size()+"");
			if (allMatches.isEmpty()) {
				line.add("");
				line.add("");
				line.add("");
				line.add("");
			} else {
				line.add(twoDigit.format(allTrack.getAverage()));
				line.add((int)allTrack.getMin()+"");
				line.add((int)allTrack.getMax()+"");
				line.add((int)allTrack.getSum()+"");
			}
			line.add("");
		}
		
		return line;
	}
	
	public static void plotCountHist(File outputDir, String prefix,
			Map<Integer, List<ASK_EventData>> data) throws IOException {
		plotMultiCountHist(outputDir, prefix, "Event Recordings Distribution",
				data.values(), null, null, null);
	}
	
	private static HistogramFunction countsHist(Collection<List<ASK_EventData>> data, int max) {
		double delta;
		if (max > 100)
			delta = 10;
		else if (max < 50)
			delta = 5d;
		else if (max > 10)
			delta = 2d;
		else
			delta = 1d;
		HistogramFunction hist = HistogramFunction.getEncompassingHistogram(1d, (double)max, delta);
		
		for (List<ASK_EventData> list : data)
			hist.add(hist.getClosestXIndex(list.size()), 1d);
		
		return hist;
	}
	
	public static void plotMultiCountHist(File outputDir, String prefix, String title,
			Collection<List<ASK_EventData>> data1, String label1,
			Collection<List<ASK_EventData>> data2, String label2) throws IOException {
		int max = 0;
		for (List<ASK_EventData> data : data1)
			max = Integer.max(max, data.size());
		if (data2 != null)
			for (List<ASK_EventData> data : data2)
				max = Integer.max(max, data.size());
		max = (int)(10d*Math.ceil(max/10d));
		
		HistogramFunction hist1 = countsHist(data1, max);
		PlotCurveCharacterstics pChar = new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK);
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(hist1);
		chars.add(pChar);
		
		PlotSpec spec1 = new PlotSpec(funcs, chars, title, "Recordings", "Count");
		
		Font labelFont = new Font(Font.SANS_SERIF, Font.BOLD, 20);
		double annX = max*0.95;

		List<org.jfree.data.Range> yRanges = new ArrayList<>();
		double maxY1 = hist1.getMaxY()*1.2;
		yRanges.add(new org.jfree.data.Range(0d, maxY1));
		if (label1 != null) {
			double annY = maxY1*0.95;
			XYTextAnnotation ann = new XYTextAnnotation(label1, annX, annY);
			ann.setFont(labelFont);
			ann.setTextAnchor(TextAnchor.TOP_RIGHT);
			List<XYTextAnnotation> anns = new ArrayList<>();
			anns.add(ann);
			spec1.setPlotAnnotations(anns);
		}
		
		List<PlotSpec> specs = new ArrayList<>();
		specs.add(spec1);
		
		if (data2 != null) {
			HistogramFunction hist2 = countsHist(data2, max);
			funcs = new ArrayList<>();
			chars = new ArrayList<>();
			
			funcs.add(hist2);
			chars.add(pChar);
			
			PlotSpec spec2 = new PlotSpec(funcs, chars, title, "Recordings", "Count");
			
			double maxY2 = hist2.getMaxY()*1.2;
			yRanges.add(new org.jfree.data.Range(0d, maxY2));
			if (label2 != null) {
				double annY = maxY2*0.95;
				XYTextAnnotation ann = new XYTextAnnotation(label2, annX, annY);
				ann.setFont(labelFont);
				ann.setTextAnchor(TextAnchor.TOP_RIGHT);
				List<XYTextAnnotation> anns = new ArrayList<>();
				anns.add(ann);
				spec2.setPlotAnnotations(anns);
			}
			
			specs.add(spec2);
		}
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(28);
		gp.setBackgroundColor(Color.WHITE);
		
		org.jfree.data.Range xRange = new org.jfree.data.Range(0, (double)max);
		List<org.jfree.data.Range> xRanges = new ArrayList<>();
		xRanges.add(xRange);
		
		gp.drawGraphPanel(specs, false, false, xRanges, yRanges);
		
		File file = new File(outputDir, prefix);
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
	}
	
	private static final DecimalFormat twoDigit = new DecimalFormat("0.00");
	
	public static void main(String[] args) throws IOException {
		double[] periods = { 1d, 2d, 3d };
		
		List<Range<Double>> distRanges = new ArrayList<>();
		distRanges.add(Range.closed(10d, 30d));
		distRanges.add(Range.closed(40d, 60d));
		distRanges.add(Range.closed(80d, 120d));
		
		List<Range<Double>> magRanges = new ArrayList<>();
//		magRanges.add(Range.closed(6.5, 6.7));
		magRanges.add(Range.closed(6.4, 6.8));
//		magRanges.add(Range.closed(7.1, 7.3));
		magRanges.add(Range.closed(7.0, 7.4));
		
		List<FaultStyle> faultStyles = new ArrayList<>();
		faultStyles.add(FaultStyle.STRIKE_SLIP);
		faultStyles.add(FaultStyle.REVERSE);
		double rakeTolerance = 30;
		
		CSVFile<String> csv = new CSVFile<>(false);
		List<String> header = new ArrayList<>();
		header.add("");
		for (double period : periods) {
			String periodStr = Math.round(period) == period ? (int)period+"s" : (float)period+"s";
			header.add(periodStr+" Events");
			header.add(periodStr+" Ave. Rec Per Event");
			header.add(periodStr+" Min Rec Per Event");
			header.add(periodStr+" Max Rec Per Event");
			header.add(periodStr+" Total # Records");
			header.add("");
		}
		csv.addLine(header);
		
		Map<Double, Map<Integer, List<ASK_EventData>>> periodData = new HashMap<>();
		for (double period : periods) {
			System.out.println("Loading data for "+(float)period+"s");
			Map<Integer, List<ASK_EventData>> data = load(period);
			System.out.println("Loaded data for "+data.size()+" events");
			int totNum = 0;
			for (int eqID : data.keySet())
				totNum += data.get(eqID).size();
			System.out.println("Loaded "+totNum+" total records");
			System.out.println();
			periodData.put(period, data);
		}
		
		for (Range<Double> magRange : magRanges) {
			for (FaultStyle style : faultStyles) {
				csv.addLine(csvLine("M"+magRange.lowerEndpoint().floatValue()+"-"+magRange.upperEndpoint().floatValue()+", "+style,
						periods, periodData, magRange, null, style, rakeTolerance));
				
				for (Range<Double> distRange : distRanges) {
					csv.addLine(csvLine("Distance range ["+distRange.lowerEndpoint().intValue()+" "+distRange.upperEndpoint().intValue()+"] km",
							periods, periodData, magRange, distRange, style, rakeTolerance));
				}
				
				csv.addLine("");
			}
		}
		
		csv.writeToFile(new File("/tmp/ask_event_counts.csv"));
	}

}
