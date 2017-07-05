package scratch.kevin.simulators.synch;

import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JPanel;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.ChartPanel;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.GraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.RuptureIdentifier;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.lowagie.text.Document;
import com.lowagie.text.DocumentException;
import com.lowagie.text.HeaderFooter;
import com.lowagie.text.Phrase;
import com.lowagie.text.pdf.DefaultFontMapper;
import com.lowagie.text.pdf.PdfContentByte;
import com.lowagie.text.pdf.PdfTemplate;
import com.lowagie.text.pdf.PdfWriter;

import scratch.kevin.markov.EmpiricalMarkovChain;
import scratch.kevin.simulators.MarkovChainBuilder;
import scratch.kevin.simulators.PeriodicityPlotter;
import scratch.kevin.simulators.SimAnalysisCatLoader;
import scratch.kevin.simulators.SynchIdens;

public class OccupancyCopulaCalculator {
	
	private RuptureIdentifier iden1;
	private List<? extends SimulatorEvent> events1;
	private RuptureIdentifier iden2;
	private List<? extends SimulatorEvent> events2;
	
	private double highResTimeDelta = 0.02; // ~7 day time discretization
//	private double highResTimeDelta = 0.1; // coarse day time discretization
	private double lowResTimeDelta = 5; // 10 years for plotting occupancy
	
	private List<int[]> highResStatesPath;
	private List<double[]> highResYearsPath;
	
	private EmpiricalMarkovChain lowResChain;
	private EvenlyDiscrXYZ_DataSet lowResOccupancy;
	
	private HistogramFunction density1;
	private HistogramFunction density2;
	
	private EvenlyDiscretizedFunc cumulative1;
	private EvenlyDiscretizedFunc cumulative2;
	
	private DefaultXY_DataSet copulaScatter;
	
	private static double[] diags_to_plot = { 0, 50, 100, 150 };
	
	public OccupancyCopulaCalculator(List<? extends SimulatorEvent> events, RuptureIdentifier iden1, RuptureIdentifier iden2) {
		this(iden1, iden1.getMatches(events), iden2, iden2.getMatches(events));
	}
	
	public OccupancyCopulaCalculator(RuptureIdentifier iden1, List<? extends SimulatorEvent> events1,
			RuptureIdentifier iden2, List<? extends SimulatorEvent> events2) {
		this.iden1 = iden1;
		this.events1 = events1;
		this.iden2 = iden2;
		this.events2 = events2;
		
		List<List<? extends SimulatorEvent>> eventLists = Lists.newArrayList();
		eventLists.add(events1);
		eventLists.add(events2);
		
		System.out.println("Calculating high res states path");
		highResStatesPath = MarkovChainBuilder.getStatesPath(highResTimeDelta, eventLists, 0d);
		System.out.println("Path has "+highResStatesPath.size()+" states");
		
		System.out.println("Calculating low res Markov/Occupancy");
		lowResChain = MarkovChainBuilder.build(lowResTimeDelta, eventLists);
		StateSpacePlotter statePlotter = new StateSpacePlotter(lowResChain, Lists.newArrayList(iden1, iden2), null);
		lowResOccupancy = statePlotter.getOccupancy(0, 1);
		
		System.out.println("Calculating densities");
		calcDensities();
		System.out.println("Marginals have "+density1.size()+" points");
		
		System.out.println("Calculating cumulatives");
		cumulative1 = getCumulativeMarginal(density1);
		cumulative2 = getCumulativeMarginal(density2);
		
		System.out.println("Calculating coplua scatter");
		calcCopulaScatter();
		System.out.println("Copula scatter has "+copulaScatter.size()+" points");
	}
	
	private void calcDensities() {
		int maxState = 0;
		for (int[] state : highResStatesPath) {
			for (int val : state)
				if (val > maxState)
					maxState = val;
		}
		double maxTime = (double)maxState*highResTimeDelta;
		
		double minX = 0.5*highResTimeDelta;
		int num = (int)(maxTime/highResTimeDelta) + 10;
		density1 = new HistogramFunction(minX, num, highResTimeDelta);
		double maxX = density1.getMaxX();
		Preconditions.checkState(maxX > maxTime);
		density2 = new HistogramFunction(minX, maxX, num);
		
		highResYearsPath = Lists.newArrayList();
		for (int[] state : highResStatesPath) {
			double time1 = state[0]*highResTimeDelta;
			double time2 = state[1]*highResTimeDelta;
			density1.add(time1, 1d);
			density2.add(time2, 1d);
			highResYearsPath.add(new double[] {time1, time2});
		}
		
		density1.normalizeBySumOfY_Vals();
		density2.normalizeBySumOfY_Vals();
	}
	
	private void calcCopulaScatter() {
		copulaScatter = new DefaultXY_DataSet();
		
		for (double[] times : highResYearsPath) {
			double t1 = times[0];
			double t2 = times[1];
			
			copulaScatter.set(getCopulaPoint(t1, t2));
		}
	}
	
	private Point2D getCopulaPoint(double t1, double t2) {
		double copulaY = cumulative1.getInterpolatedY(t1);
		double copulaX = cumulative2.getInterpolatedY(t2);
		
		return new Point2D.Double(copulaX, copulaY);
	}
	
	private ArbitrarilyDiscretizedFunc getCopulaDiagonal(double startX, double startY, double delta) {
		ArbitrarilyDiscretizedFunc diag = new ArbitrarilyDiscretizedFunc();
		double x = startX;
		double y = startY;
		
		while (x < cumulative1.getMaxX() && y < cumulative1.getMaxX()) {
			diag.set(getCopulaPoint(x, y));
			
			x += delta;
			y += delta;
		}
		
		return diag;
	}
	
	private ArbitrarilyDiscretizedFunc getOccDiagonal(double startX, double startY, double delta) {
		ArbitrarilyDiscretizedFunc diag = new ArbitrarilyDiscretizedFunc();
		double x = startX;
		double y = startY;
		
		while (x < lowResOccupancy.getMaxX() && y < lowResOccupancy.getMaxY()) {
			diag.set(x, y);
			
			x += delta;
			y += delta;
		}
		
		return diag;
	}
	
	public XY_DataSet getCopulaScatter() {
		return copulaScatter;
	}
	
	public EvenlyDiscrXYZ_DataSet getCopulaDensity(int num) {
		Preconditions.checkArgument(num > 2);
		
		double gridSpacing = 1d/num;
		
		double min = 0.5*gridSpacing;
		EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(num, num, min, min, gridSpacing);
		
		for (Point2D pt : copulaScatter) {
			double x = pt.getX();
			int xInd;
			if ((float)x == 1f)
				xInd = num-1;
			else
				xInd = xyz.getXIndex(pt.getX());
			Preconditions.checkState(xInd<num, "Bad x index: x=%s, xInd=%s, numX=%s", x, xInd, num);
			double y = pt.getY();
			int yInd;
			if ((float)y == 1f)
				yInd = num-1;
			else
				yInd = xyz.getYIndex(pt.getY());
			Preconditions.checkState(yInd<num, "Bad y index: y=%s, yInd=%s, numY=%s", y, yInd, num);
			double prev = xyz.get(xInd, yInd);
			xyz.set(xInd, yInd, 1d+prev);
		}
		
		xyz.scale(1d/xyz.getSumZ());
		
		return xyz;
	}
	
	public JPanel getCopulaComboPlot(int numCopulaBins) {
		// ranges
		Range oiRange = new Range(0d, lowResOccupancy.getMaxX()+0.5*lowResOccupancy.getGridSpacingX());
		Range copulaRange = new Range(0d, 1d);
		
		String name1 = iden1.getName();
		String name2 = iden2.getName();
		
		// Copula
		EvenlyDiscrXYZ_DataSet copula = getCopulaDensity(numCopulaBins);
		CPT cpt;
		try {
			cpt = GMT_CPT_Files.MAX_SPECTRUM.instance();
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		XYZPlotSpec copulaSpec = new XYZPlotSpec(copula, cpt.rescale(0d, copula.getMaxZ()),
				"Copula", null, null, null);
		List<XY_DataSet> copulaDiags = Lists.newArrayList();
		List<PlotCurveCharacterstics> diagChars = Lists.newArrayList();
		List<XY_DataSet> occDiags = Lists.newArrayList();
		if (diags_to_plot != null && diags_to_plot.length > 0) {
			CPT diagCPT = new CPT(0d, StatUtils.max(diags_to_plot), Color.BLACK, Color.GRAY);
			for (double diag : diags_to_plot) {
				copulaDiags.add(getCopulaDiagonal(0d, diag, 10d));
				occDiags.add(getOccDiagonal(0d, diag, 10d));
				Color c = diagCPT.getColor((float)diag);
				diagChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, c));
				if (diag != 0d) {
					copulaDiags.add(getCopulaDiagonal(diag, 0d, 10d));
					occDiags.add(getOccDiagonal(diag, 0d, 10d));
					diagChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, c));
				}
			}
			copulaSpec.setXYElems(copulaDiags);
			copulaSpec.setXYChars(diagChars);
		}
		copulaSpec.setCPTVisible(false);
		
		XYZGraphPanel copulaGP = new XYZGraphPanel();
		copulaGP.drawPlot(copulaSpec, false, false, copulaRange, copulaRange);
		ChartPanel copulaChart = copulaGP.getChartPanel();
		
		// Occupancy
		XYZPlotSpec occSpec = new XYZPlotSpec(lowResOccupancy, cpt.rescale(0d, lowResOccupancy.getMaxZ()),
				"State Occupancy", name1+" OI (years)", name2+" OI (years)", null);
		if (!occDiags.isEmpty()) {
			occSpec.setXYElems(occDiags);
			occSpec.setXYChars(diagChars);
		}
		occSpec.setCPTVisible(false);
		
		XYZGraphPanel occGP = new XYZGraphPanel();
		occGP.drawPlot(occSpec, false, false, oiRange, oiRange);
		ChartPanel occChart = occGP.getChartPanel();
		
		// Cumulative along X axis
		List<XY_DataSet> cum1Funcs = Lists.newArrayList();
		cum1Funcs.add(cumulative1);
		List<PlotCurveCharacterstics> cumChars = Lists.newArrayList();
		cumChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		PlotSpec cum1Spec = new PlotSpec(cum1Funcs, cumChars, name1+" Marginal", null, "Cumulative Density");
		PlotPreferences prefs = PlotPreferences.getDefault();
		prefs.setBackgroundColor(Color.WHITE);
		
		GraphPanel cumGP1 = new GraphPanel(prefs);
		cumGP1.drawGraphPanel(cum1Spec, false, false, oiRange, copulaRange);
		ChartPanel cumChart1 = cumGP1.getChartPanel();
		
		// Cumulative along Y axis (rotated)
		List<XY_DataSet> cum2Funcs = Lists.newArrayList();
		XY_DataSet cum2Rotated = new DefaultXY_DataSet();
		for (Point2D pt : cumulative2)
			cum2Rotated.set(pt.getY(), pt.getX());
		cum2Funcs.add(cum2Rotated);
		PlotSpec cum2Spec = new PlotSpec(cum2Funcs, cumChars, name2+" Marginal", "Cumulative Density", null);
		
		GraphPanel cumGP2 = new GraphPanel(prefs);
		cumGP2.drawGraphPanel(cum2Spec, false, false, copulaRange, oiRange);
		ChartPanel cumChart2 = cumGP2.getChartPanel();
		
		JPanel xPanel = new JPanel();
		xPanel.setLayout(new BoxLayout(xPanel, BoxLayout.X_AXIS));
		
		JPanel y1Panel = new JPanel();
		y1Panel.setLayout(new BoxLayout(y1Panel, BoxLayout.Y_AXIS));
		xPanel.add(y1Panel);
		
		JPanel y2Panel = new JPanel();
		y2Panel.setLayout(new BoxLayout(y2Panel, BoxLayout.Y_AXIS));
		xPanel.add(y2Panel);
		
		y1Panel.add(cumChart2);
		y1Panel.add(copulaChart);
		y2Panel.add(occChart);
		y2Panel.add(cumChart1);
		
		int totalSize = 900;
		int bigSize = 550;
		int smallSize = totalSize - bigSize;
		
		xPanel.setPreferredSize(new Dimension(totalSize, totalSize));
		occChart.setPreferredSize(new Dimension(bigSize, bigSize));
		copulaChart.setPreferredSize(new Dimension(smallSize, smallSize));
		cumChart1.setPreferredSize(new Dimension(bigSize, smallSize));
		cumChart2.setPreferredSize(new Dimension(smallSize, bigSize));
		
		xPanel.setSize(xPanel.getPreferredSize());
		xPanel.setVisible(true);
//		xPanel.doLayout();
//		for (Component c : xPanel.getComponents())
//			c.doLayout();
		layoutRecursive(xPanel);
		xPanel.setSize(xPanel.getPreferredSize());
		
		return xPanel;
	}
	
	private static void layoutRecursive(Container c) {
		c.doLayout();
		for (Component c1 : c.getComponents()) {
			c1.doLayout();
			c1.setSize(c1.getPreferredSize());
			if (c1 instanceof Container)
				layoutRecursive((Container)c1);
		}
	}
	
	public void plotCopulaCombo(int numCopulaBins) {
		JPanel panel = getCopulaComboPlot(numCopulaBins);
		JFrame frame = new JFrame();
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setContentPane(panel);
		frame.setSize(panel.getPreferredSize());
		frame.setVisible(true);
		
		while ("".isEmpty()) {
			try {
				Thread.sleep(1000);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	public static void writePlotPDF(JPanel panel, File outputFile) throws FileNotFoundException {
		int width = (int)panel.getPreferredSize().getWidth();
		int height = (int)panel.getPreferredSize().getHeight();
		Document metadataDocument = new Document(new com.lowagie.text.Rectangle(
				width, height));
		metadataDocument.addAuthor("OpenSHA");
		metadataDocument.addCreationDate();
		HeaderFooter footer = new HeaderFooter(new Phrase("Powered by OpenSHA"), true);
		metadataDocument.setFooter(footer);
		try {
			// step 2
			PdfWriter writer;

			writer = PdfWriter.getInstance(metadataDocument,
					new FileOutputStream(outputFile));
			// step 3
			metadataDocument.open();
			// step 4
			PdfContentByte cb = writer.getDirectContent();
			PdfTemplate tp = cb.createTemplate(width, height);
			Graphics2D g2d = tp.createGraphics(width, height,
					new DefaultFontMapper());
			panel.print(g2d);
			g2d.dispose();
			cb.addTemplate(tp, 0, 0);
		}
		catch (DocumentException de) {
			de.printStackTrace();
		}
		// step 5
		metadataDocument.close();
	}
	
	public static void writePlotPNG(JPanel panel, File outputFile) throws IOException {
		int width = (int)panel.getPreferredSize().getWidth();
		int height = (int)panel.getPreferredSize().getHeight();
		BufferedImage bi = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB); 
		Graphics g = bi.createGraphics();
		panel.paint(g);
//		panel.print(g);
		g.dispose();
		ImageIO.write(bi, "png", outputFile);
	}
	
	private static EvenlyDiscretizedFunc getCumulativeMarginal(EvenlyDiscretizedFunc marginal) {
		double delta = marginal.getDelta();
		double min = marginal.getMinX()-0.5*delta;
//		if ((float)min == 0f)
//			min = 0d;
		EvenlyDiscretizedFunc cumulative = new EvenlyDiscretizedFunc(min, marginal.size(), delta);
		Preconditions.checkState(min == 0d, "Bad min: %s, marginal min: %s, delta: %s", min, marginal.getMinX(), delta);
		
		double val = 0;
		for (int i=0; i<marginal.size(); i++) {
			val += marginal.getY(i);
			cumulative.set(i, val);
		}
		
		cumulative.scale(1d/cumulative.getMaxY());
		
		return cumulative;
	}
	
	public static void main(String[] args) throws IOException {
		File outputDir = new File("/home/kevin/Simulators/copula");
		List<RuptureIdentifier> idens = SynchIdens.getStandardSoCal();
		
		List<? extends SimulatorEvent> events = new SimAnalysisCatLoader(true, idens, true).getEvents();
		
		int numCopulaBins = 50;
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance();
		
		for (int i=0; i<idens.size(); i++) {
			RuptureIdentifier iden1 = idens.get(i);
			String name1 = iden1.getName();
			for (int j=i+1; j<idens.size(); j++) {
				RuptureIdentifier iden2 = idens.get(j);
				String name2 = iden2.getName();
				
				System.out.println("Doing "+name1+" and "+name2);
				
				OccupancyCopulaCalculator calc = new OccupancyCopulaCalculator(events, iden1, iden2);
//				calc.plotCopulaCombo();
				JPanel comboPlot = calc.getCopulaComboPlot(numCopulaBins);
				String prefix = "copula_"+PeriodicityPlotter.getFileSafeString(name1)
						+"_"+PeriodicityPlotter.getFileSafeString(name2);
				writePlotPNG(comboPlot, new File(outputDir, prefix+".png"));
				writePlotPDF(comboPlot, new File(outputDir, prefix+".pdf"));
				
//				EvenlyDiscrXYZ_DataSet copula = calc.getCopulaDensity(numCopulaBins);
//				
//				double maxZ = copula.getMaxZ();
//				
//				StateSpacePlotter.plot2D(copula, cpt.rescale(0d, maxZ), name1, name2, "Copula", "copula", false, outputDir);
			}
		}
	}

}
