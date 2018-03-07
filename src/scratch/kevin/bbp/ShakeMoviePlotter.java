package scratch.kevin.bbp;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipException;

import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.DocumentException;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.axis.TickUnit;
import org.jfree.chart.axis.TickUnits;
import org.jfree.chart.title.PaintScaleLegend;
import org.jfree.data.Range;
import org.jfree.ui.RectangleEdge;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.Vertex;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;

import com.google.common.base.Preconditions;

import scratch.kevin.bbp.BBP_SimZipLoader.BBP_ShakeMapSimZipLoader;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class ShakeMoviePlotter {
	
	public static int PLOT_WIDTH = 900;
	public static int PLOT_HEIGHT = 800;
	public static double TIME_BUFFER = 45; // seconds after end of rupture
	
	public static Color CA_OUTLINE_COLOR = Color.DARK_GRAY;
	
	public static void plotShakeMovie(BBP_ShakeMapSimZipLoader loader, GriddedRegion gridReg, List<BBP_Site> sites, RSQSimEvent event,
			RSQSimEventSlipTimeFunc slipTimeFunc, String title, File outputDir, String prefix, boolean acceleration, double maxFPS)
					throws IOException {
		Preconditions.checkState(gridReg.getNodeCount() == sites.size());
		
		DiscretizedFunc[] seismograms = new DiscretizedFunc[sites.size()];
		
		double maxZ = 0;
		
		for (int i=0; i<sites.size(); i++) {
			DiscretizedFunc[] mySeis;
			if (acceleration)
				mySeis = loader.readAccelSeis(i);
			else
				mySeis = loader.readVelSeis(i);
			
			DiscretizedFunc nsFunc = mySeis[0];
			DiscretizedFunc ewFunc = mySeis[1];
			Preconditions.checkState(nsFunc.size() == ewFunc.size());
			
			double[] xs = new double[nsFunc.size()];
			double[] ys = new double[nsFunc.size()];
			
			for (int j=0; j<xs.length; j++) {
				xs[j] = nsFunc.getX(j);
				// geometric mean horizontal shaking
				ys[j] = StatUtils.geometricMean(new double[] { Math.abs(nsFunc.getY(j)), Math.abs(ewFunc.getY(j)) });
//				ys[j] = Math.abs(nsFunc.getY(j));
			}
			
			seismograms[i] = new LightFixedXFunc(xs, ys);
			
			maxZ = Math.max(maxZ, seismograms[i].getMaxY());
		}
		
		// decide on actual framerate
		double dataDeltaX = seismograms[0].getX(1) - seismograms[0].getX(0);
		double dataFramerate = 1d/dataDeltaX;
		System.out.println("Data framerate: "+dataFramerate);
		double dataModDouble = dataFramerate/maxFPS;
		int dataMod = (int)Math.round(dataModDouble);
		if (dataMod <= 0)
			dataMod = 1;
		System.out.println("Data mod: "+dataModDouble+" = "+dataMod);
		double framerate = dataFramerate / dataMod;
		System.out.println("Output framerate: "+framerate);
		
		boolean logZ = true;
		
		CPT seisCPT = GMT_CPT_Files.GMT_HOT.instance().reverse();
		if (acceleration) {
			seisCPT = seisCPT.rescale(-2, 3);
		} else {
			seisCPT = seisCPT.rescale(-3, 2);
		}
		seisCPT.setNanColor(seisCPT.getMinColor());
		seisCPT.setBelowMinColor(seisCPT.getMinColor());
		seisCPT.setAboveMaxColor(seisCPT.getMaxColor());
		
		List<RenderThread> threads = new ArrayList<>();
		
		slipTimeFunc = slipTimeFunc.asRelativeTimeFunc();
		double maxTime = slipTimeFunc.getEndTime() + TIME_BUFFER;
		
		int numFrames = 0;
		ArrayDeque<Integer> frameDeque = new ArrayDeque<>();
		for (int i=0; i<seismograms[0].size(); i+=dataMod) {
			double time = seismograms[0].getX(i);
			if (time > maxTime)
				break;
			frameDeque.add(numFrames++);
		}
		
		for (int t=0; t<Runtime.getRuntime().availableProcessors(); t++)
			threads.add(new RenderThread(seismograms, dataMod, seisCPT, logZ, outputDir, prefix,
					acceleration, title, gridReg, event, slipTimeFunc, numFrames, frameDeque));
		
		for (RenderThread t : threads)
			t.start();
		
		for (RenderThread t : threads) {
			try {
				t.join();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		FileWriter fw = new FileWriter(new File(outputDir, prefix+"_build_movie.sh"));
		
		fw.write("#!/bin/bash\n");
		fw.write("\n");
		fw.write("ffmpeg -framerate "+(float)framerate+" -pattern_type glob -i '"+prefix+"*.png' "
				+ "-c:v libx264 -r 30 -pix_fmt yuv420p "+prefix+".mp4\n");
		
		fw.close();
	}
	
	private static class RenderThread extends Thread {
		
		private int numFrames;
		private ArrayDeque<Integer> frameDeque;
		private DiscretizedFunc[] seismograms;
		private int dataMod;
		private CPT seisCPT;
		private boolean logZ;
		private File outputDir;
		private String prefix;
		private boolean acceleration;
		private String title;
		private GriddedRegion gridReg;
		private RSQSimEvent event;
		private RSQSimEventSlipTimeFunc slipTimeFunc;

		public RenderThread(DiscretizedFunc[] seismograms, int dataMod, CPT seisCPT, boolean logZ,
				File outputDir, String prefix, boolean acceleration, String title, GriddedRegion gridReg,
				RSQSimEvent event, RSQSimEventSlipTimeFunc slipTimeFunc, int numFrames, ArrayDeque<Integer> frameDeque)
						throws IOException {
			this.seismograms = seismograms;
			this.dataMod = dataMod;
			this.seisCPT = seisCPT;
			this.logZ = logZ;
			this.outputDir = outputDir;
			this.prefix = prefix;
			this.acceleration = acceleration;
			this.title = title;
			this.gridReg = gridReg;
			this.event = event;
			this.slipTimeFunc = slipTimeFunc;
			this.numFrames = numFrames;
			this.frameDeque = frameDeque;
		}
		
		private int popFrame() {
			synchronized(frameDeque) {
				if (frameDeque.isEmpty())
					return -1;
				return frameDeque.pop();
			}
		}
		
		@Override
		public void run() {
			GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);
			
			Range xRange = new Range(xyz.getMinLon() - 0.5*gridReg.getLonSpacing(), xyz.getMaxLon() + 0.5*gridReg.getLonSpacing());
			Range yRange = new Range(xyz.getMinLat() - 0.5*gridReg.getLatSpacing(), xyz.getMaxLat() + 0.5*gridReg.getLatSpacing());
			
			TickUnits tus = new TickUnits();
			TickUnit tu;
			if (xRange.getLength() > 3)
				tu = new NumberTickUnit(1d);
			else if (xRange.getLength() > 1.5)
				tu = new NumberTickUnit(0.5);
			else if (xRange.getLength() > 0.75)
				tu = new NumberTickUnit(0.25);
			else
				tu = new NumberTickUnit(0.1);
			tus.add(tu);
			
			String zAxisLabel;
			if (acceleration)
				zAxisLabel = "Horizontal Acceleration (cm/s/s)";
			else
				zAxisLabel = "Horizontal Velocity (cm/s)";
			
			if (logZ)
				zAxisLabel = "Log10 "+zAxisLabel;
			
			XYZPlotSpec spec = new XYZPlotSpec(xyz, seisCPT, title, "Longitude", "Latitude", zAxisLabel);
			XYTextAnnotation timeAnnotation = new XYTextAnnotation("", xRange.getLowerBound()+0.9*xRange.getLength(),
					yRange.getLowerBound()+0.9*yRange.getLength());
			timeAnnotation.setTextAnchor(TextAnchor.TOP_RIGHT);
			timeAnnotation.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 20));
			List<XYTextAnnotation> anns = new ArrayList<>();
			anns.add(timeAnnotation);
			spec.setPlotAnnotations(anns);
			
			PlotPreferences plotPrefs = PlotPreferences.getDefault();
			plotPrefs.setTickLabelFontSize(18);
			plotPrefs.setAxisLabelFontSize(20);
			plotPrefs.setPlotLabelFontSize(21);
			plotPrefs.setBackgroundColor(Color.WHITE);
			XYZGraphPanel xyzGP = new XYZGraphPanel();
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			if (CA_OUTLINE_COLOR != null) {
				DefaultXY_DataSet[] outlines;
				try {
					outlines = RupturePlotGenerator.loadCAOutlines();
				} catch (IOException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
				PlotCurveCharacterstics outlineChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, CA_OUTLINE_COLOR);
				
				for (DefaultXY_DataSet outline : outlines) {
					funcs.add(outline);
					chars.add(outlineChar);
				}
			}
			
			spec.setXYElems(funcs);
			spec.setXYChars(chars);
			
			PlotLineType lineType = PlotLineType.SOLID;
			float lineThickness = 1f;
			
			Map<SimulatorElement, PlotCurveCharacterstics> elemCharMap = new HashMap<>();
			List<SimulatorElement> elems = event.getAllElements();
			
			double maxSlip = slipTimeFunc.getMaxCumulativeSlip();
			CPT slipCPT = new CPT(0d, maxSlip, new Color(200, 255, 200, 30), new Color(100, 255, 100, 80),
					new Color(0, 255, 0, 110), new Color(0, 150, 0, 140));
			PaintScaleLegend slipCPTbar = XYZGraphPanel.getLegendForCPT(slipCPT, "Slip (m)",
					plotPrefs.getAxisLabelFontSize(), plotPrefs.getTickLabelFontSize(), 1d, RectangleEdge.BOTTOM);
			
			int maxDigits = ((numFrames-1)+"").length();
			
			int frameNum;
			
			while ((frameNum = popFrame()) >= 0) {
				System.out.println("Frame "+frameNum);
				
				int i = frameNum * dataMod;
				
				double time = seismograms[0].getX(i);
				timeAnnotation.setText(timeDF.format(time)+" s");
				
				// update XYZ
				for (int l=0; l<seismograms.length; l++) {
					double val = seismograms[l].getY(i);
					if (logZ)
						val = Math.log10(val);
					xyz.set(l, val);
				}
				
				// update source
				for (SimulatorElement elem : elems) {
					double slip = slipTimeFunc.getCumulativeEventSlip(elem.getID(), time);
					if (slip > 0) {
						PlotCurveCharacterstics pChar = elemCharMap.get(elem);
						if (pChar == null) {
							Vertex[] verts = elem.getVertices();
							double[] xs = new double[verts.length+1];
							double[] ys = new double[verts.length+1];
							for (int v=0; v<xs.length; v++) {
								xs[v] = verts[v % verts.length].getLongitude();
								ys[v] = verts[v % verts.length].getLatitude();
							}
							funcs.add(new LightFixedXFunc(xs, ys));
							pChar = new PlotCurveCharacterstics(lineType, lineThickness, Color.WHITE);
							chars.add(pChar);
							elemCharMap.put(elem, pChar);
						}
						pChar.setColor(slipCPT.getColor((float)slip));
					}
				}
				
				String frameStr = (frameNum++)+"";
				while (frameStr.length() < maxDigits)
					frameStr = "0"+frameStr;
				
				String frameName = prefix+"_"+frameStr+".png";
				
				xyzGP.drawPlot(spec, false, false, xRange, yRange);
				xyzGP.getChartPanel().getChart().addSubtitle(slipCPTbar);
				xyzGP.getYAxis().setStandardTickUnits(tus);
				xyzGP.getXAxis().setStandardTickUnits(tus);
				xyzGP.getChartPanel().setSize(PLOT_HEIGHT, PLOT_WIDTH);
				
				try {
					xyzGP.saveAsPNG(new File(outputDir, frameName).getAbsolutePath());
				} catch (IOException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
			}
		}
		
	}
	
	private static Color getWithAlpha(Color c, int alpha) {
		return new Color(c.getRed(), c.getGreen(), c.getBlue(), alpha);
	}
	
	private static DecimalFormat timeDF = new DecimalFormat("0.0#");

	public static void main(String[] args) throws ZipException, IOException, DocumentException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		RSQSimCatalog catalog = Catalogs.BRUCE_2585.instance(baseDir);
		int eventID = 2637969;
		File dir = new File("/data/kevin/bbp/parallel/2018_03_01-rundir2585-event2637969-shakemap-noHF");
//		int eventID = 1670183;
//		File dir = new File("/data/kevin/bbp/parallel/2018_03_01-rundir2585-event1670183-shakemap-noHF");
//		int eventID = 81854;
//		File dir = new File("/data/kevin/bbp/parallel/2018_03_02-rundir2585-event81854-shakemap-noHF");
		File zipFile = new File(dir, "results.zip");
		List<BBP_Site> sites = BBP_Site.readFile(dir);
		BBP_ShakeMapSimZipLoader loader = new BBP_ShakeMapSimZipLoader(zipFile, sites);
		
		File sitesXML = new File(dir, "sites.xml");
		GriddedRegion gridReg = ShakemapPlotter.loadGriddedRegion(sitesXML);
		
		boolean accel = false;
		double maxFPS = 20;
		
		File outputDir = new File("/tmp/shakemove");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		String title = "ShakeMovie #"+eventID;
		String prefix = "shake_movie_"+eventID;
		
		if (accel)
			prefix += "_accel";
		else
			prefix += "_vel";
		
		RSQSimEvent event = catalog.loader().byID(eventID);
		RSQSimEventSlipTimeFunc slipTimeFunc = catalog.getSlipTimeFunc(event);
		
		plotShakeMovie(loader, gridReg, sites, event, slipTimeFunc, title, outputDir, prefix, accel, maxFPS);
	}

}
