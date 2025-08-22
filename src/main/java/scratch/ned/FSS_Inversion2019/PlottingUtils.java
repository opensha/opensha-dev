package scratch.ned.FSS_Inversion2019;

import java.awt.Color;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.data.Range;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;

public class PlottingUtils {
	
	final static double defaultWidthInches = 3.5;
	final static double defaultHeightInches = 3.0;
	
	final static boolean integerYaxisTickLabeIncrements = false;
	
	/**
	 * The general x-y plotting method with default width and height
	 * @param funcs
	 * @param plotChars
	 * @param plotName
	 * @param xAxisLabel
	 * @param yAxisLabel
	 * @param xAxisRange
	 * @param yAxisRange
	 * @param logX
	 * @param logY
	 * @param fileNamePrefix - set a null if you don't want to save to files
	 * @param popupWindow - set as false if you don't want a pop-up windows with the plots
	 */
	public static void writeAndOrPlotFuncs(
			ArrayList<XY_DataSet> funcs, 
			ArrayList<PlotCurveCharacterstics> plotChars, 
			String plotName,
			String xAxisLabel,
			String yAxisLabel,
			Range xAxisRange,
			Range yAxisRange,
			boolean logX,
			boolean logY,
			String fileNamePrefix, 
			boolean popupWindow) {
		
		writeAndOrPlotFuncs(funcs, plotChars, plotName, xAxisLabel, yAxisLabel, xAxisRange, yAxisRange, 
				logX, logY, defaultWidthInches, defaultHeightInches, fileNamePrefix,  popupWindow);
	}
	
	
	/**
	 * Plotting method for a list of UncertainArbDiscFunc with default width and height
	 * @param funcs
	 * @param plotChars
	 * @param plotName
	 * @param xAxisLabel
	 * @param yAxisLabel
	 * @param xAxisRange
	 * @param yAxisRange
	 * @param logX
	 * @param logY
	 * @param fileNamePrefix - set a null if you don't want to save to files
	 * @param popupWindow - set as false if you don't want a pop-up windows with the plots
	 */
	public static void writeAndOrPlotUncertFuncs(
			ArrayList<UncertainArbDiscFunc> uncertFuncs, 
			ArrayList<PlotCurveCharacterstics> plotChars, 
			String plotName,
			String xAxisLabel,
			String yAxisLabel,
			Range xAxisRange,
			Range yAxisRange,
			boolean logX,
			boolean logY,
			String fileNamePrefix, 
			boolean popupWindow) {
		// make the douplicates for shaded regions
		ArrayList<XY_DataSet> funcs = new ArrayList<XY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars2 = new ArrayList<PlotCurveCharacterstics>();
		for(int i=0; i<uncertFuncs.size(); i++) {
			UncertainArbDiscFunc func = uncertFuncs.get(i);
			funcs.add(func);
			funcs.add(func);
			plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, plotChars.get(i).getColor()));
			plotChars2.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, plotChars.get(i).getColor()));
		}
		
		writeAndOrPlotFuncs(funcs, plotChars2, plotName, xAxisLabel, yAxisLabel, xAxisRange, yAxisRange, 
				logX, logY, defaultWidthInches, defaultHeightInches, fileNamePrefix,  popupWindow);
	}

	
	
	/**
	 * The general x-y plotting method with default width and height, for a single function
	 * @param func
	 * @param plotChar
	 * @param plotName
	 * @param xAxisLabel
	 * @param yAxisLabel
	 * @param xAxisRange
	 * @param yAxisRange
	 * @param logX
	 * @param logY
	 * @param fileNamePrefix - set a null if you don't want to save to files
	 * @param popupWindow - set as false if you don't want a pop-up windows with the plots
	 */
	public static void writeAndOrPlotFuncs(
			XY_DataSet func, 
			PlotCurveCharacterstics plotChar, 
			String plotName,
			String xAxisLabel,
			String yAxisLabel,
			Range xAxisRange,
			Range yAxisRange,
			boolean logX,
			boolean logY,
			String fileNamePrefix, 
			boolean popupWindow) {
		
		ArrayList<XY_DataSet> funcs = new ArrayList<XY_DataSet>(); 
		funcs.add(func);
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(plotChar);
		writeAndOrPlotFuncs(funcs, plotChars, plotName, xAxisLabel, yAxisLabel, xAxisRange, yAxisRange, 
				logX, logY, defaultWidthInches, defaultHeightInches, fileNamePrefix,  popupWindow);
	}


	/**
	 * The general x-y plotting method with default width and height, for a single function
	 * @param func
	 * @param plotChar
	 * @param plotName
	 * @param xAxisLabel
	 * @param yAxisLabel
	 * @param xAxisRange
	 * @param yAxisRange
	 * @param logX
	 * @param logY
	 * @param fileNamePrefix - set a null if you don't want to save to files
	 * @param popupWindow - set as false if you don't want a pop-up windows with the plots
	 */
	public static void writeAndOrPlotFuncs(
			XY_DataSet func, 
			PlotCurveCharacterstics plotChar, 
			String plotName,
			String xAxisLabel,
			String yAxisLabel,
			Range xAxisRange,
			Range yAxisRange,
			boolean logX,
			boolean logY,
			double widthInches,
			double heightInches,
			String fileNamePrefix, 
			boolean popupWindow) {
		
		ArrayList<XY_DataSet> funcs = new ArrayList<XY_DataSet>(); 
		funcs.add(func);
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(plotChar);
		writeAndOrPlotFuncs(funcs, plotChars, plotName, xAxisLabel, yAxisLabel, xAxisRange, yAxisRange, 
				logX, logY, widthInches, heightInches, fileNamePrefix,  popupWindow);
	}

	
	/**
	 * The general x-y plotting method
	 * @param funcs
	 * @param plotChars
	 * @param plotName
	 * @param xAxisLabel
	 * @param yAxisLabel
	 * @param xAxisRange
	 * @param yAxisRange
	 * @param logX
	 * @param logY
	 * @param widthInches
	 * @param heightInches
	 * @param fileNamePrefix - set a null if you don't want to save to files
	 * @param popupWindow - set as false if you don't want a pop-up windows with the plots
	 */
	public static PlotSpec writeAndOrPlotFuncs(
			ArrayList<XY_DataSet> funcs, 
			ArrayList<PlotCurveCharacterstics> plotChars, 
			String plotName,
			String xAxisLabel,
			String yAxisLabel,
			Range xAxisRange,
			Range yAxisRange,
			boolean logX,
			boolean logY,
			double widthInches,
			double heightInches,
			String fileNamePrefix, 
			boolean popupWindow) {
		
		if(popupWindow) {
			
			GraphWindow graph = new GraphWindow(funcs, plotName);

			if(xAxisRange != null)
				graph.setX_AxisRange(xAxisRange.getLowerBound(),xAxisRange.getUpperBound());
			if(yAxisRange != null)
				graph.setY_AxisRange(yAxisRange.getLowerBound(),yAxisRange.getUpperBound());
			graph.setXLog(logX);
			graph.setYLog(logY);
			graph.setPlotChars(plotChars);
			graph.setX_AxisLabel(xAxisLabel);
			graph.setY_AxisLabel(yAxisLabel);
			graph.setTickLabelFontSize(18);
			graph.setAxisLabelFontSize(20);

		}
		
		PlotSpec spec = new PlotSpec(funcs, plotChars, plotName, xAxisLabel, yAxisLabel);
		
		if (fileNamePrefix != null){
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setUserBounds(xAxisRange, yAxisRange);
			gp.setTickLabelFontSize(9);
			gp.setAxisLabelFontSize(11);
			gp.setPlotLabelFontSize(9);
			gp.setBackgroundColor(Color.WHITE);
			gp.drawGraphPanel(spec, logX, logY); // spec can be a list
			int width = (int)(widthInches*72.);
			int height = (int)(heightInches*72.);
			gp.getChartPanel().setSize(width, height); 

			if(integerYaxisTickLabeIncrements)
				gp.getChartPanel().getChart().getXYPlot().getRangeAxis().setStandardTickUnits(NumberAxis.createIntegerTickUnits());
			
//			XYTextAnnotation annotation = new XYTextAnnotation("here",xAxisRange.getCentralValue(),yAxisRange.getCentralValue());
//			gp.getChartPanel().getChart().getXYPlot().addAnnotation(annotation);	
//
			// None of these worked for setting the x-axis width:
//			gp.getXAxis().setFixedDimension(width-72);
//			gp.getYAxis().setFixedDimension(height-36);
//			gp.getPlot().getDomainAxis().setFixedDimension(width-72);;
//			gp.getPlot().getRangeAxis().setFixedDimension(height-72);;
			try {
				gp.saveAsPNG(fileNamePrefix+".png");
				gp.saveAsPDF(fileNamePrefix+".pdf");
				gp.saveAsTXT(fileNamePrefix+".txt");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		return spec;
	}
	
	

	
	
	/**
	 * This is the general 2D plotting method; the color scale is for log rates (-1 to -9).
	 * @param xyzData
	 * @param title
	 * @param xAxisLabel
	 * @param yAxisLabel
	 * @param zAxisLabel
	 * @param fileNamePrefix - set as null if no files are to be saved
	 * @param popUpWindow - this tells whether to make a pop-up plot and save it
	 * @param xAddOn - amount to add to min and max data values for plotting range
	 * @param yAddOn - amount to add to min and max data values for plotting range
	 * @param widthInches - for written out plot
	 * @param heightInches - for written out plot

	 */
	protected static XYZPlotSpec make2D_plot(EvenlyDiscrXYZ_DataSet xyzData, String title,
			String xAxisLabel, String yAxisLabel, String zAxisLabel, String fileNamePrefix, 
			boolean popupWindow, Range xRange, Range yRange, Range zRange, double widthInches,
			double heightInches) {
		
//		cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(zRange.getLowerBound(),zRange.getUpperBound());

		CPT cpt = new CPT();;
		Color purple = new Color(128,0,128);
		Color orange = new Color(255,165,0); // this looks better than the default			
		cpt.add(new CPTVal(-9f, Color.BLACK, -8f, purple));
		cpt.add(new CPTVal(-8f, purple, -7f, Color.BLUE));
		cpt.add(new CPTVal(-7f, Color.BLUE, -6f, Color.CYAN));
		cpt.add(new CPTVal(-6f, Color.CYAN, -5f, Color.GREEN));
		cpt.add(new CPTVal(-5f, Color.GREEN, -4f, Color.YELLOW));
		cpt.add(new CPTVal(-4f, Color.YELLOW, -3f, orange));
		cpt.add(new CPTVal(-3f, orange, -2f, Color.RED));
		cpt.add(new CPTVal(-2f, Color.RED, -1f, Color.MAGENTA));
		cpt.setBelowMinColor(Color.BLACK);
		cpt.setAboveMaxColor(Color.MAGENTA);
		cpt.setNanColor(Color.white);

		if(zRange == null) {
			double ratio = Math.pow(10, xyzData.getMaxZ()-xyzData.getMinZ());
//			System.out.println("ratio\t"+ratio);
			if(ratio>1.01) {	// more than 1% difference
				zRange = new Range(xyzData.getMinZ(), xyzData.getMaxZ());
			}
			else {
				zRange = new Range(xyzData.getMinZ()-Math.log10(2d), xyzData.getMinZ()+Math.log10(2d));  // factor of 2 above and below
			}
	
		}
			
		CPT cptAlt = cpt.rescale(zRange.getLowerBound(), zRange.getUpperBound());
			
		XYZPlotSpec spec = new XYZPlotSpec(xyzData, cptAlt, title, xAxisLabel, yAxisLabel, zAxisLabel);
				
//		System.out.println(xyzData.getMinZ()+"\t"+xyzData.getMaxZ());
		
//		if(xyzData.getMaxX() > xyzData.getMinX()+xyzData.getGridSpacingX()/2)
//			xRange = new Range(xyzData.getMinX(),xyzData.getMaxX());
//		else
//			xRange = new Range(xyzData.getMinX(),xyzData.getMaxX()+xyzData.getGridSpacingX());
//		
//		if(xyzData.getMaxY() > xyzData.getMinY()+xyzData.getGridSpacingY()/2)
//			yRange = new Range(xyzData.getMinY(),xyzData.getMaxY());
//		else
//			yRange = new Range(xyzData.getMinY(),xyzData.getMaxY()+xyzData.getGridSpacingY());

		
//System.out.println("xRange\t"+xRange.getLowerBound()+"\t"+xRange.getUpperBound());
//System.out.println("yRange\t"+yRange.getLowerBound()+"\t"+yRange.getUpperBound());
		
		try {
			if(popupWindow) {
				GraphWindow window = new GraphWindow(spec, xRange, yRange);
				window.setVisible(true);
			}
			if (fileNamePrefix != null) {
				int width = (int)(widthInches*72.);
				int height = (int)(heightInches*72.);
				PlotPreferences plotPrefs = PlotPreferences.getDefault();
				plotPrefs.setAxisLabelFontSize(11);
				plotPrefs.setTickLabelFontSize(9);
				plotPrefs.setBackgroundColor(Color.WHITE);
				HeadlessGraphPanel xyzGP = new HeadlessGraphPanel(plotPrefs);
				xyzGP.drawGraphPanel(spec, false, false, xRange, yRange);
				xyzGP.getChartPanel().setSize(width, height);
				xyzGP.saveAsPNG(fileNamePrefix+".png");
				xyzGP.saveAsPDF(fileNamePrefix+".pdf");
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return spec;
	}


}
