package scratch.kevin.ucerf3.eal;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.zip.ZipException;

import org.dom4j.DocumentException;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotWindow;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.imr.AttenRelRef;

import com.google.common.collect.Lists;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.FaultSystemIO;

public class FatalityDollarsScatterPlot {

	public static void main(String[] args) throws ZipException, IOException, DocumentException {
		// true mean FSS which includes rupture mapping information. this must be the exact file used to calulate EALs
		File trueMeanSolFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_"
				+ "COMPOUND_SOL_TRUE_HAZARD_MEAN_SOL_WITH_MAPPING.zip");

		// directory which contains EAL data

		// 90%
//		File dollarDataDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_05_28-ucerf3-99percent-wills-smaller");
		// CEA proxy
		File dollarDataDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_05_28-ucerf3-99percent-wills-smaller");
		String dollarLabel = "$ (Billions)";
		final double dollarScale = 1d/1e6; // portfolio units are in thousands (1e3), so convert to billions by dividing by 1e6
		double maxDollars = 230;
		double deltaDollars = 2;
		int nx = (int)(maxDollars/deltaDollars + 0.5);

		// Fatality portfolio
		File fatalityDataDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_05_28-ucerf3-fatality-smaller");
		String fatalityLabel = "Fatalities";
		final double fatalityScale = 1d;
		double maxFatalities = 3500;
		double deltaFatalities = 30;
		int ny = (int)(maxFatalities/deltaFatalities + 0.5);
		
		boolean includeFault = true;
		boolean includeGridded = false;
		boolean plotScatter = false;
		boolean plotXYZ = true;

		// IMR for which EAL data has already been computed
		AttenRelRef attenRelRef = AttenRelRef.BSSA_2014;

		// Fault model of interest
		FaultModels fm = FaultModels.FM3_1;

		// Branch averaged FSS
		FaultSystemSolution baSol = FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
						+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_"
						+ "COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));

		// Compound fault system solution
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
						+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip"));

		UCERF3_BranchAvgLossFetcher fatalityFetch = new UCERF3_BranchAvgLossFetcher(trueMeanSolFile, cfss, fatalityDataDir);
		UCERF3_BranchAvgLossFetcher dollarFetch = new UCERF3_BranchAvgLossFetcher(trueMeanSolFile, cfss, dollarDataDir);
		
		DefaultXY_DataSet xy = new DefaultXY_DataSet();
		
		// fault based
		if (includeFault) {
			DiscretizedFunc[] fatalityFaults = fatalityFetch.getFaultLosses(attenRelRef, fm, true);
			DiscretizedFunc[] dollarFaults = dollarFetch.getFaultLosses(attenRelRef, fm, true);
			for (int i=0; i<fatalityFaults.length; i++)
				xy.set(calcMeanLoss(dollarFaults[i])*dollarScale, calcMeanLoss(fatalityFaults[i])*fatalityScale);
		}
		
		// gridded
		if (includeGridded) {
			DiscretizedFunc[] fatalityGridded = fatalityFetch.getGriddedMagLossDists(attenRelRef, null);
			DiscretizedFunc[] dollarGridded = dollarFetch.getGriddedMagLossDists(attenRelRef, null);
			for (int i=0; i<fatalityGridded.length; i++)
				for (int j=0; j<fatalityGridded[i].size(); j++)
					xy.set(dollarGridded[i].getY(j)*dollarScale, fatalityGridded[i].getY(j)*fatalityScale);
		}
		
		System.out.println("Max Dollars: "+xy.getMaxX());
		System.out.println("Max Fatalities: "+xy.getMaxY());
		
		if (plotScatter) {
			List<PlotElement> elems = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			
			elems.add(xy);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLACK));
			
			PlotSpec spec = new PlotSpec(elems, chars, "Dollars vs Fatalities Per Rupture", dollarLabel, fatalityLabel);
			
			GraphWindow gw = new GraphWindow(spec);
			
			gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		}
		
		if (plotXYZ) {
			EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(nx, ny,
					0.5*deltaDollars, 0.5*deltaFatalities, deltaDollars, deltaFatalities);
			for (Point2D pt : xy) {
				int xInd = xyz.getXIndex(pt.getX());
				int yInd = xyz.getYIndex(pt.getY());
				if (xInd >= 0 && xInd < nx && yInd >= 0 && yInd < ny)
					xyz.set(xInd, yInd, xyz.get(xInd, yInd)+1);
			}
			xyz.log10();
			CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, xyz.getMaxZ());
			cpt.setBelowMinColor(Color.WHITE);
			cpt.setNanColor(Color.WHITE);
			XYZPlotSpec spec = new XYZPlotSpec(xyz, cpt, "Dollars vs Fatalities Per Rupture", dollarLabel, fatalityLabel, "Log10(Number)");
			XYZPlotWindow gw = new XYZPlotWindow(spec);
			gw.setDefaultCloseOperation(XYZPlotWindow.EXIT_ON_CLOSE);
		}
	}
	
	private static final double calcMeanLoss(DiscretizedFunc faultLoss) {
		double loss = 0;
		for (Point2D pt : faultLoss)
			// x=loss, y=weight
			loss += pt.getX()*pt.getY();
		return loss;
	}

}
