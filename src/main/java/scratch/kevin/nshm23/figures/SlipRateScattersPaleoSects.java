package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoProbabilityModel;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.UncertainDataConstraint.SectMappedUncertainDataConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.PaleoseismicConstraintData;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionSlipRates;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

public class SlipRateScattersPaleoSects {

	public static void main(String[] args) throws IOException {
		File mainDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		
		File solFile = new File(mainDir, "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip");
		
		File evenFile = new File(mainDir, "node_branch_averaged/PaleoUncert_EvenFitPaleo.zip");
		File underFile = new File(mainDir, "node_branch_averaged/PaleoUncert_UnderFitPaleo.zip");
		File overFile = new File(mainDir, "node_branch_averaged/PaleoUncert_OverFitPaleo.zip");
		
		File outputDir = new File(mainDir, "misc_plots/slip_rate_scatters_paleo_sects");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		List<String> prefixes = new ArrayList<>();
		List<String> names = new ArrayList<>();
		List<FaultSystemSolution> sols = new ArrayList<>();
		List<Color> zScatterColors = new ArrayList<>();
		
		prefixes.add("full_ba");
		names.add("Full Branch Average");
		sols.add(FaultSystemSolution.load(solFile));
		zScatterColors.add(Color.BLACK);
		
		prefixes.add("even_fit");
		names.add("Even-Fit Paleo Branch");
		sols.add(FaultSystemSolution.load(evenFile));
		PaleoseismicConstraintData evenPaleoData = sols.get(sols.size()-1).getRupSet().requireModule(PaleoseismicConstraintData.class);
		zScatterColors.add(Color.BLUE.darker());
		
		prefixes.add("under_fit");
		names.add("Under-Fit Paleo Branch");
		sols.add(FaultSystemSolution.load(underFile));
		zScatterColors.add(Color.RED.darker());
		
		prefixes.add("over_fit");
		names.add("Over-Fit Paleo Branch");
		sols.add(FaultSystemSolution.load(overFile));
		zScatterColors.add(Color.GREEN.darker());
		
//		int tieTogetherID = 4262; // SJC Mystic Lake
//		int tieTogetherID = 1313; // Elsinore Glen Ivy
//		int tieTogetherID = 5281; // Wasatch Spring Lake
		int tieTogetherID = 2779; // SAF Thousand Palms
		Point2D tieTogetherEven = null;
		Point2D tieTogetherOver = null;
		Point2D tieTogetherUnder = null;
		
		List<DefaultXY_DataSet> zScoreScatters = new ArrayList<>();
		NSHM23_DeformationModels.HARDCODED_FRACTIONAL_STD_DEV = Double.NaN;
		NSHM23_DeformationModels.HARDCODED_FRACTIONAL_STD_DEV_UPPER_BOUND = Double.NaN;
		List<? extends FaultSection> sectsFullSD = NSHM23_DeformationModels.AVERAGE.build(NSHM23_FaultModels.NSHM23_v2);
		
		for (int i=0; i<prefixes.size(); i++) {
			String prefix = prefixes.get(i);
			FaultSystemSolution sol = sols.get(i);
			FaultSystemRupSet rupSet = sol.getRupSet();
			String name = names.get(i);
			
			CSVFile<String> csv = new CSVFile<>(false);
			csv.addLine("Section Index", "Section Name", "Target Slip Rate (mm/yr)", "Solution Slip Rate (mm/yr)", "Parent Section Mapped Paleo Site(s)");
			
			SectSlipRates targets = rupSet.getSectSlipRates();
			SolutionSlipRates solRates = sol.requireModule(SolutionSlipRates.class);
			
			DefaultXY_DataSet scatter = new DefaultXY_DataSet();
			DefaultXY_DataSet scatterNoPaleo = new DefaultXY_DataSet();
			DefaultXY_DataSet scatterWithPaleo = new DefaultXY_DataSet();
			
			PaleoseismicConstraintData paleoData = rupSet.requireModule(PaleoseismicConstraintData.class);
			Map<Integer, List<SectMappedUncertainDataConstraint>> paleoParents = new HashMap<>();
			for (SectMappedUncertainDataConstraint constr : paleoData.getPaleoRateConstraints()) {
				int parentID = rupSet.getFaultSectionData(constr.sectionIndex).getParentSectionId();
				List<SectMappedUncertainDataConstraint> sectConstraints = paleoParents.get(parentID);
				if (sectConstraints == null) {
					sectConstraints = new ArrayList<>();
					paleoParents.put(parentID, sectConstraints);
				}
				sectConstraints.add(constr);
			}
			
			for (int s=0; s<rupSet.getNumSections(); s++) {
				double target = targets.getSlipRate(s)*1e3;
				double solSlip = solRates.get(s)*1e3;
				
				scatter.set(target, solSlip);
				
				List<SectMappedUncertainDataConstraint> sectConstraints = paleoParents.get(rupSet.getFaultSectionData(s).getParentSectionId());
				if (sectConstraints == null)
					scatterNoPaleo.set(target, solSlip);
				else
					scatterWithPaleo.set(target, solSlip);
				
				List<String> line = new ArrayList<>();
				line.add(s+"");
				line.add(rupSet.getFaultSectionData(s).getSectionName());
				line.add((float)target+"");
				line.add((float)solSlip+"");
				if (sectConstraints != null)
					for (SectMappedUncertainDataConstraint constr : sectConstraints)
						line.add(constr.name);
				
				csv.addLine(line);
			}
			
			csv.writeToFile(new File(outputDir, prefix+".csv"));
			
			for (boolean log : new boolean[] {false,true}) {
				Range range = log ? new Range(1e-3, 1e2) : new Range(0, 32);
				for (int j=0; j<3; j++) {
					DefaultXY_DataSet myScatter;
					String myPrefix, title;
					if (j == 0) {
						myScatter = scatter;
						myPrefix = prefix+"_all_sects";
						title = name+", All Subsections";
					} else if (j == 1) {
						myScatter = scatterNoPaleo;
						myPrefix = prefix+"_no_paleo";
						title = name+", Subsections Without Paleo Data";
					} else {
						myScatter = scatterWithPaleo;
						myPrefix = prefix+"_with_paleo";
						title = name+", Subsections With Paleo Data";
					}
					
					if (log)
						myPrefix += "_log";
					
					List<XY_DataSet> funcs = new ArrayList<>();
					List<PlotCurveCharacterstics> chars = new ArrayList<>();
					
					DefaultXY_DataSet oneToOne = new DefaultXY_DataSet();
					oneToOne.set(range.getLowerBound(), range.getLowerBound());
					oneToOne.set(range.getUpperBound(), range.getUpperBound());
					
					funcs.add(oneToOne);
					chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.GRAY));
					
					funcs.add(myScatter);
					chars.add(new PlotCurveCharacterstics(PlotSymbol.X, 3f, Color.BLACK));
					
					PlotSpec spec = new PlotSpec(funcs, chars, title, "Target Slip Rate (mm/yr)", "Solution Slip Rate (mm/yr)");
					
					HeadlessGraphPanel gp = PlotUtils.initHeadless();
					
					gp.drawGraphPanel(spec, log, log, range, range);
					
					PlotUtils.writePlots(outputDir, myPrefix, gp, 850, false, true, true, false);
				}
			}
			
			CSVFile<String> zCSV = new CSVFile<>(true);
			zCSV.addLine("Section Index", "Section Name",
					"Target Slip Rate (mm/yr)", "Slip Rate Std Dev (mm/yr)", "Solution Slip Rate (mm/yr)", "Slip Rate z-score",
					"Paleo Site Name", "Target Paleo Rate", "Paleo Rate Std Dev", "Solution Paleo Rate", "Paleo Rate z-score");
			
			// now sect only z-score comparisons
			DefaultXY_DataSet zScatter = new DefaultXY_DataSet();
			PaleoProbabilityModel paleoProbModel = evenPaleoData.getPaleoProbModel();
			
			DefaultXY_DataSet slip68Scatter = new DefaultXY_DataSet();
			DefaultXY_DataSet slip95Scatter = new DefaultXY_DataSet();
			DefaultXY_DataSet slipOutsideScatter = new DefaultXY_DataSet();
			for (SectMappedUncertainDataConstraint constr : evenPaleoData.getPaleoRateConstraints()) {
				int sectIndex = constr.sectionIndex;
				
				FaultSection sect = sectsFullSD.get(sectIndex);
				double fractSD = sect.getOrigSlipRateStdDev()/sect.getOrigAveSlipRate();
				Preconditions.checkState(Double.isFinite(fractSD));
				
				double target = targets.getSlipRate(sectIndex)*1e3;
				double solSlip = solRates.get(sectIndex)*1e3;
				double sd = target*fractSD;
				
				double slipZ = (solSlip - target)/sd;
				
				double paleoRate = sol.calcTotPaleoVisibleRateForSect(sectIndex, paleoProbModel);
				double paleoSD = constr.getPreferredStdDev();
				double paleoZ = (paleoRate - constr.bestEstimate)/paleoSD;
				
				zScatter.set(slipZ, paleoZ);
				
				zCSV.addLine(sectIndex+"", sect.getSectionName(),
						(float)target+"", (float)sd+"", (float)solSlip+"", (float)slipZ+"",
						constr.name, (float)constr.bestEstimate+"", (float)paleoSD+"", (float)paleoRate+"", (float)paleoZ+"");
				
				if (Math.abs(slipZ) < 1d)
					slip68Scatter.set(target, solSlip);
				else if (Math.abs(slipZ) < 2d)
					slip95Scatter.set(target, solSlip);
				else
					slipOutsideScatter.set(target, solSlip);
				
				if (sectIndex == tieTogetherID) {
					Point2D pt = new Point2D.Double(slipZ, paleoZ);
					if (prefix.contains("even"))
						tieTogetherEven = pt;
					else if (prefix.contains("over"))
						tieTogetherOver = pt;
					else if (prefix.contains("under"))
						tieTogetherUnder = pt;
				}
			}
			
			zScoreScatters.add(zScatter);
			
			double minZ = Math.min(zScatter.getMinX(), zScatter.getMinY());
			double maxZ = Math.min(zScatter.getMaxX(), zScatter.getMaxY());
			minZ = 0.5*Math.floor(minZ*2d);
			maxZ = 0.5*Math.ceil(maxZ*2d);
			Range range = new Range(minZ, maxZ);
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			funcs.add(zScatter);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 3f, Color.BLACK));
			
			PlotSpec spec = new PlotSpec(funcs, chars, name, "Slip Rate Fit z-score", "Paleoseismic Rate Fit z-score");
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(spec, false, false, range, range);
			
			PlotUtils.writePlots(outputDir, prefix+"_z_scores", gp, 850, false, true, true, false);
			zCSV.writeToFile(new File(outputDir, prefix+"_z_scores.csv"));
			
			for (boolean log : new boolean[] {false,true}) {
				range = log ? new Range(1e-1, 1e2) : new Range(0, 32);
				
				funcs = new ArrayList<>();
				chars = new ArrayList<>();
				
				if (slip68Scatter.size() > 0) {
					slip68Scatter.setName("Inside 68% Conf");
					funcs.add(slip68Scatter);
					chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 7f, Color.GREEN.darker()));
				}
				
				if (slip95Scatter.size() > 0) {
					slip95Scatter.setName("Inside 95% Conf");
					funcs.add(slip95Scatter);
					chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 7f, Color.BLUE.darker()));
				}

				if (slipOutsideScatter.size() > 0) {
					slipOutsideScatter.setName("Outside 95% Conf");
					funcs.add(slipOutsideScatter);
					chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 7f, Color.RED.darker()));
				}
				
				DefaultXY_DataSet oneToOne = new DefaultXY_DataSet();
				oneToOne.set(range.getLowerBound(), range.getLowerBound());
				oneToOne.set(range.getUpperBound(), range.getUpperBound());
				
				funcs.add(0, oneToOne);
				chars.add(0, new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
				
				spec = new PlotSpec(funcs, chars, name, "Target Slip Rate (mm/yr)", "Solution Slip Rate (mm/yr)");
				
				gp.drawGraphPanel(spec, log, log, range, range);
				
				PlotUtils.writePlots(outputDir, prefix+"_with_paleo_slip_z_colored"+(log ? "_log" : ""),
						gp, 800, false, true, true, false);
			}
		}
		
		// now consolidated z scatter plot
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		Range range = new Range(-4d, 4d);
		
		if (tieTogetherEven != null && tieTogetherOver != null && tieTogetherUnder != null) {
			tieTogetherEven = ptWithinRange(tieTogetherEven, range);
			tieTogetherOver = ptWithinRange(tieTogetherOver, range);
			tieTogetherUnder = ptWithinRange(tieTogetherUnder, range);
			
			DefaultXY_DataSet xy = new DefaultXY_DataSet();
			xy.set(tieTogetherEven);
			xy.set(tieTogetherOver);
			
			PlotCurveCharacterstics connectChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.LIGHT_GRAY);
//			PlotCurveCharacterstics connectChar = new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, Color.GRAY);
			
			funcs.add(xy);
			chars.add(connectChar);
			
			xy = new DefaultXY_DataSet();
			xy.set(tieTogetherEven);
			xy.set(tieTogetherUnder);
			
			funcs.add(xy);
			chars.add(connectChar);
		}
		
		for (int i=0; i<names.size(); i++) {
			DefaultXY_DataSet scatter = zScoreScatters.get(i);
			
			if (i == 0)
				// skip full BA for plot
				continue;
			
			int numXwithin1 = 0;
			int numYwithin1 = 0;
			int numXYwithin1 = 0;
			
			DefaultXY_DataSet boundedScatter = new DefaultXY_DataSet();
			boundedScatter.setName(names.get(i));
			for (Point2D pt : scatter) {
				boolean xWithin1 = Math.abs(pt.getX()) < 1d;
				boolean yWithin1 = Math.abs(pt.getY()) < 1d;
				if (xWithin1)
					numXwithin1++;
				if (yWithin1)
					numYwithin1++;
				if (xWithin1 && yWithin1)
					numXYwithin1++;
				
				boundedScatter.set(ptWithinRange(pt, range));
			}
			
			System.out.println(names.get(i));
			System.out.println("Slip rates with |z|<1: "+numXwithin1+"/"+scatter.size());
			System.out.println("Paleo rates with |z|<1: "+numYwithin1+"/"+scatter.size());
			System.out.println("Slip & Paleo rates with |z|<1: "+numXYwithin1+"/"+scatter.size());
			
			Color color = zScatterColors.get(i);
			
			funcs.add(boundedScatter);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 5f, color));
			
			System.out.println("Regressing for best fit line");
			
			//create points in a double array
			double[] xs = new double[scatter.size()];
			double[] ys = new double[scatter.size()];
			for (int j=0; j<xs.length; j++) {
				xs[j] = scatter.getX(j);
				ys[j] = scatter.getY(j);
			}
			double avgX = StatUtils.mean(xs);
			double avgY = StatUtils.mean(ys);
			double[][] pointsArray = new double[scatter.size()][2];
			for (int j=0; j<scatter.size(); j++) {
				pointsArray[j][0] = xs[j]-avgX;
				pointsArray[j][1] = ys[j]-avgY;
			}

			//create real matrix
			RealMatrix realMatrix = MatrixUtils.createRealMatrix(pointsArray);

			//create covariance matrix of points, then find eigenvectors
			//see https://stats.stackexchange.com/questions/2691/making-sense-of-principal-component-analysis-eigenvectors-eigenvalues

			Covariance covariance = new Covariance(realMatrix);
			RealMatrix covarianceMatrix = covariance.getCovarianceMatrix();
			EigenDecomposition ed = new EigenDecomposition(covarianceMatrix);
			
			System.out.println("Center: "+(float)avgX+", "+(float)avgY);
			
			System.out.print("Eigenvalues:");
			for (double val : ed.getRealEigenvalues())
				System.out.print("\t"+(float)val);
			System.out.println();
			RealVector vect1 = ed.getEigenvector(0);
			System.out.println("Eigenvector 1: "+vect1);
			
			double vectX = vect1.getEntry(0);
			double vectY = vect1.getEntry(1);
			double slope = vectY/vectX;
			
			EvenlyDiscretizedFunc centeredFit = new EvenlyDiscretizedFunc(range.getLowerBound()-1d, range.getUpperBound()+1d, 100);
			// use one to one for x values
			for (int j=0; j<centeredFit.size(); j++) {
				double x = centeredFit.getX(j);
				double y = slope*x;
				centeredFit.set(j, y);
			}
			ArbitrarilyDiscretizedFunc fit = new ArbitrarilyDiscretizedFunc();
			for (Point2D pt : centeredFit)
				fit.set(pt.getX()+avgX, pt.getY()+avgY);
			
			funcs.add(0, fit);
			chars.add(0, new PlotCurveCharacterstics(PlotLineType.DASHED, 3f,
					new Color(color.getRed(), color.getGreen(), color.getBlue(), 120)));
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, " ", "Slip Rate Fit z-score", "Paleoseismic Rate Fit z-score");
		spec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.getPlotPrefs().scaleFontSizes(1.15d);
		
//		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		gp.drawGraphPanel(spec, false, false, range, range);
		
		PlotUtils.setXTick(gp, 1d);
		PlotUtils.setYTick(gp, 1d);
		
		PlotUtils.writePlots(outputDir, "combined_z_scores", gp, 850, false, true, true, false); 
	}
	
	private static Point2D ptWithinRange(Point2D pt, Range range) {
		double x = Math.max(pt.getX(), range.getLowerBound());
		x = Math.min(x, range.getUpperBound());
		double y = Math.max(pt.getY(), range.getLowerBound());
		y = Math.min(y, range.getUpperBound());
		return new Point2D.Double(x, y);
	}
	
//	private static class BruteForceBestFitLine {
//		
//		private XY_DataSet scatter;
//		private EvenlyDiscretizedFunc xFunc;
//		private EvenlyDiscretizedFunc yFunc;
//		
//		private double bestSlope = Double.NaN;
//		private double bestIntercept = Double.NaN;
//		private double bestLSQ = Double.POSITIVE_INFINITY;
//
//		public BruteForceBestFitLine(XY_DataSet scatter, Range xRange, Range yRange, int numX, int numY) {
//			this.scatter = scatter;
//			this.xFunc = new EvenlyDiscretizedFunc(xRange.getLowerBound(), xRange.getUpperBound(), numX);
//			this.yFunc = new EvenlyDiscretizedFunc(yRange.getLowerBound(), yRange.getUpperBound(), numY);
//		}
//		
//		public void calc() {
//			List<Point2D> allEdgePoints = new ArrayList<>();
//			
//			double minX = xFunc.getMinX();
//			double maxX = xFunc.getMaxX();
//			double minY = yFunc.getMinX();
//			double maxY = yFunc.getMaxX();
//			
//			// bottom and top edges
//			for (int i=0; i<xFunc.size(); i++) {
//				allEdgePoints.add(new Point2D.Double(xFunc.getX(i), minY));
//				allEdgePoints.add(new Point2D.Double(xFunc.getX(i), maxY));
//			}
//			
//			// left and right edges
//			for (int i=0; i<yFunc.size(); i++) {
//				allEdgePoints.add(new Point2D.Double(minX, yFunc.getX(i)));
//				allEdgePoints.add(new Point2D.Double(maxX, yFunc.getX(i)));
//			}
//			
//			List<CompletableFuture<Void>> futures = new ArrayList<>(allEdgePoints.size());
//			
//			for (Point2D pt : allEdgePoints) {
//				futures.add(CompletableFuture.runAsync(new Runnable() {
//					
//					@Override
//					public void run() {
//						// TODO Auto-generated method stub
//						
//					}
//				}))
//			}
//		}
//		
//		private class PointCalcRunnable implements Runnable {
//			
//			private Point2D point1;
//			private List<Point2D> allEdgePoints;
//
//			PointCalcRunnable(Point2D point1, List<Point2D> allEdgePoints) {
//				this.point1 = point1;
//				this.allEdgePoints = allEdgePoints;
//			}
//
//			@Override
//			public void run() {
//				for (Point2D point2 : allEdgePoints) {
//					// see if it's on the same edge
//					if (point2.getX() == point1.getX()) {
//						if (point)
//					}
//				}
//			}
//			
//		}
//	}

}
