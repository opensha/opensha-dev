package scratch.kevin.pointSources;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.FocalMech;

import com.google.common.base.Preconditions;

import scratch.kevin.pointSources.PointSourceHazardComparison.PointSourceCalcERF;
import scratch.kevin.pointSources.PointSourceHazardComparison.PointSourceType;

public class InvCDF_RJBTableWriter {

	public static void main(String[] args) throws IOException {
		// calculate the distribution of rJBs for a given mag and horizontal dist
		
		Location centerLoc = new Location(0d, 0d);
//		Location centerLoc = new Location(41d, 0d); // approx. middle latitude for US
		
		GutenbergRichterMagFreqDist mfd = new GutenbergRichterMagFreqDist(5.05, 31, 0.1);
		mfd.setAllButTotMoRate(mfd.getMinX(), mfd.getMaxX(), 1d, 1d);
		double[] plotHighlightMags = {5.05, 6.05, 7.05, 8.05};
		
		EvenlyDiscretizedFunc distFunc = new EvenlyDiscretizedFunc(0d, 301, 1d);
		double[] plotDists = {10d, 50d, 100d, 200d};
		int numEachDist = 10;
		List<List<Location>> distLocs = new ArrayList<>();
		for (int i=0; i<distFunc.size(); i++) {
			double dist = distFunc.getX(i);
			if (dist == 0d) {
				distLocs.add(List.of(centerLoc));
			} else {
				List<Location> locs = new ArrayList<>(numEachDist);
				double deltaEach = 2d*Math.PI/numEachDist;
				for (int j=0; j<numEachDist; j++)
					locs.add(LocationUtils.location(centerLoc, deltaEach*j, dist));
				distLocs.add(locs);
			}
		}
		
		CPT magCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(mfd.getMinX(), mfd.getMaxX());
		
//		PointSourceType sourceType = PointSourceType.OCT_QUAD;
//		int numPerStochastic = 1000;
//		PointSourceType sourceType = PointSourceType.OCT_QUAD_RAND_DAS_DD;
//		int numPerStochastic = 2000;
//		PointSourceType sourceType = PointSourceType.OCT_QUAD_RAND_CELL;
//		int numPerStochastic = 2000;
		PointSourceType sourceType = PointSourceType.OCT_QUAD_RAND_DAS_DD_CELL;
		int numPerStochastic = 5000;
		
		File mainDir = new File("/data/kevin/markdown/nshm23-misc/point_source_corr/");
		File subDir = new File(mainDir, sourceType.name());
		Preconditions.checkState(subDir.exists() || subDir.mkdir());
		File outputDir = new File(subDir, "rjb_dists");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		Random rand = new Random(123456789l);
		
		FocalMech[] mechs = { FocalMech.STRIKE_SLIP, FocalMech.REVERSE };
		int numDistPts = 100;
		
		DecimalFormat oDF = new DecimalFormat("0.##");
		
		for (FocalMech mech : mechs) {
			System.out.println("Building ERFs for "+mech);
			PointSourceCalcERF erf = new PointSourceCalcERF(sourceType, centerLoc, mfd,
					mech.rake(), mech.dip(), numPerStochastic, rand);
			
			System.out.println("Calculating rJB as a function of center location for "+sourceType+", "+mech);
			EvenlyDiscretizedFunc[][] rJBDists = calcRJBDists(erf, distFunc, distLocs, mfd, numDistPts);
			
			File csvFile = new File(outputDir, mech.name()+"_rJB_cumDists.csv");
			System.out.println("Writing "+csvFile.getAbsolutePath());
			CSVFile<String> csv = new CSVFile<>(true);
			List<String> header = new ArrayList<>(mfd.size()+1);
			header.add("Horizontal Distance");
			header.add("Magnitude");
			for (int i=0; i<numDistPts; i++)
				header.add((float)rJBDists[0][0].getX(i)+"");
			csv.addLine(header);
			for (int d=0; d<distFunc.size(); d++) {
				for (int m=0; m<mfd.size(); m++) {
					List<String> line = new ArrayList<>(header.size());
					line.add((float)distFunc.getX(d)+"");
					line.add((float)mfd.getX(m)+"");
					for (int i=0; i<numDistPts; i++)
						line.add((float)rJBDists[d][m].getY(i)+"");
					csv.addLine(line);
				}
			}
			
			csv.writeToFile(csvFile);
			
			for (double plotDist : plotDists) {
				System.out.println("Plotting distributions at Repi="+(float)plotDist);
				
				List<XY_DataSet> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				int d = distFunc.getClosestXIndex(plotDist);
				
				List<XY_DataSet> highlightFuncs = new ArrayList<>();
				List<PlotCurveCharacterstics> highlightChars = new ArrayList<>();
				
				for (int m=0; m<mfd.size(); m++) {
					double mag = mfd.getX(m);
					Color color = magCPT.getColor((float)mag);
					
					EvenlyDiscretizedFunc dist =  rJBDists[d][m];
					boolean highlight = false;
					for (double testMag : plotHighlightMags)
						if ((float)testMag == (float)mag)
							highlight = true;
					
					if (highlight) {
						dist.setName("M"+oDF.format(mag));
						highlightFuncs.add(dist);
						highlightChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 5f, color));
					} else {
						color = new Color((255+color.getRed())/2, (255+color.getGreen())/2, (255+color.getBlue())/2);
						funcs.add(dist);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, color));
					}
				}
				// put these on top
				funcs.addAll(highlightFuncs);
				chars.addAll(highlightChars);
				
				PlotSpec spec = new PlotSpec(funcs, chars, "rEPI="+oDF.format(plotDist)+" km", "Inverse Cumulative Probability", "rJB (km)");
				spec.setLegendInset(RectangleAnchor.BOTTOM_RIGHT);
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				gp.drawGraphPanel(spec, false, false, new Range(0d, 1d), new Range(0d, plotDist+10d));
				
				PlotUtils.writePlots(outputDir, mech.name()+"_rJB_cumDists_"+oDF.format(plotDist)+"km", gp, 800, 800, true, true, false);
			}
		}
	}
	
	private static EvenlyDiscretizedFunc[][] calcRJBDists(PointSourceCalcERF erf, EvenlyDiscretizedFunc distFunc,
			List<List<Location>> distLocs, IncrementalMagFreqDist mfd, int numDistPts) {
		EvenlyDiscretizedFunc[][] ret = new EvenlyDiscretizedFunc[distFunc.size()][mfd.size()];
		
		System.out.println("Calculating for "+distLocs.size()+" distances");
		List<Integer> doneIndexes = new ArrayList<>();
//		IntStream.range(0, distFunc.size()).forEach(d -> {
		IntStream.range(0, distFunc.size()).parallel().forEach(d -> {
			List<List<Double>> magRJBLists = new ArrayList<>(mfd.size());
			for (int m=0; m<mfd.size(); m++)
				magRJBLists.add(new ArrayList<>());
			
			for (Location loc : distLocs.get(d)) {
				for (ProbEqkSource source : erf) {
					for (ProbEqkRupture rup : source) {
						double mag = rup.getMag();
						int m = mfd.getClosestXIndex(mag);
						RuptureSurface surf = rup.getRuptureSurface();
						magRJBLists.get(m).add(surf.getDistanceJB(loc));
					}
				}
			}
			
			for (int m=0; m<mfd.size(); m++) {
				List<Double> rJBList = magRJBLists.get(m);
				Preconditions.checkState(rJBList.size() > 1);
				
				LightFixedXFunc normCDF = ArbDiscrEmpiricalDistFunc.calcQuickNormCDF(rJBList, null);
//				System.out.println("NormCDF for horzDist="+(float)distFunc.getX(d)+", M"+(float)mfd.getX(m)
//					+": "+rJBList.size()+" rJBs, normCDF.size()="+normCDF.size());
//				System.out.println("\tFirst point: "+normCDF.get(0));
//				System.out.println("\tLast point: "+normCDF.get(normCDF.size()-1));
				
				double minDist = normCDF.getX(0);
				double minProb = normCDF.getY(0);
				double maxDist = normCDF.getX(normCDF.size()-1);
				double maxProb = normCDF.getY(normCDF.size()-1);
				
				EvenlyDiscretizedFunc invCDF = new EvenlyDiscretizedFunc(0d, 1d, numDistPts);
				for (int i=0; i<invCDF.size(); i++) {
					double prob = invCDF.getX(i);
					if (i == 0 || prob <= minProb)
						invCDF.set(i, minDist);
					else if (i == numDistPts-1 || prob >= maxProb)
						invCDF.set(i, maxDist);
					else
						invCDF.set(i, normCDF.getFirstInterpolatedX(prob));
				}
				ret[d][m] = invCDF;
			}
			synchronized (doneIndexes) {
				System.out.print(".");
				doneIndexes.add(d);
				if (doneIndexes.size() % 100 == 0)
					System.out.println(" "+doneIndexes.size());
			}
		});
		
		if (doneIndexes.size() % 100 != 0)
			System.out.println(" "+doneIndexes.size());
		
		return ret;
	}

}
