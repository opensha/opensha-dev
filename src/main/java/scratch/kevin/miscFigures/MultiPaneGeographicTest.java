package scratch.kevin.miscFigures;

import java.awt.Color;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.ml.neuralnet.SquareNeighbourhood;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartRenderingInfo;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.Range;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;

import com.google.common.base.Preconditions;

public class MultiPaneGeographicTest {
	
	public static void main(String[] args) throws IOException {
		Region reg = new Region(new Location(42, -125), new Location(44, -114));
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(reg);
		
		Location center = new Location(0.5*(reg.getMinLat() + reg.getMaxLat()), 0.5*(reg.getMinLon() + reg.getMaxLon()));
		double minSpan = Math.min(
				LocationUtils.horzDistanceFast(new Location(reg.getMinLat(), reg.getMinLon()),
						new Location(reg.getMinLat(), reg.getMaxLon())),
				LocationUtils.horzDistanceFast(new Location(reg.getMinLat(), reg.getMinLon()),
						new Location(reg.getMaxLat(), reg.getMinLon())));
		
		double[] squareFracts = {0.1, 0.3, 0.5, 0.7};
		
		List<LocationList> lines = new ArrayList<>();
		for (double fract : squareFracts) {
			double width = minSpan * fract;
			double half = 0.5*width;
			double minLat = LocationUtils.location(center, Math.PI, half).lat;
			double maxLat = LocationUtils.location(center, 0d, half).lat;
			double minLon = LocationUtils.location(center, 0.5*Math.PI, half).lon;
			double maxLon = LocationUtils.location(center, 1.5*Math.PI, half).lon;
			
			lines.add(LocationList.of(new Location(minLat, minLon),
					new Location(maxLat, minLon),
					new Location(maxLat, maxLon),
					new Location(minLat, maxLon),
					new Location(minLat, minLon)));
		}
		
		for (LocationList list : lines)
			for (Location loc : list)
				Preconditions.checkState(reg.contains(loc), "Not in region? %s", loc);
		mapMaker.plotLines(lines, Color.black, 1f);
		PlotSpec spec = mapMaker.buildPlot(" ");
		Range xRange = mapMaker.getXRange();
		Range yRange = mapMaker.getYRange();
		
		List<PlotSpec> specs = new ArrayList<>();
		specs.add(spec);
		specs.add(spec);
		
		List<Range> xRanges = List.of(xRange);
		List<Range> yRanges = new ArrayList<>();
		for (int i=0; i<specs.size(); i++)
			yRanges.add(yRange);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(specs, false, false, xRanges, yRanges);
		
		double expectedAR = PlotUtils.calcAspectRatio(xRange, yRange, true);
		System.out.println("Expected top panel AR: "+(float)expectedAR);
		
//		PlotUtils.setSubPlotWeights(gp, 10001, 1000);
		PlotUtils.setSubPlotWeights(gp, 1, 2);
//		PlotUtils.setSubPlotWeights(gp, 2, 1);
		
		int width = 800;
		ChartPanel cp = gp.getChartPanel();
		int height = PlotUtils.calcHeight(cp, width, expectedAR);
		System.out.println("Calculated full plot dimensions: "+width+" x "+height);
		
		ChartRenderingInfo chartInfo = new ChartRenderingInfo();
		cp.getChart().createBufferedImage(width, height, chartInfo);
		Rectangle2D plotArea = chartInfo.getPlotInfo().getDataArea();
		double myWidth = plotArea.getWidth();
		double myHeight = plotArea.getHeight();
		System.out.println("Plot area: "+(float)myWidth+" x "+(float)myHeight+"; AR="+(float)(myWidth/myHeight));
		
		if (specs.size() > 1) {
			CombinedDomainXYPlot combPlot = (CombinedDomainXYPlot)gp.getPlot();
			double gap = combPlot.getGap();
			System.out.println("Gap is "+(float)gap);
			double heightAfterGaps = myHeight - gap*(specs.size()-1);
			int weightSum = 0;
			List<XYPlot> subPlots = PlotUtils.getSubPlots(gp);
			for (XYPlot subPlot : subPlots)
				weightSum += subPlot.getWeight();
			
			for (int i=0; i<subPlots.size(); i++) {
				int weight = subPlots.get(i).getWeight();
				System.out.println("Subplot "+i+" with weight "+weight);
				double subplotHeight = heightAfterGaps * (double)weight / (double)weightSum;
				System.out.println("\t"+(float)myWidth+" x "+(float)subplotHeight+"; AR="+(float)(myWidth/subplotHeight));
			}
		}
		
		PlotUtils.writePlots(new File("/tmp"), "map_ratio_test", gp, 800, height, true, false, false);
	}

}
