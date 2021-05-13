package scratch.kevin;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.comcat.ComcatAccessor;
import org.opensha.commons.data.comcat.ComcatRegion;
import org.opensha.commons.data.comcat.ComcatRegionAdapter;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

public class ComcatBetweenEventsMFDPlot {

	public static void main(String[] args) throws IOException {
		// Ridgecrest
//		String[] eventIDs = {
//				"ci38443095",
//				"ci38443183",
//				"ci38457511"
//		};
//		Region region = new Region(new Location(35.25, -118.25), new Location(36.25, -117));
		
		String[] eventIDs = {
				"ci3019681",
				"ci3031111"
		};
		Region region = new Region(new Location(34.5, -117.2), new Location(33.8, -115.8));
		
		File outputDir = new File("/tmp");
		String prefix = "mfd_plot";
		
		Color[] colors = {
				Color.BLACK,
				Color.BLUE,
				Color.RED
		};
		
		double minDepth = -10;
		double maxDepth = 30;
		
		double minMag = 1d;
		double maxMag = 8d;
		double mfdDelta = 0.1;
		double mfdMinX = minMag + 0.5*mfdDelta;
		int mfdNum = (int)((maxMag - minMag)/mfdDelta)+1;
		
		ComcatRegion cReg = new ComcatRegionAdapter(region);
		
		ComcatAccessor access = new ComcatAccessor();
		
		ObsEqkRupture[] events = new ObsEqkRupture[eventIDs.length];
		
		for (int i=0; i<eventIDs.length; i++)
			events[i] = access.fetchEvent(eventIDs[i], false);
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		List<DiscretizedFunc> cumFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> cumChars = new ArrayList<>();
		
		double maxY = 0;
		
		for (int i=0; i<events.length; i++) {
			long startTime = events[i].getOriginTime()+1000l;
			long endTime = i == events.length-1 ? System.currentTimeMillis() : events[i+1].getOriginTime()-1000l;
			ObsEqkRupList aftershocks = access.fetchEventList(eventIDs[i], startTime, endTime, minDepth, maxDepth,
						cReg, false, false, minMag);
			
			IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(mfdMinX, mfdNum, mfdDelta);
			
			maxY = Math.max(maxY, aftershocks.size());
			for (ObsEqkRupture rup : aftershocks)
				mfd.add(mfd.getClosestXIndex(rup.getMag()), 1d);
			
			mfd.setName("After "+eventIDs[i]+" (M"+(float)events[i].getMag()+")");
			funcs.add(mfd);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 4f, colors[i]));
			
			EvenlyDiscretizedFunc cumMFD = mfd.getCumRateDistWithOffset();
			cumMFD.setName(mfd.getName());
			cumFuncs.add(cumMFD);
			cumChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, colors[i]));
		}
		
		for (boolean cumulative : new boolean[] {false, true}) {
			PlotSpec spec;
			if (cumulative)
				spec = new PlotSpec(cumFuncs, cumChars, "MFDs Between Events", "Magnitude", "Cumulative Count");
			else
				spec = new PlotSpec(funcs, chars, "MFDs Between Events", "Magnitude", "Incremental Count");
			spec.setLegendVisible(true);
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setTickLabelFontSize(18);
			gp.setAxisLabelFontSize(24);
			gp.setPlotLabelFontSize(24);
			gp.setLegendFontSize(28);
			gp.setBackgroundColor(Color.WHITE);
			
			maxY = Math.pow(10, Math.ceil(Math.log10(maxY)));
			
			gp.drawGraphPanel(spec, false, true, new Range(minMag, maxMag), new Range(1d, maxY));
//			gp.getYAxis().setTickLabelsVisible(false);
			
			gp.getChartPanel().setSize(800, 600);
			String myPrefix = prefix;
			if (cumulative)
				myPrefix += "_cumulative";
			File outputFile = new File(outputDir, myPrefix+".png");
			gp.saveAsPNG(outputFile.getAbsolutePath());
		}
		
		
	}

}
