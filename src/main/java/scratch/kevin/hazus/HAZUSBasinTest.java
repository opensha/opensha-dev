package scratch.kevin.hazus;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.impl.CVM4BasinDepth;
import org.opensha.commons.data.siteData.impl.WillsMap2006;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.sha.calc.hazus.parallel.HazusJobWriter;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;

public class HAZUSBasinTest {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		System.out.println("Loading locs");
		LocationList locs = HazusJobWriter.loadCSV(new File("/home/kevin/OpenSHA/hazus/05grid.csv"));
//		CVM4BasinDepth cvm = new CVM4BasinDepth(SiteData.TYPE_DEPTH_TO_1_0);
		CVM4BasinDepth cvm = new CVM4BasinDepth(SiteData.TYPE_DEPTH_TO_2_5);
		WillsMap2006 wills = new WillsMap2006();
		
		DefaultXY_DataSet scatter = new DefaultXY_DataSet();
		HistogramFunction vs30Hist = new HistogramFunction(180, 30, 20d);
		HistogramFunction basinHist = new HistogramFunction(0, 5*8, 0.2d);
		
		System.out.println("Getting vals");
		ArrayList<Double> vals = cvm.getValues(locs);
		ArrayList<Double> vs30Vals = wills.getValues(locs);
		int numNan = 0;
		int tot = vals.size();
		int cnt = 0;
		for (int i=0; i<locs.size(); i++) {
			Location loc = locs.get(i);
			double val = vals.get(i);
			double vs30 = vs30Vals.get(i);
			vs30Hist.add(vs30Hist.getClosestXIndex(vs30), 1d);
			if (Double.isNaN(val))
				numNan++;
			else {
				scatter.set(vs30, val);
				basinHist.add(basinHist.getClosestXIndex(val), 1);
			}
			if (val > DepthTo2pt5kmPerSecParam.MAX) {
				cnt++;
				System.out.println(loc.getLatitude() + ", " + loc.getLongitude() + ": " + val);
			}
		}
		System.out.println("Num above: " + cnt);
		
		GraphWindow gw = new GraphWindow(scatter, "Vs30 vs Z2.5",
				new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLACK));
		gw.setX_AxisLabel("Vs30 (m/s)");
		gw.setY_AxisLabel("Z2.5 (km/s)");
		
		gw = new GraphWindow(vs30Hist, "Vs30 Histogram",
				new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 2f, Color.BLACK));
		gw.setX_AxisLabel("Vs30 (m/s)");
		gw.setY_AxisLabel("Number");
		
		gw = new GraphWindow(basinHist, "Z2.5 Histogram ("+numNan+"/"+tot+" are NaN)",
				new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 2f, Color.BLACK));
		gw.setX_AxisLabel("Z2.5 (km/s)");
		gw.setY_AxisLabel("Number");
	}

}
