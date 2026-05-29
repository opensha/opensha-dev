package scratch.kevin;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.comcat.ComcatAccessor;
import org.opensha.commons.data.comcat.ComcatRegion;
import org.opensha.commons.data.comcat.ComcatRegionAdapter;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.AnalysisRegions;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import net.mahdilamb.colormap.Colors;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;

public class ComcatMFDPlot {

	public static void main(String[] args) throws IOException {
		boolean doLinear = false;
		
//		double minMag = 4d;
//		double maxMag = 8d;
//		int startYear = 1900;
////		int startYear = 1981;
//		int endYear = 2025;
		
//		double minMag = 2d;
//		double maxMag = 8d;
//		int startYear = 2026;
//		int endYear = 2027;
		
		double minMag = 3d;
		double maxMag = 8d;
		int startYear = 2000;
		int endYear = 2027;
		doLinear = true;
		
//		double minMag = 3d;
//		double maxMag = 8d;
//		int startYear = 1980;
//		int endYear = 2027;
//		doLinear = true;
		
		Date startDate = new GregorianCalendar(startYear, 0, 1).getTime();
		Date endDate = new GregorianCalendar(endYear, 0, 1).getTime();
		
		Region region = AnalysisRegions.CONUS_U3_RELM.load();
		String title = "ComCat Events, CA, "+startYear+"-"+endYear;
		String prefix = "ca_comcat_mfd_plot_"+startYear+"_"+endYear;
		
//		Region region = new CaliforniaRegions.RELM_SOCAL();
//		String title = "ComCat Events, Southern CA, "+startYear+"-"+endYear;
//		String prefix = "socal_comcat_mfd_plot_"+startYear+"_"+endYear;
		
//		Region region = NSHM23_RegionLoader.loadFullConterminousWUS();
//		String title = "ComCat Events, Western US, "+startYear+"-"+endYear;
//		String prefix = "wus_comcat_mfd_plot_"+startYear+"_"+endYear;
		
		File outputDir = new File("/tmp");
		
		double minDepth = -10;
		double maxDepth = 200;
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(minMag+0.01, maxMag-0.01);
		
		ComcatRegion cReg = new ComcatRegionAdapter(region);
		
		ComcatAccessor access = new ComcatAccessor();
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		ObsEqkRupList events = access.fetchEventList(null, startDate.getTime(), endDate.getTime(), minDepth, maxDepth,
				cReg, false, false, minMag);
		
		System.out.println("Loaded "+events.size()+" events");

		IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		mfd.setName("Incremental");
		
		for (ObsEqkRupture rup : events)
			mfd.add(mfd.getClosestXIndex(rup.getMag()), 1d);

		funcs.add(mfd);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 4f, Colors.tab_red));

		EvenlyDiscretizedFunc cmlMFD = mfd.getCumRateDistWithOffset();
		cmlMFD.setName("Cumulative");
		funcs.add(cmlMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		long startTime = startDate.getTime();
		long endTime = Long.min(System.currentTimeMillis(), endDate.getTime());
		long durationMillis = endTime - startTime;
		double durationDays = durationMillis / ProbabilityModelsCalc.MILLISEC_PER_DAY;
		double durationYears = durationMillis / ProbabilityModelsCalc.MILLISEC_PER_YEAR;
		for (int i=0; i<cmlMFD.size(); i++) {
			double count = cmlMFD.getY(i);
			double perDay = count/durationDays;
			double perYear = count/durationYears;
			System.out.println("M>"+(float)cmlMFD.getX(i)+": "+(int)count+" ("+(float)perDay+" /day, "+(float)perYear+" /year)");
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, title, "Magnitude", "Event Count");
		spec.setLegendInset(true);
		
		HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
		
		double maxY = cmlMFD.getY(0);
		
		maxY = Math.pow(10, Math.ceil(Math.log10(maxY)));
		
		gp.drawGraphPanel(spec, false, true, new Range(minMag, maxMag), new Range(0.9d, maxY));
//		gp.getYAxis().setTickLabelsVisible(false);
		
		PlotUtils.writePlots(outputDir, prefix, gp, 800, 750, true, true, false);
		
		if (doLinear) {
			maxY = cmlMFD.getY(0);
			
			if (maxY > 10000)
				maxY = Math.ceil(maxY/5000d)*5000d;
			else if (maxY > 5000)
				maxY = Math.ceil(maxY/1000d)*1000d;
			else if (maxY > 1000)
				maxY = Math.ceil(maxY/500d)*500d;
			else
				maxY = Math.ceil(maxY/100d)*100d;
			
			gp.drawGraphPanel(spec, false, false, new Range(minMag, maxMag), new Range(0.9d, maxY));
//			gp.getYAxis().setTickLabelsVisible(false);
			
			PlotUtils.writePlots(outputDir, prefix+"_linear", gp, 800, 750, true, true, false);
		}
	}

}