package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipFile;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.iden.CatalogLengthLoadIden;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simCompare.IMT;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.ruptures.CatalogGMPE_Compare.EventComparison;
import scratch.kevin.simulators.ruptures.rotation.RSQSimRotatedRupVariabilityConfig;
import scratch.kevin.simulators.ruptures.rotation.RotatedRupVariabilityConfig.Quantity;

public class ZScoreEvolutionPlot {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_4983_STITCHED.instance();
		double minMag = 6.5d;
		int skipYears = 0;
		double maxYears = 100000;
		double binDeltaYears = 5000d;
		double maxDist = MPJ_BBP_CatalogSim.CUTOFF_DIST_DEFAULT;
		
		File bbpDir = new File("/data/kevin/bbp/parallel/"
				+ "2020_09_01-rundir4983-all-m6.5-skipYears0-noHF-vmLA_BASIN_500-cs500Sites");
		File bbpFile = new File(bbpDir, "results_rotD.zip");
		
		AttenRelRef gmpeRef = AttenRelRef.ASK_2014;
		
		IMT[] imts = { IMT.SA3P0 };
		
		VelocityModel vm = VelocityModel.LA_BASIN_500;
		
		List<BBP_Site> sites = BBP_Site.readFile(bbpDir);
		
		File gmpeCacheDir = new File(bbpDir, "gmpe_cache");
		Preconditions.checkState(gmpeCacheDir.exists() || gmpeCacheDir.mkdir());
		
		CatalogLengthLoadIden lenIden = new CatalogLengthLoadIden(maxYears);
		
		CatalogGMPE_Compare gmpePageGen = new CatalogGMPE_Compare(catalog, new ZipFile(bbpFile),
				sites, minMag, skipYears, gmpeCacheDir, lenIden, vm, maxDist);
		
		List<EventComparison> comps = gmpePageGen.loadCalcComps(gmpeRef, imts);
		
		DefaultXY_DataSet[] zScoreVsTimes = new DefaultXY_DataSet[imts.length];
		for (int i=0; i<imts.length; i++)
			zScoreVsTimes[i] = new DefaultXY_DataSet();
		
		BBP_CatalogSimZipLoader bbpLoader = gmpePageGen.getSimProv();
		
		for (EventComparison comp : comps) {
			RSQSimEvent rupture = comp.getRupture();
			double time = rupture.getTimeInYears();
			for (Site site : comp.getApplicableSites()) {
				for (int i=0; i<imts.length; i++) {
					IMT imt = imts[i];
					double simVal = bbpLoader.getValue(site, rupture, imt, 0);
					// in log space
					simVal = Math.log(simVal);
					
					double gmpeVal = comp.getLogMean(site, imt);
					double z = (simVal - gmpeVal)/comp.getStdDev(site, imt);
					zScoreVsTimes[i].set(time, z);
				}
			}
		}
		System.out.println("Calculated "+zScoreVsTimes[0].size()+" z-scores");
		
		for (int i=0; i<imts.length; i++) {
			System.out.println("Plotting for "+imts[i]);
			DefaultXY_DataSet xy = zScoreVsTimes[i];
			Mean mean = new Mean();
			StandardDeviation stdDev = new StandardDeviation();
			
			HistogramFunction meanFunc = HistogramFunction.getEncompassingHistogram(
					0d, maxYears, binDeltaYears);
			HistogramFunction stdDevFunc = HistogramFunction.getEncompassingHistogram(
					0d, maxYears, binDeltaYears);
			List<List<Double>> binnedVals = new ArrayList<>();
			for (int j=0; j<meanFunc.size(); j++)
				binnedVals.add(new ArrayList<>());
			for (Point2D pt : xy) {
				mean.increment(pt.getY());
				int index = meanFunc.getClosestXIndex(pt.getX());
				binnedVals.get(index).add(pt.getY());
				stdDev.increment(pt.getY());
			}
			double overallMean = mean.getResult();
			DefaultXY_DataSet meanLine = new DefaultXY_DataSet();
			meanLine.set(0d, overallMean);
			meanLine.set(maxYears, overallMean);
			double overallSD = stdDev.getResult();
			DefaultXY_DataSet stdDevLine = new DefaultXY_DataSet();
			stdDevLine.set(0d, overallSD);
			stdDevLine.set(maxYears, overallSD);
			
			for (int j=0; j<meanFunc.size(); j++) {
				double[] vals = Doubles.toArray(binnedVals.get(j));
				double binMean = StatUtils.mean(vals);
				double binStdDev = Math.sqrt(StatUtils.variance(vals));
				
				meanFunc.set(j, binMean);
				stdDevFunc.set(j, binStdDev);
			}
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			xy.setName("z-scores");
			funcs.add(xy);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.X, 1f, Color.GRAY));
			
			meanFunc.setName("binned mean z-scores");
			funcs.add(meanFunc);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 10f, Color.BLACK));
			
			meanLine.setName("mean z-score");
			funcs.add(meanLine);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.BLACK));
			
			stdDevFunc.setName("binned sigma-fracts");
			funcs.add(stdDevFunc);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 10f, Color.BLUE));
			
			stdDevLine.setName("mean sigma-fract");
			funcs.add(stdDevLine);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.BLUE));
			
			PlotSpec spec = new PlotSpec(funcs, chars, imts[i].getDisplayName()+" z-score evolution",
					"Time (years)", " ");
			spec.setLegendVisible(true);
			
			Range xRange = new Range(0, maxYears);
			Range yRange = new Range(-3d, 3d);
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setTickLabelFontSize(18);
			gp.setAxisLabelFontSize(24);
			gp.setPlotLabelFontSize(24);
			gp.setLegendFontSize(24);
			gp.setBackgroundColor(Color.WHITE);
			
			gp.drawGraphPanel(spec, false, false, xRange, yRange);
			
			String prefix = "z_score_evolution_"+imts[i].getPrefix();
			
			File file = new File(bbpDir, prefix);
			gp.getChartPanel().setSize(1200, 800);
			gp.saveAsPNG(file.getAbsolutePath()+".png");
//			gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		}
		System.out.println("DONE");
		
		// unrelated, but look at rotations and see when the earliest event was
		File csRotDir = new File("/home/kevin/Simulators/catalogs/"
				+ "rundir4983_stitched/cybershake_rotation_inputs");
		Map<Integer, RSQSimEvent> eventsMap = new HashMap<>();
		for (RSQSimEvent event : gmpePageGen.getEvents())
			eventsMap.put(event.getID(), event);
		RSQSimSubSectionMapper mapper = catalog.getSubSectMapper();
		for (File csvFile : csRotDir.listFiles()) {
			if (!csvFile.getName().endsWith(".csv"))
				continue;
			System.out.println("Loading "+csvFile.getName());
			RSQSimRotatedRupVariabilityConfig config = RSQSimRotatedRupVariabilityConfig.loadCSV(
					catalog, csvFile);
			RSQSimEvent earliestEvent = null;
			for (int eventID : config.getValues(Integer.class, Quantity.EVENT_ID)) {
				RSQSimEvent event = eventsMap.get(eventID);
				if (event == null)
					// after our maxYears
					continue;
				if (earliestEvent == null || event.getTime() < earliestEvent.getTime())
					earliestEvent = event;
				if (event.getTimeInYears() < 75000) {
					System.out.println("\tEvent "+event.getID()+" is at "+(float)event.getTimeInYears());
					
					List<List<SubSectionMapping>> mappings = mapper.getAllSubSectionMappings(event);
					for (List<SubSectionMapping> parentMaps : mappings) {
						String name = parentMaps.get(0).getSubSect().getParentSectionName();
						System.out.println("\t\tInclues section: "+name);
					}
					
				}
			}
			
			System.out.println("Earliest event is at "+(float)earliestEvent.getTimeInYears());
		}
		
		System.exit(0);
	}

}
