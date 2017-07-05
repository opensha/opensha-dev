package scratch.kevin.cybershake;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.siteData.impl.CVM_CCAi6BasinDepth;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.sha.cybershake.calc.HazardCurveComputation;
import org.opensha.sha.cybershake.calc.RuptureProbabilityModifier;
import org.opensha.sha.cybershake.db.CachedPeakAmplitudesFromDB;
import org.opensha.sha.cybershake.db.CybershakeIM;
import org.opensha.sha.cybershake.db.CybershakeRun;
import org.opensha.sha.cybershake.db.CybershakeSite;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.db.HazardCurve2DB;
import org.opensha.sha.cybershake.db.MeanUCERF2_ToDB;
import org.opensha.sha.cybershake.db.Runs2DB;
import org.opensha.sha.cybershake.db.SiteInfo2DB;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class CCA_WithAndWithoutNSAF_Test implements RuptureProbabilityModifier {

	public static void main(String[] args) throws IOException {
//		int runID = 4675; // PARK, 1-D
//		int runID = 4715; // PARK, 3-D
//		int runID = 4717; // s1207, 3-D
//		int runID = 4719; // s1252, CCA
//		int runID = 4723; // s1252, S4.26
		int runID = 4724; // s1250, S4.26
//		int runID = 4725; // s1250, CCA
		int[] imTypeIDs = {167, 162, 158, 152};
		
		File outputDir = new File("/home/kevin/CyberShake/cca_without_nsaf");
		
		outputDir = new File(outputDir, "run"+runID);
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		DBAccess db = Cybershake_OpenSHA_DBApplication.getDB(Cybershake_OpenSHA_DBApplication.PRODUCTION_HOST_NAME);
		
		HazardCurveComputation calc = new HazardCurveComputation(db);
		MeanUCERF2 erf = (MeanUCERF2) MeanUCERF2_ToDB.createUCERF2ERF();
		calc.setPeakAmpsAccessor(new CachedPeakAmplitudesFromDB(db, null, erf));
		Runs2DB runs2db = new Runs2DB(db);
		HazardCurve2DB curves2db = new HazardCurve2DB(db);
		
		CybershakeRun run = runs2db.getRun(runID);
		CybershakeSite site = new SiteInfo2DB(db).getSiteFromDB(run.getSiteID());
		
		List<Double> xVals = Lists.newArrayList();
		for (Point2D pt : new IMT_Info().getDefaultHazardCurve(SA_Param.NAME))
			xVals.add(pt.getX());
		
		CCA_WithAndWithoutNSAF_Test probMod = new CCA_WithAndWithoutNSAF_Test(erf);
		
		Location loc = site.createLocation();
		double minDist = probMod.calcDistance(loc);
		double z10 = new CVM_CCAi6BasinDepth(CVM_CCAi6BasinDepth.TYPE_DEPTH_TO_1_0).getValue(loc);
		double z25 = new CVM_CCAi6BasinDepth(CVM_CCAi6BasinDepth.TYPE_DEPTH_TO_2_5).getValue(loc);
		System.out.println("Min dist: "+minDist);
		System.out.println("Z1.0: "+z10);
		System.out.println("Z2.5: "+z25);
		
		for (int imTypeID : imTypeIDs) {
			CybershakeIM imType = curves2db.getIMFromID(imTypeID);
			
			calc.setRupProbModifier(null);
			DiscretizedFunc curveWith = calc.computeHazardCurve(xVals, run, imType);
			calc.setRupProbModifier(probMod);
			DiscretizedFunc curveWithout = calc.computeHazardCurve(xVals, run, imType);
			
			curveWith.setName("With N.SAF");
			curveWithout.setName("Without N.SAF");
			
			List<DiscretizedFunc> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			
			funcs.add(curveWith);
			chars.add(new PlotCurveCharacterstics(
					PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 4f, Color.BLACK));
			funcs.add(curveWithout);
			chars.add(new PlotCurveCharacterstics(
					PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 4f, Color.BLUE));
			
			String siteName = site.short_name;
			String xAxisLabel = (int)imType.getVal()+"sec SA (g), "+imType.getComponent().name();
			
			double minX = 1e-2;
			double maxX = 3d;
			
			// add 2% in 50
			DiscretizedFunc uhsLine = new ArbitrarilyDiscretizedFunc();
			uhsLine.setName("2% in 50yr");
			uhsLine.set(minX, 4e-4);
			uhsLine.set(maxX, 4e-4);
			
			funcs.add(uhsLine);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.GRAY));
			
			PlotSpec spec = new PlotSpec(funcs, chars, siteName+" Comparison (Run "+runID+")", xAxisLabel, "1-year POE");
			spec.setLegendVisible(true);
			GraphWindow gw = new GraphWindow(spec);
			gw.getGraphWidget().setBackgroundColor(Color.WHITE);
			gw.setXLog(true);
			gw.setYLog(true);
			gw.setX_AxisRange(minX, maxX);
			gw.setY_AxisRange(1e-8, 1e-1);
			String prefix = siteName+"_run"+runID+"_"+(int)imType.getVal()+"s_"+imType.getComponent().name();
			gw.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
			gw.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
			gw.saveAsTXT(new File(outputDir, prefix+".txt").getAbsolutePath());
			gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		}
		
		db.destroy();
	}
	
	private HashSet<Integer> nsafSources;
	private MeanUCERF2 erf;
	
	private CCA_WithAndWithoutNSAF_Test(MeanUCERF2 erf) {
		this.erf = erf;
//		127	N. San Andreas
//		33	N. San Andreas;SAN
//		34	N. San Andreas;SAN+SAP
//		35	N. San Andreas;SAN+SAP+SAS
//		36	N. San Andreas;SAO
//		37	N. San Andreas;SAO+SAN
//		38	N. San Andreas;SAO+SAN+SAP
//		39	N. San Andreas;SAO+SAN+SAP+SAS
//		40	N. San Andreas;SAP
//		41	N. San Andreas;SAP+SAS
//		42	N. San Andreas;SAS
		nsafSources = new HashSet<Integer>();
		nsafSources.add(127);
		for (int sourceID=33; sourceID<=42; sourceID++)
			nsafSources.add(sourceID);
	}

	@Override
	public double getModifiedProb(int sourceID, int rupID, double origProb) {
		if (nsafSources.contains(sourceID))
			return 0d;
		return origProb;
	}
	
	public double calcDistance(Location loc) {
		// use the full N. SAF char source, source 39
		Site site = new Site(loc);
		return erf.getSource(39).getMinDistance(site);
	}

}
