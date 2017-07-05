package scratch.kevin.cybershake;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Collections;
import java.util.List;

import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.cybershake.calc.HazardCurveComputation;
import org.opensha.sha.cybershake.calc.RuptureVariationProbabilityModifier;
import org.opensha.sha.cybershake.db.CachedPeakAmplitudesFromDB;
import org.opensha.sha.cybershake.db.CybershakeIM;
import org.opensha.sha.cybershake.db.CybershakeRun;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.db.HazardCurve2DB;
import org.opensha.sha.cybershake.db.MeanUCERF2_ToDB;
import org.opensha.sha.cybershake.db.Runs2DB;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class RVDownsampleModTest implements RuptureVariationProbabilityModifier {
	
	private double fractToKeep;
	private CachedPeakAmplitudesFromDB amps2db;
	
	public RVDownsampleModTest(double fractToKeep, CachedPeakAmplitudesFromDB amps2db) {
		this.fractToKeep = fractToKeep;
		this.amps2db = amps2db;
	}

	@Override
	public List<Double> getVariationProbs(int sourceID, int rupID,
			double originalProb, CybershakeRun run, CybershakeIM im) {
		int num;
		try {
			num = amps2db.getIM_Values(run.getRunID(), sourceID, rupID, im).size();
		} catch (SQLException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		int numToKeep = (int)((double)num*fractToKeep+0.5);
		if (numToKeep < 1)
			numToKeep = 1;
		if (numToKeep == num)
			return null;
		
		double ratePer = originalProb/(double)numToKeep;
		
		List<Integer> indexes = Lists.newArrayList();
		for (int i=0; i<num; i++)
			indexes.add(i);
		Collections.shuffle(indexes);
		
		List<Double> probs = Lists.newArrayList();
		for (int i=0; i<num; i++)
			probs.add(0d);
		
		for (int i=0; i<numToKeep; i++)
			probs.set(indexes.get(i), ratePer);
		
		return probs;
	}
	
	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp/cs_rv_downsample");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		DBAccess db = Cybershake_OpenSHA_DBApplication.getDB();
		ERF erf = MeanUCERF2_ToDB.createUCERF2ERF();
		CachedPeakAmplitudesFromDB amps2db = new CachedPeakAmplitudesFromDB(
				db, new File("/home/kevin/CyberShake/MCER/.amps_cache"), erf);
		HazardCurveComputation calc = new HazardCurveComputation(db);
		calc.setPeakAmpsAccessor(amps2db);
		
		double[] levels = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
		
		List<Integer> runIDs = Lists.newArrayList();
		List<String> siteNames = Lists.newArrayList();
		
		runIDs.add(3870);
		siteNames.add("LADT");
		
		runIDs.add(3970);
		siteNames.add("USC");
		
		runIDs.add(3873);
		siteNames.add("STNI");
		
		runIDs.add(3880);
		siteNames.add("SBSM");
		
		runIDs.add(3878);
		siteNames.add("PAS");
		
		int imTypeID = 21;
		
		Runs2DB runs2db = new Runs2DB(db);
		List<CybershakeRun> runs = Lists.newArrayList();
		for (int runID : runIDs)
			runs.add(runs2db.getRun(runID));
		CybershakeIM imType = new HazardCurve2DB(db).getIMFromID(imTypeID);
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, 1d);
		
		for (int i=0; i<runs.size(); i++) {
			CybershakeRun run = runs.get(i);
			String siteName = siteNames.get(i);
			System.out.println("Now doing curve test for "+run.getRunID()+" ("+siteName+")");
			calc.setPeakAmpsAccessor(new CachedPeakAmplitudesFromDB(
					db, new File("/home/kevin/CyberShake/MCER/.amps_cache"), erf));
			List<Double> xVals = Lists.newArrayList();
			for (Point2D pt : new IMT_Info().getDefaultHazardCurve(SA_Param.NAME))
				xVals.add(pt.getX());
			System.out.println("Calculating Original Curve");
			calc.setRupVarProbModifier(null);
			DiscretizedFunc origCurve = calc.computeHazardCurve(xVals, run, imType);
			origCurve.setName("Original");
			
			List<DiscretizedFunc> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			
			funcs.add(origCurve);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			
			for (double level : levels) {
				calc.setRupVarProbModifier(new RVDownsampleModTest(level, amps2db));
				System.out.println("Calculating Modified Curve");
				DiscretizedFunc modCurve = calc.computeHazardCurve(xVals, run, imType);
				modCurve.setName((int)(level*100d)+"%");
				System.out.println("Done calculating");
				
				funcs.add(0, modCurve);
				chars.add(0, new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, cpt.getColor((float)level)));
			}
			
			PlotSpec spec = new PlotSpec(funcs, chars, siteName+" Comparison", "3sec SA (g)", "1-year POE");
			spec.setLegendVisible(true);
			GraphWindow gw = new GraphWindow(spec);
			gw.setXLog(true);
			gw.setYLog(true);
			gw.setX_AxisRange(1e-2, 3d);
			gw.setY_AxisRange(1e-8, 1e-1);
			gw.saveAsPNG(new File(outputDir, "curve_comparison_"+siteName+".png").getAbsolutePath());
			gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		}
	}

}
