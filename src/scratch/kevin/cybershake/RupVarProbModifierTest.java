package scratch.kevin.cybershake;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.cybershake.calc.HazardCurveComputation;
import org.opensha.sha.cybershake.calc.RuptureVariationProbabilityModifier;
import org.opensha.sha.cybershake.db.CachedPeakAmplitudesFromDB;
import org.opensha.sha.cybershake.db.CybershakeIM;
import org.opensha.sha.cybershake.db.CybershakeRun;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.db.HazardCurve2DB;
import org.opensha.sha.cybershake.db.MeanUCERF2_ToDB;
import org.opensha.sha.cybershake.db.PeakAmplitudesFromDB;
import org.opensha.sha.cybershake.db.Runs2DB;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.gcim.ui.infoTools.IMT_Info;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class RupVarProbModifierTest implements RuptureVariationProbabilityModifier {
	
	private static final boolean do_rate_space = false;
	
	private PeakAmplitudesFromDB amps2db;
	private double bundleFactor;
	private int minAmps;
	
	public RupVarProbModifierTest(PeakAmplitudesFromDB amps2db, double bundleFactor, int minAmps) {
		Preconditions.checkArgument(bundleFactor > 0d && bundleFactor <= 1d);
		this.bundleFactor = bundleFactor;
		this.amps2db = amps2db;
		this.minAmps = minAmps;
	}

	@Override
	public List<Double> getVariationProbs(int sourceID,
			int rupID, double originalProb, CybershakeRun run, CybershakeIM im) {
		// get the number of amps
		int numAmps;
		try {
			numAmps = amps2db.getIM_Values(run.getRunID(), sourceID, rupID, im).size();
		} catch (SQLException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		if (numAmps < minAmps)
			return null;
		int index = 0;
		int ampsPerBundle = (int)(bundleFactor * (double)numAmps + 0.5);
		if (ampsPerBundle < 1)
			ampsPerBundle = 1;
		
		if (do_rate_space)
			originalProb = -Math.log(1-originalProb)/1d; // 1 year
		
		double probPerRV = originalProb / (double)numAmps;
		double perturb_scale = originalProb * 0.0000001;
		
		Map<Double, List<Integer>> ret = Maps.newHashMap();
		
		while (index < numAmps) {
			List<Integer> indexes = Lists.newArrayList();
			for (int i=0; i<ampsPerBundle && index<numAmps; i++)
				indexes.add(index++);
			double prob = probPerRV * (double)indexes.size();
			// now randomly perturb
			double perturb = perturb_scale*(Math.random()-0.5);
			prob += perturb;
			
			Preconditions.checkState(!ret.containsKey(prob));
			
			ret.put(prob, indexes);
		}
		
		double totalProb = 0d;
		int runningCount = 0;
		
		for (double prob : ret.keySet()) {
			totalProb += prob;
			runningCount += ret.get(prob).size();
		}
		
//		System.out.println("Num Amps: "+numAmps);
		
		Preconditions.checkState(runningCount == numAmps);
		Preconditions.checkState(DataUtils.getPercentDiff(totalProb, originalProb) < 0.01);
		
		if (do_rate_space) {
			// convert back to probabilities
			Map<Double, List<Integer>> ret2 = Maps.newHashMap();
			for (double rate : ret.keySet()) {
				double prob = 1-Math.exp(-rate*1d); // 1 year
				ret2.put(prob, ret.get(rate));
			}
			ret = ret2;
		}
		
		// convert back to lits of probabilities
		List<Double> rvProbs = Lists.newArrayList();
		for (int i=0; i<numAmps; i++)
			rvProbs.add(0d);
		
		for (Double prob : ret.keySet()) {
			List<Integer> indexes = ret.get(prob);
			double probPer = prob/(double)indexes.size();
			for (int rvIndex : indexes)
				rvProbs.set(rvIndex, probPer);
		}
		
//		System.out.println("Num amps: "+numAmps);
		
		return rvProbs;
	}

	public static void main(String[] args) {
		double bundleFactor = 0.01d;
		int minAmps = 0;
		
		DBAccess db = Cybershake_OpenSHA_DBApplication.getDB();
		ERF erf = MeanUCERF2_ToDB.createUCERF2ERF();
		CachedPeakAmplitudesFromDB amps2db = new CachedPeakAmplitudesFromDB(db,
				new File("/home/kevin/CyberShake/MCER/.amps_cache"), erf);
		
		RupVarProbModifierTest mod = new RupVarProbModifierTest(amps2db, bundleFactor, minAmps);
		HazardCurveComputation calc = new HazardCurveComputation(db);
		calc.setPeakAmpsAccessor(amps2db);
		
		ArrayList<Double> xVals = Lists.newArrayList();
		for (Point2D pt : IMT_Info.getUSGS_SA_Function())
			xVals.add(pt.getX());
		
		Runs2DB runs2db = new Runs2DB(db);
		CybershakeRun run = runs2db.getRun(3852);
		CybershakeIM im = new HazardCurve2DB(db).getIMFromID(21);
//		DiscretizedFunc origCurve = calc.computeDeterministicCurve(xVals, run, 128, 1296, im);
		DiscretizedFunc origSourceCurve = calc.computeHazardCurve(xVals, run, im, Lists.newArrayList(128));
		DiscretizedFunc origTotalCurve = calc.computeHazardCurve(xVals, run, im);
		calc.setRupVarProbModifier(mod);
//		DiscretizedFunc curve = calc.computeDeterministicCurve(xVals, run, 128, 1296, im);
		DiscretizedFunc sourceCurve = calc.computeHazardCurve(xVals, run, im, Lists.newArrayList(128));
		DiscretizedFunc totalCurve = calc.computeHazardCurve(xVals, run, im);
		
		// plot curves
		List<DiscretizedFunc> curves = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		curves.add(origSourceCurve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_SQUARE, 4f, Color.BLACK));
		origSourceCurve.setName("Original");
		curves.add(sourceCurve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 4f, Color.BLUE));
		sourceCurve.setName("Modified");
		
		PlotSpec spec = new PlotSpec(curves, chars, "RV Mod Source Test", "SA", "Exceed. Prob");
		spec.setLegendVisible(true);
		
		GraphWindow gw = new GraphWindow(spec);
		gw.setYLog(true);
		gw.setAxisRange(new Range(0, 0.4), new Range(1e-8, 1e-2));
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		
		curves = Lists.newArrayList();
		
		curves.add(origTotalCurve);
		origTotalCurve.setName("Original");
		curves.add(totalCurve);
		totalCurve.setName("Modified");
		
		spec = new PlotSpec(curves, chars, "RV Mod Total Test", "SA", "Exceed. Prob");
		spec.setLegendVisible(true);
		
		gw = new GraphWindow(spec);
		gw.setYLog(true);
		gw.setAxisRange(new Range(0, 1.2), new Range(1e-7, 1e-1));
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
	}

}
