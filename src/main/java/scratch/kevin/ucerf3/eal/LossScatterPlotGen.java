package scratch.kevin.ucerf3.eal;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.ETAS.ETAS_MultiSimAnalysisTools;
import scratch.UCERF3.griddedSeismicity.GridSourceProvider;
import scratch.UCERF3.utils.FaultSystemIO;

public class LossScatterPlotGen {

	public static void main(String[] args) throws IOException, DocumentException {
		// true mean FSS which includes rupture mapping information. this must be the exact file used to calulate EALs
		File trueMeanSolFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_TRUE_HAZARD_MEAN_SOL_WITH_MAPPING.zip");

		// directory which contains EAL data
		File dataDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_05_28-ucerf3-99percent-wills-smaller");

		// Fault model of interest
		FaultModels fm = FaultModels.FM3_1;

		// Branch averaged FSS
		FaultSystemSolution baSol = FaultSystemIO.loadSol(new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		GridSourceProvider prov = baSol.getGridSourceProvider();

		// Compound fault system solution
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
						+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip"));

		UCERF3_BranchAvgLossFetcher fetch = new UCERF3_BranchAvgLossFetcher(trueMeanSolFile, cfss, dataDir);
		
		Map<AttenRelRef, Double> imrWeightsMap = Maps.newHashMap();
		imrWeightsMap.put(AttenRelRef.CB_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.CY_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.ASK_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.BSSA_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.IDRISS_2014, 0.12);
		
		double[] faultRupLosses = new double[baSol.getRupSet().getNumRuptures()];
		DiscretizedFunc[] griddedLossDists = new DiscretizedFunc[prov.getGriddedRegion().getNodeCount()];
		
		for (AttenRelRef imr : imrWeightsMap.keySet()) {
			double weight = imrWeightsMap.get(imr);
			DiscretizedFunc[] myFaultLosses = fetch.getFaultLosses(imr, fm, true);
			Preconditions.checkState(myFaultLosses.length == faultRupLosses.length);
			
			for (int i=0; i<faultRupLosses.length; i++) {
				double loss = 0;
				double sumWeight = 0;
				for (Point2D pt : myFaultLosses[i]) {
					sumWeight += pt.getY();
					loss += pt.getX()*pt.getY();
				}
				if (loss == 0)
					continue;
				Preconditions.checkState((float)sumWeight == 1f, "Weights don't sum to 1: "+(float)sumWeight);
				faultRupLosses[i] += weight*loss;
			}
			
			DiscretizedFunc[] myGriddedLosses = fetch.getGriddedMagLossDists(imr, prov.getGriddedRegion());
			Preconditions.checkState(griddedLossDists.length == myGriddedLosses.length);
			
			for (int i=0; i<myGriddedLosses.length; i++) {
				if (myGriddedLosses[i] == null)
					continue;
				myGriddedLosses[i].scale(weight);
				if (griddedLossDists[i] == null) {
					griddedLossDists[i] = myGriddedLosses[i];
				} else {
					Preconditions.checkState(griddedLossDists[i].size() == myGriddedLosses[i].size());
					for (int j=0; j<griddedLossDists[i].size(); j++)
						griddedLossDists[i].set(j, griddedLossDists[i].getY(j) + myGriddedLosses[i].getY(j));
				}
			}
		}
		
		String yAxisLabel = "$ (Billions)";
		double thousandsToBillions = 1d/1e6; // portfolio units are in thousands (1e3), so convert to billions by dividing by 1e6
		
		double inflationScalar = 1d/0.9d;
		
		double lossScale = thousandsToBillions*inflationScalar;
		
		DefaultXY_DataSet faultLossVsMag = new DefaultXY_DataSet();
		faultLossVsMag.setName("Fault Ruptures");
		DefaultXY_DataSet griddedLossVsMag = new DefaultXY_DataSet();
		griddedLossVsMag.setName("Gridded Ruptures");
		DefaultXY_DataSet faultLossVsRate = new DefaultXY_DataSet();
		faultLossVsRate.setName("Fault Ruptures");
		DefaultXY_DataSet griddedLossVsRate = new DefaultXY_DataSet();
		griddedLossVsRate.setName("Gridded Ruptures");
		
		FaultSystemRupSet rupSet = baSol.getRupSet();
		for (int i=0; i<faultRupLosses.length; i++) {
			double loss = faultRupLosses[i]*lossScale;
			if (loss == 0)
				continue;
			double mag = rupSet.getMagForRup(i);
			double rate = baSol.getRateForRup(i);
			
			faultLossVsMag.set(mag, loss);
			if (rate > 0)
				faultLossVsRate.set(rate, loss);
		}
		
		for (int i=0; i<griddedLossDists.length; i++) {
			IncrementalMagFreqDist mfd = prov.getNodeMFD(i);
			DiscretizedFunc losses = griddedLossDists[i];
			
			if (losses == null)
				continue;
			
//			Preconditions.checkState(losses.size() == mfd.size(), "Size mismatch! %s != %s", losses.size(), mfd.size());
			
			for (int j=0; j<losses.size(); j++) {
				double loss = losses.getY(j)*lossScale;
				if (loss == 0)
					continue;
				double mag = losses.getX(j);
				double rate = mfd.getY(mag);
				
				griddedLossVsMag.set(mag, loss);
				if (rate > 0)
					griddedLossVsRate.set(rate, loss);
			}
		}
		
		List<XY_DataSet> magFuncs = Lists.newArrayList();
		List<XY_DataSet> rateFuncs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		magFuncs.add(griddedLossVsMag);
		rateFuncs.add(griddedLossVsRate);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.GREEN.darker()));
		
		magFuncs.add(faultLossVsMag);
		rateFuncs.add(faultLossVsRate);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.RED.darker()));
		
		PlotSpec magSpec = new PlotSpec(magFuncs, chars, "Loss vs Mag", "Magnitude", yAxisLabel);
		magSpec.setLegendVisible(true);
		PlotSpec rateSpec = new PlotSpec(rateFuncs, chars, "Loss vs Rate", "Annual Rate", yAxisLabel);
		rateSpec.setLegendVisible(true);
		
		File outputDir = new File("/tmp");
		
		// now plot
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		// set user bounds?
//		gp.setUserBounds(new Range(0d, durationYears), null);

		ETAS_MultiSimAnalysisTools.setFontSizes(gp);

		gp.drawGraphPanel(magSpec, false, true);
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(new File(outputDir, "mag_losses.png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, "mag_losses.pdf").getAbsolutePath());
		
		gp.drawGraphPanel(rateSpec, true, true);
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(new File(outputDir, "rate_losses.png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, "rate_losses.pdf").getAbsolutePath());
		
		// now remove gridded
		rateFuncs.remove(0);
		magFuncs.remove(0);
		chars.remove(0);
		chars.set(0, new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLACK));
		magSpec.setLegendVisible(false);
		rateSpec.setLegendVisible(false);

		gp.drawGraphPanel(magSpec, false, true);
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(new File(outputDir, "mag_losses_fault_only.png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, "mag_losses_fault_only.pdf").getAbsolutePath());
		
		gp.drawGraphPanel(rateSpec, true, true);
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(new File(outputDir, "rate_losses_fault_only.png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, "rate_losses_fault_only.pdf").getAbsolutePath());
	}

}
