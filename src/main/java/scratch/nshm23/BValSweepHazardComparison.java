package scratch.nshm23;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportMetadata;
import org.opensha.sha.earthquake.faultSysSolution.reports.RupSetMetadata;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.HazardMapPlot;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.imr.AttenRelRef;

import com.google.common.base.Preconditions;

public class BValSweepHazardComparison {

	public static void main(String[] args) throws IOException {
		File inversionDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
//		try {
//			Thread.sleep(1000l*900l); // 15 m
//		} catch (InterruptedException e) {}
		
		double[] bVals = { 0d, 0.1d, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1d };
//		double[] bVals = { 0d, 0.3, 0.8, 1d };
		
//		File compDir = new File(inversionDir,
////				"2021_10_18-reproduce-ucerf3-ref_branch-uniform-new_anneal-5x_avg-try_zero-var_perturb-u3WL-5h");
//				"2021_11_21-reproduce-ucerf3-ref_branch-uniform-u3_constr-new_anneal-no_u2_ss_mfds-single_mfd_region-u3WL-5h");
////				"2021_11_21-reproduce-ucerf3-ref_branch-uniform-u3_constr-new_anneal-single_mfd_region-u3WL-5h");
////		String sweepPrefix = "2021_11_18-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_";
////		String sweepSuffix = "-adj_ucert_for_data-paleo_wt_5-parkfield_wt_100-mfd_wt_10-sect_wt_0.5-smooth_paleo_wt_10000-skipBelow-2h";
//		String sweepPrefix = "2021_11_19-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_";
//		String sweepSuffix = "-adj_ucert_for_data-u3_supra_reduction-paleo_wt_5-parkfield_wt_100-mfd_wt_10-sect_wt_0.5-smooth_paleo_wt_10000-skipBelow-2h";
////		String sweepPrefix = "2021_11_19-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_";
////		String sweepSuffix = "-adj_ucert_for_data-u3_supra_reduction-paleo_wt_5-parkfield_wt_100-mfd_wt_10-no_sect_rate-smooth_paleo_wt_10000-skipBelow-10h";
		
		File compDir = new File(inversionDir,
//				"2021_11_20-reproduce-ucerf3-ref_branch-uniform-u3_constr-new_anneal-no_paleo-no_parkfield-no_smooth-u3WL-5h");
//				"2021_11_20-reproduce-ucerf3-ref_branch-uniform-u3_constr-new_anneal-no_paleo-no_parkfield-no_u2_ss_mfds-noWL-5h");
//				"2021_11_21-reproduce-ucerf3-ref_branch-uniform-u3_constr-new_anneal-no_paleo-no_parkfield-no_smooth-single_mfd_region-u3WL-5h");
//				"2021_11_22-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_NaN-no_paleo-no_parkfield-no_mfd-no_sect_rate-u3_target_mfds-single_mfd_region-skipBelow-2h");
				"2021_11_22-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_NaN-no_paleo-no_parkfield-no_mfd-no_sect_rate-u3_target_mfds-skipBelow-2h");
		String sweepPrefix = "2021_11_19-reproduce-ucerf3-ref_branch-uniform-nshm23_draft-supra_b_";
		String sweepSuffix = "-u3_supra_reduction-no_paleo-no_parkfield-mfd_wt_10-sect_wt_0.5-skipBelow-2h";
		
		File sweepDir = new File(inversionDir, sweepPrefix+"sweep"+sweepSuffix);
		
		FaultSystemSolution compSol = FaultSystemSolution.load(new File(compDir, "mean_solution.zip"));
		
		double[] periods = { 0d, 0.2d, 1d };
		AttenRelRef gmpeRef = AttenRelRef.ASK_2014;;
		double gridSpacing = 0.1;
		
		Region region = new ReportMetadata(new RupSetMetadata(null, compSol)).region;
		
		GriddedRegion gridReg = new GriddedRegion(region, gridSpacing, GriddedRegion.ANCHOR_0_0);
		
		// calculate comparison hazard
		System.out.println("Calculating comparison hazard...");
		SolHazardMapCalc compCalc = getCalcCurves(new File(compDir, "hazard"), compSol, gmpeRef, periods, gridReg);
		
		ReturnPeriods[] rps = ReturnPeriods.values();
		
		GriddedGeoDataSet[][] compMaps = new GriddedGeoDataSet[periods.length][rps.length];
		System.out.println("Calculating comparison maps...");
		for (int p=0; p<periods.length; p++)
			for (int r=0; r<rps.length; r++)
				compMaps[p][r] = compCalc.buildMap(periods[p], rps[r]);
		
		GriddedGeoDataSet[][][] sweepMaps = new GriddedGeoDataSet[periods.length][rps.length][bVals.length];
		
		for (int b=0; b<bVals.length; b++) {
			double bVal = bVals[b];
			File subSweepDir = new File(sweepDir, sweepPrefix+(float)bVal+sweepSuffix);
			Preconditions.checkState(subSweepDir.exists());
			FaultSystemSolution sweepSol = FaultSystemSolution.load(new File(subSweepDir, "mean_solution.zip"));
			
			File hazardDir = new File(subSweepDir, "hazard");
			
			System.out.println("Calculating b="+(float)bVal+" sweep hazard...");
			SolHazardMapCalc sweepCalc = getCalcCurves(hazardDir, sweepSol, gmpeRef, periods, gridReg);
			
			System.out.println("Calculating maps...");
			for (int p=0; p<periods.length; p++)
				for (int r=0; r<rps.length; r++)
					sweepMaps[p][r][b] = sweepCalc.buildMap(periods[p], rps[r]);
			
			HazardMapPlot plot = new HazardMapPlot(gmpeRef, gridReg.getSpacing(), periods);
			
			File resourcesDir = new File(hazardDir, "resources");
			Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
			
			List<String> lines = new ArrayList<>();
			
			lines.add("# Hazard Maps");
			lines.add("");
			
			lines.addAll(plot.plot(resourcesDir, resourcesDir.getName(), "", gridReg, sweepCalc, compCalc));
			
			MarkdownUtils.writeReadmeAndHTML(lines, hazardDir);
		}
		
		File sweepHazDir = new File(sweepDir, "sweep_hazard");
		Preconditions.checkState(sweepHazDir.exists() || sweepHazDir.mkdir());
		
		File outputDir = new File(sweepHazDir, "comp_"+compDir.getName());
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		double minB = StatUtils.min(bVals);
		double maxB = StatUtils.max(bVals);
		CPT bValCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(minB, maxB);
		bValCPT.setBelowMinColor(bValCPT.getMinColor().darker());
		bValCPT.setAboveMaxColor(bValCPT.getMaxColor().darker());
		for (CPTVal c : bValCPT) {
			c.minColor = saturate(c.minColor);
			c.maxColor = saturate(c.maxColor);
		}
		bValCPT.add(0, new CPTVal((float)(minB-0.02), bValCPT.getBelowMinColor(), (float)minB, bValCPT.getBelowMinColor()));
		bValCPT.add(new CPTVal((float)maxB, bValCPT.getAboveMaxColor(), (float)(maxB+0.02), bValCPT.getAboveMaxColor()));
		CPT pDiffCPT = GMT_CPT_Files.GMT_POLAR.instance().rescale(-100d, 100d);
		
		DecimalFormat optionalDigitDF = new DecimalFormat("0.#");
		DecimalFormat percentDF = new DecimalFormat("0.0%");
		for (int p=0; p<periods.length; p++) {
			String perLabel, perUnits;
			if (periods[p] == -1d) {
				perLabel = "PGV";
				perUnits = "cm/s";
			} else if (periods[p] == 0d) {
				perLabel = "PGA";
				perUnits = "g";
			} else {
				Preconditions.checkState(periods[p] > 0);
				perLabel = optionalDigitDF.format(periods[p])+"s SA";
				perUnits = "g";
			}
			String perPrefix = perLabel.toLowerCase().replaceAll(" ", "_");
			for (int r=0; r<rps.length; r++) {
				GriddedGeoDataSet compMap = compMaps[p][r];
				GriddedGeoDataSet compWithin = new GriddedGeoDataSet(gridReg, false);
				GriddedGeoDataSet compPDiff = new GriddedGeoDataSet(gridReg, false);
				
				System.out.println("Computing for p="+(float)periods[p]+", rp="+rps[r]);
				
				int numWithin = 0;
				
				for (int i=0; i<compWithin.size(); i++) {
					ArbitrarilyDiscretizedFunc dist = new ArbitrarilyDiscretizedFunc();
					for (int b=0; b<bVals.length; b++)
						dist.set(bVals[b], sweepMaps[p][r][b].get(i));
					double compVal = compMap.get(i);
					double position, pDiff;
					if (compVal > dist.getMinY() && compVal < dist.getMaxY()) {
						position = dist.getFirstInterpolatedX(compVal);
						numWithin++;
						pDiff = 0d;
					} else {
						double lowBVal = dist.getY(0);
						double highBVal = dist.getY(dist.size()-1);
						double distLower = Math.abs(compVal - lowBVal);
						double distUpper = Math.abs(compVal - highBVal);
						if (distLower < distUpper) {
							position = Double.NEGATIVE_INFINITY;
							pDiff = 100d*(compVal-lowBVal)/lowBVal;
						} else {
							position = Double.POSITIVE_INFINITY;
							pDiff = 100d*(compVal-highBVal)/highBVal;
						}
					}
					compWithin.set(i, position);
					compPDiff.set(i, pDiff);
				}
				
				System.out.println("\t"+numWithin+"/"+compWithin.size()+" within b-value sweep hazard ("
						+percentDF.format((double)numWithin/(double)compWithin.size())+")");
				
				String prefix = "b_val_seep_map_"+perPrefix+"_"+rps[r].name().toLowerCase();
				String zLabel = "b-value w/ Closest U3 Hazard, "+perLabel+", "+rps[r].label;
				compCalc.plotMap(outputDir, prefix+"_b_dist", compWithin, bValCPT, "b-Value Sweep/UCERF3 Comparison", zLabel);
				
				zLabel = "% Difference from Closest Sweep Map, "+perLabel+", "+rps[r].label;
				compCalc.plotMap(outputDir, prefix+"_pDiff", compPDiff, pDiffCPT, "b-Value Sweep/UCERF3 Comparison", zLabel);
			}
		}
	}
	
	private static SolHazardMapCalc getCalcCurves(File hazardDir, FaultSystemSolution sol, AttenRelRef gmpeRef,
			double[] periods, GriddedRegion gridReg) throws IOException {
		// see if we already have curves
		SolHazardMapCalc calc = null;
		if (hazardDir.exists()) {
			try {
				calc = SolHazardMapCalc.loadCurves(sol, gridReg, periods, hazardDir, "curves");
				System.out.println("Loaded existing curves!");
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		if (calc == null) {
			if (!hazardDir.exists())
				hazardDir.mkdir();
			
			// need to calculate
			calc = new SolHazardMapCalc(sol, gmpeRef, gridReg, periods);
			
			calc.calcHazardCurves(FaultSysTools.defaultNumThreads());
			
			calc.writeCurvesCSVs(hazardDir, "curves", gridReg.getSpacing() < 0.1d);
		}
		return calc;
	}
	
	private static final int saturation_steps = 1;
	
	private static Color saturate(Color c) {
		int r = c.getRed();
		int g = c.getGreen();
		int b = c.getBlue();
		
		for (int i=0; i<saturation_steps; i++) {
			r = (int)(0.5d*(r + 255d)+0.5);
			g = (int)(0.5d*(g + 255d)+0.5);
			b = (int)(0.5d*(b + 255d)+0.5);
		}
		
		return new Color(r, g, b, c.getAlpha());
	}

}
