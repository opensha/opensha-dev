package scratch.kevin.mfdInversion;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.UncertainDataConstraint.SectMappedUncertainDataConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.NamedFaults;
import org.opensha.sha.earthquake.faultSysSolution.modules.PaleoseismicConstraintData;
import org.opensha.sha.earthquake.faultSysSolution.modules.PosteriorSectionBValueDistributions;
import org.opensha.sha.earthquake.nshmp.inversion.mfdPreInversion.PaleoBValueEstimator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SingleStates;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;
import scratch.kevin.nshm23.figures.PaleoZTablesBuilder;
import scratch.kevin.nshm23.figures.SlipZTablesBuilder;
import scratch.kevin.nshm23.figures.SlipZTablesBuilder.SlipZRecord;

public class PaleoZComparisons {

	public static void main(String[] args) throws IOException {
		List<File> solFiles = new ArrayList<>();
		List<String> names = new ArrayList<>();
		List<String> prefixes = new ArrayList<>();
		
		Region[] regions = {
				NSHM23_RegionLoader.AnalysisRegions.CONUS_U3_RELM.load(),
				NSHM23_RegionLoader.LocalRegions.CONUS_SF_BAY.load(),
				NSHM23_RegionLoader.LocalRegions.CONUS_LA_BASIN.load(),
				NSHM23_SingleStates.UT.loadRegion(),
				NSHM23_RegionLoader.LocalRegions.CONUS_WASATCH.load(),
		};
		String[] regionNames = {
				"CA",
				"SF",
				"LA",
				"UT",
				"Wasatch"
		};
		File outputDir = new File("/home/kevin/OpenSHA/nshm23/paleo_b_value/run_comparisons");
		
		File invDir = new File("/home/kevin/OpenSHA/fss_inversions");
		
		solFiles.add(new File(invDir, "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged.zip"));
		names.add("NSHM23 Full BA");
		prefixes.add("nshm23_full_ba");
		
		solFiles.add(new File(invDir, "2026_08_14-nshm23_branches-WUS_FM_v3-bUniform10x-1000samples/results_WUS_FM_v3_branch_averaged.zip"));
		names.add("Uniform-b BA");
		prefixes.add("uniform_b");
		
//		solFiles.add(new File(invDir, "2024_02_02-nshm23_branches-WUS_FM_v3/node_branch_averaged/PaleoUncert_EvenFitPaleo.zip"));
//		names.add("NSHM23 Even-Fit Paleo");
		
		solFiles.add(new File(invDir, "2026_08_18-nshm23_branches-WUS_FM_v3-bPosterior10x-students_t_5-comb_likelihood-paleo_ris-1000samples/results_WUS_FM_v3_branch_averaged.zip"));
		names.add("Posterior-b Full BA");
		prefixes.add("posterior_b");
		
		solFiles.add(new File(invDir, "2026_08_19-nshm23_branches-WUS_FM_v3-bUniform10x-smooth_paleo_visible-1000samples/results_WUS_FM_v3_branch_averaged.zip"));
		names.add("Uniform-b BA (paleo-visible smooth)");
		prefixes.add("uniform_b_smooth");
		
		solFiles.add(new File(invDir, "2026_08_19-nshm23_branches-WUS_FM_v3-bPosterior10x-students_t_5-comb_likelihood-paleo_ris-smooth_paleo_visible-1000samples/results_WUS_FM_v3_branch_averaged.zip"));
		names.add("Posterior-b Full BA (paleo-visible smooth)");
		prefixes.add("posterior_b_smooth");
		
		List<CSVFile<String>> csvs = new ArrayList<>();
		List<String> header = new ArrayList<>();
		header.add("Model");
		header.add("Average slip z");
		header.add("Improvement");
		header.add("Mo-weighted slip z");
		header.add("Improvement");
		header.add("Paleo-parent slip z");
		header.add("Improvement");
		header.add("Paleo z");
		header.add("Paleo |z|");
		header.add("Improvement");
		for (int r=0; r<regions.length; r++) {
			CSVFile<String> csv = new CSVFile<>(true);
			csv.addLine(header);
			csvs.add(csv);
		}
		
//		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-10, 10).expandAt(-2d, 2d);
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-5, 5).expandAt(-1d, 1d);
//		CPT riDiffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().rescale(-20, 20);
		
//		double[] minMags = {0d, 7d};
//		String[] minMagPrefixes = {"supra", "m7"};
//		String[] minMagLabels = {"Supra-seis", "M>7"};
		double[] minMags = {0d};
		String[] minMagPrefixes = {"supra"};
		String[] minMagLabels = {"Supra-seis"};
		double[][] prevRates = null;
		
		Double[] prevSlipZs = new Double[regions.length];
		Double[] prevMoWeightedSlipZs = new Double[regions.length];
		Double[] prevSlipPaleoParentZs = new Double[regions.length];
		Double[] prevPaleoZs = new Double[regions.length];
		for (int s=0; s<solFiles.size(); s++) {
			FaultSystemSolution sol = FaultSystemSolution.load(solFiles.get(s));
			
			FaultSystemRupSet rupSet = sol.getRupSet();
			HashSet<Integer> paleoParents = new HashSet<>();
			List<? extends SectMappedUncertainDataConstraint> paleoConstraints =
					rupSet.requireModule(PaleoseismicConstraintData.class).getPaleoRateConstraints();
			for (SectMappedUncertainDataConstraint data : paleoConstraints)
				paleoParents.add(rupSet.getFaultSectionData(data.sectionIndex).getParentSectionId());
			NamedFaults namedFaults = rupSet.requireModule(NamedFaults.class);
			// add in neighbors on named faults
			for (List<Integer> parents : namedFaults.get().values()) {
				boolean any = false;
				for (int parentID : parents)
					any |= paleoParents.contains(parentID);
				if (any)
					paleoParents.addAll(parents);
			}
			
			double[][] magRates = new double[minMags.length][];
			for (int m=0; m<minMags.length; m++) {
				magRates[m] = sol.calcParticRateForAllSects(minMags[m], Double.POSITIVE_INFINITY);
				for (int i=0; i<magRates[m].length; i++)
					if (!paleoParents.contains(rupSet.getFaultSectionData(i).getParentSectionId()))
						magRates[m][i] = Double.NaN;
			}
			
			for (int r=0; r<regions.length; r++) {
				List<String> line = new ArrayList<>();
				line.add(names.get(s));
				boolean[] sectsInRegion = new boolean[rupSet.getNumSections()];
				double[] fractsInRegion = rupSet.getFractSectsInsideRegion(regions[r], false);
				for (int i=0; i<fractsInRegion.length; i++)
					sectsInRegion[i] = fractsInRegion[i] > 0;
				SlipZRecord slips = SlipZTablesBuilder.calcSlipZ(sol, NSHM23_FaultModels.WUS_FM_v3,
						NSHM23_DeformationModels.AVERAGE, true, sectsInRegion);
				line.add(zDF.format(slips.average()));
				line.add(improvementStr(slips.average(), prevSlipZs[r]));
				line.add(zDF.format(slips.momentWeightedAverage()));
				line.add(improvementStr(slips.momentWeightedAverage(), prevMoWeightedSlipZs[r]));
				prevSlipZs[r] = slips.average();
				prevMoWeightedSlipZs[r] = slips.momentWeightedAverage();
				
				// now paleo parents only
				for (int i=0; i<sectsInRegion.length; i++)
					sectsInRegion[i] &= paleoParents.contains(rupSet.getFaultSectionData(i).getParentSectionId());
				slips = SlipZTablesBuilder.calcSlipZ(sol, NSHM23_FaultModels.WUS_FM_v3,
						NSHM23_DeformationModels.AVERAGE, true, sectsInRegion);
				line.add(zDF.format(slips.average()));
				line.add(improvementStr(slips.average(), prevSlipPaleoParentZs[r]));
				prevSlipPaleoParentZs[r] = slips.average();
				
				double paleoZ = PaleoZTablesBuilder.avgZScore(sol, false, regions[r]);
				double paleoAbsZ = PaleoZTablesBuilder.avgZScore(sol, true, regions[r]);
				line.add(zDF.format(paleoZ));
				line.add(zDF.format(paleoAbsZ));
				line.add(improvementStr(paleoAbsZ, prevPaleoZs[r]));
				prevPaleoZs[r] = paleoAbsZ;
				csvs.get(r).addLine(line);
			}
			
			if (prevRates != null) {
				GeographicMapMaker mapMaker = new GeographicMapMaker(rupSet.getFaultSectionDataList());
				mapMaker.setPrintSizingDefault(true);
				mapMaker.setDefaultPrintDPI(300);
				mapMaker.setWriteGeoJSON(false);
				mapMaker.setSkipNaNs(true);
				mapMaker.setSectNaNChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 0.5f, Color.GRAY));
				mapMaker.setSectOutlineChar(null);
				
				List<Location> scatterLocs = new ArrayList<>();
				for (SectMappedUncertainDataConstraint data : paleoConstraints)
					scatterLocs.add(data.dataLocation);
//				mapMaker.plotScatters(scatterLocs, Colors.tab_green);
//				mapMaker.setScatterSymbol(PlotSymbol.TRIANGLE, 4f);
				mapMaker.plotScatters(scatterLocs, Color.DARK_GRAY);
//				mapMaker.setScatterSymbol(PlotSymbol.FILLED_CIRCLE, 2f);
				mapMaker.setScatterSymbol(PlotSymbol.FILLED_CIRCLE, 5f, PlotSymbol.CIRCLE, Color.WHITE);
				mapMaker.setAbsoluteSort(true);
				
				String modelsPrefix = prefixes.get(s)+"_vs_"+prefixes.get(s-1);
				
				for (int m=0; m<minMags.length; m++) {
					double[] riDiffs = new double[rupSet.getNumSections()];
					double[] ratePDiffs = new double[rupSet.getNumSections()];
					
					for (int i=0; i<riDiffs.length; i++) {
						double prevRate = prevRates[m][i];
						double rate = magRates[m][i];
						
						ratePDiffs[i] = 100d*(rate - prevRate)/prevRate;
						riDiffs[i] = (1d/rate) - (1d/prevRate);
					}
					
					String magPrefix = modelsPrefix+"_"+minMagPrefixes[m];
					
					for (int r=0; r<regions.length; r++) {
						mapMaker.setRegion(regions[r]);
						
						String prefix = magPrefix+"_"+regionNames[r];
						
						mapMaker.plotSectScalars(ratePDiffs, pDiffCPT, minMagLabels[m]+" participation rate, % change");
						mapMaker.plot(outputDir, prefix+"_rates", " ");
//						mapMaker.plotSectScalars(riDiffs, riDiffCPT, minMagLabels[m]+" partic RI change (years)");
//						mapMaker.plot(outputDir, prefix+"_ris", " ");
					}
				}
			}
			
			prevRates = magRates;
			
			if (rupSet.hasModule(PosteriorSectionBValueDistributions.class)) {
				PosteriorSectionBValueDistributions posteriors = rupSet.requireModule(PosteriorSectionBValueDistributions.class);
				File solDir = new File(outputDir, prefixes.get(s)+"_posteriors");
				Preconditions.checkState(solDir.exists() || solDir.mkdir());
				PaleoBValueEstimator.plotSectDistributions(solDir, rupSet, posteriors);
				PaleoBValueEstimator.plotFaultSiteAverages(solDir, "fault_avg", rupSet, posteriors, namedFaults, false);
				PaleoBValueEstimator.plotSolAlongStrike(solDir, "along_strike", sol, posteriors);
			}
		}
		for (int r=0; r<regions.length; r++)
			csvs.get(r).writeToFile(new File(outputDir, "z_scores_"+regionNames[r]+".csv"));
	}
	private static DecimalFormat zDF = new DecimalFormat("0.0000");
	private static DecimalFormat pDF = new DecimalFormat("0.00%");
	
	private static String improvementStr(double absZ, Double prevAbsZ) {
		if (prevAbsZ == null)
			return "";
		double diff = prevAbsZ - absZ;
		return zDF.format(diff)+" ("+pDF.format(diff/prevAbsZ)+")";
	}

}
