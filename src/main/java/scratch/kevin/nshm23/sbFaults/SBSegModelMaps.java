package scratch.kevin.nshm23.sbFaults;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SlipRateSegmentationConstraint.RateCombiner;
import org.opensha.sha.earthquake.faultSysSolution.modules.AveSlipModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectNuclMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SolMFDPlot;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.ClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SegmentationCalculator;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

import scratch.kevin.prvi25.figures.SlipRateFigures;

public class SBSegModelMaps {

	public static void main(String[] args) throws IOException {
		File sbDir = new File("/home/kevin/OpenSHA/2025_sb_faults");
		Region region = Region.fromFeature(Feature.read(new File(sbDir, "region.geojson")));
		
		File particOutputDir = new File(sbDir, "fault_partic");
		Preconditions.checkState(particOutputDir.exists() || particOutputDir.mkdir());
		File segOutputDir = new File(sbDir, "segmentation");
		Preconditions.checkState(segOutputDir.exists() || segOutputDir.mkdir());
		
		List<FaultSystemSolution> sols = new ArrayList<>();
		List<String> titles = new ArrayList<>();
		List<String> prefixes = new ArrayList<>();
		
		File origDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3");
		File nodeBADir = new File(origDir, "node_branch_averaged");
		
		sols.add(FaultSystemSolution.load(new File(origDir, "results_WUS_FM_v3_branch_averaged_mod_gridded.zip")));
		titles.add("Branch-averaged");
		prefixes.add("full_ba");
		
		sols.add(FaultSystemSolution.load(new File(nodeBADir, "SegModel_Classic.zip")));
		titles.add("Segmentation: Classic");
		prefixes.add("seg_classic");
		
		sols.add(FaultSystemSolution.load(new File(nodeBADir, "SegModel_High.zip")));
		titles.add("Segmentation: High");
		prefixes.add("seg_high");
		
		sols.add(FaultSystemSolution.load(new File(nodeBADir, "SegModel_Middle.zip")));
		titles.add("Segmentation: Middle");
		prefixes.add("seg_middle");
		
		sols.add(FaultSystemSolution.load(new File(nodeBADir, "SegModel_Low.zip")));
		titles.add("Segmentation: Low");
		prefixes.add("seg_low");
		
		sols.add(FaultSystemSolution.load(new File(nodeBADir, "SegModel_None.zip")));
		titles.add("Segmentation: None");
		prefixes.add("seg_none");
		
		CPT particCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-5, -2);
		particCPT.setBelowMinColor(particCPT.getMinColor());
		particCPT.setNanColor(Color.GRAY);
		particCPT.setLog10(true);
		
		CPT diffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-100, 100);
		diffCPT.setNanColor(Color.GRAY);
		
		CPT segCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, 1d);
		segCPT.setPreferredTickInterval(0.1);
		
		double[] minMags = { 0, 6, 6.7, 7, 7.5, 7.8 };
		
		DecimalFormat oDF = new DecimalFormat("0.#");
		DecimalFormat pDF = new DecimalFormat("0.0%");
		
		CSVFile<String> ratesCSV = new CSVFile<>(true);
		List<String> header = new ArrayList<>();
		header.add("Segmentation Model");
		for (double minMag : minMags) {
			String magLabel;
			if (minMag > 0d)
				magLabel = "M≥"+oDF.format(minMag);
			else
				magLabel = "Supra-seismogenic";
			header.add(magLabel+" rate");
			header.add(magLabel+" multifault rupture %");
		}
		header.add("Multifault rupture moment %");
		ratesCSV.addLine(header);
		
		FaultSystemSolution sol0 = sols.get(0);
		GeographicMapMaker mapMaker = new RupSetMapMaker(sol0.getRupSet(), region);
		mapMaker.setWriteGeoJSON(false);
		mapMaker.setSectOutlineChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, new Color(210, 210, 210)));
		mapMaker.setSectTraceChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, Color.GRAY));
		mapMaker.setPoliticalBoundaryChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		mapMaker.setJumpLineThickness(3f);
		
		List<double[]> rates0 = new ArrayList<>();
		
		ClusterRuptures cRups = sol0.getRupSet().requireModule(ClusterRuptures.class);
		PlausibilityConfiguration config = sol0.getRupSet().requireModule(PlausibilityConfiguration.class);
		ClusterConnectionStrategy connStrat = config.getConnectionStrategy();
		
		double[] regRupFract = sol0.getRupSet().getFractRupsInsideRegion(region, false);
		
		for (int s=0; s<sols.size(); s++) {
			FaultSystemSolution sol = sols.get(s);
			RupMFDsModule mfds = sol.requireModule(RupMFDsModule.class);
			FaultSystemRupSet rupSet = sol.getRupSet();
			AveSlipModule slips = rupSet.requireModule(AveSlipModule.class);
			String title = titles.get(s);
			String solPrefix = prefixes.get(s);
			
			double totMoment = 0d;
			double totMultiMoment = 0d;
			
			List<String> csvLine = new ArrayList<>();
			csvLine.add(title);
			
			// participation
			for (int m=0; m<minMags.length; m++) {
				double minMag = minMags[m];
				String label, prefix;
				
				double totRate = 0d;
				double totMultiRate = 0d;
				
				for (int rupIndex=0; rupIndex<rupSet.getNumRuptures(); rupIndex++) {
					if (regRupFract[rupIndex] > 0d) {
						List<Double> mags = new ArrayList<>();
						List<Double> rates = new ArrayList<>();
						DiscretizedFunc mfd = mfds.getRuptureMFD(rupIndex);
						if (mfd == null) {
							mags.add(rupSet.getMagForRup(rupIndex));
							rates.add(sol.getRateForRup(rupIndex));
						} else {
							for (Point2D pt : mfd) {
								mags.add(pt.getX());
								rates.add(pt.getY());
							}
						}
						
						double rupMoment = FaultMomentCalc.getMoment(rupSet.getAreaForRup(rupIndex), slips.getAveSlip(rupIndex));
						for (int r=0; r<rates.size(); r++) {
							double mag = mags.get(r);
							double rate = rates.get(r);
							
							if (mag >= minMag) {
								double nuclRate = rate*regRupFract[rupIndex];
								boolean multifault = false;
								List<FaultSection> sects = rupSet.getFaultSectionDataForRupture(rupIndex);
								int parent0 = sects.get(0).getParentSectionId();
								for (FaultSection sect : sects) {
									if (sect.getParentSectionId() != parent0) {
										multifault = true;
										break;
									}
								}
								totRate += nuclRate;
								if (multifault)
									totMultiRate += nuclRate;
								if (minMag == 0d) {
									double momentRateInReg = nuclRate*rupMoment;
									totMoment += momentRateInReg;
									if (multifault)
										totMultiMoment += momentRateInReg;
								}
							}
						}
					}
				}
				csvLine.add(String.format("%.3g", totRate));
				csvLine.add(pDF.format(totMultiRate/totRate));
				
				if (minMag > 0d) {
					label = "M≥"+oDF.format(minMag);
					prefix = "m"+oDF.format(minMag);
				} else {
					label = "Supra-Seismogenic";
					prefix = "supra_seis";
				}
				prefix = solPrefix+"_paric_"+prefix;
				String diffLabel = label +" Participation Rate (% Change)";
				label = label+" Participation Rate (/yr)";
				
				double[] rates = sol.calcParticRateForAllSects(minMag, Double.POSITIVE_INFINITY);
				double[] pDiffs, refRates;
				if (s == 0) {
					rates0.add(Arrays.copyOf(rates, rates.length));
					pDiffs = null;
					refRates = null;
				} else {
					pDiffs = new double[rates.length];
					refRates = rates0.get(m);
				}
				for (int i=0; i<rates.length; i++) {
					if (s > 0) {
						if (refRates[i] == 0d)
							pDiffs[i] = Double.NaN;
						else
							pDiffs[i] = 100d*(rates[i]-refRates[i])/refRates[i];
					}
					if (rates[i] == 0d)
						rates[i] = Double.NaN;
				}
				
				mapMaker.plotSectScalars(rates, particCPT, label);
				
				mapMaker.plot(particOutputDir, prefix, title);
				
				if (s > 0) {
					mapMaker.plotSectScalars(pDiffs, diffCPT, diffLabel);
					
					mapMaker.plot(particOutputDir, prefix+"_pdiff", title);
				}
			}
			
			csvLine.add(pDF.format(totMultiMoment/totMoment));
			ratesCSV.addLine(csvLine);
			
			mapMaker.clearSectScalars();
			
			SegmentationCalculator segCalc = new SegmentationCalculator(sol, cRups.getAll(),
					connStrat, config.getDistAzCalc(), new double[] { 0d });
			
			segCalc.setLegendFontSize(26);
			
//			CPT cpt = SegmentationCalculator.getConnectionFractCPT();
			
			mapMaker.plotJumpScalars(segCalc.calcJumpPassthroughs(0, RateCombiner.MIN), segCPT, SegmentationCalculator.PASSTHROUGH_LABEL);
			
			mapMaker.plot(segOutputDir, solPrefix+"_passthrough_map", title, 1400);
			
			if (s == 0) {
				// plot distances
				HashSet<Integer> parentsInside = new HashSet<>();
				double[] sectFracts = rupSet.getFractSectsInsideRegion(region, false);
				for (FaultSection sect : rupSet.getFaultSectionDataList())
					if (sectFracts[sect.getSectionId()] > 0d)
						parentsInside.add(sect.getParentSectionId());
				Map<Jump, Double> jumpDists = new HashMap<>();
				List<Jump> jumpsInside = new ArrayList<>();
				for (Jump jump : segCalc.calcJumpPassthroughs(0, RateCombiner.MIN).keySet()) {
					if (jumpDists.containsKey(jump.reverse()))
						continue;
					jumpDists.put(jump, jump.distance);
					boolean contained = parentsInside.contains(jump.fromSection.getParentSectionId())
							|| parentsInside.contains(jump.toSection.getParentSectionId());
					if (contained)
						jumpsInside.add(jump);
				}
				
				HistogramFunction jumpBinning = HistogramFunction.getEncompassingHistogram(0.01, 14.99, 1d);
				int[] jumpCounts = new int[jumpBinning.size()];
				Map<Integer, Double> parentMinDists = new HashMap<>();
				for (Jump jump : jumpsInside) {
					int parent1 = jump.fromSection.getParentSectionId();
					int parent2 = jump.toSection.getParentSectionId();
					Double parent1Min = parentMinDists.get(parent1);
					Double parent2Min = parentMinDists.get(parent2);
					if (parent1Min == null || jump.distance < parent1Min)
						parentMinDists.put(parent1, jump.distance);
					if (parent2Min == null || jump.distance < parent2Min)
						parentMinDists.put(parent2, jump.distance);
					int index = jumpBinning.getClosestXIndex(jump.distance);
					jumpCounts[index]++;
				}
				int[] sectsWithDistCounts = new int[jumpBinning.size()];
				for (int parentID : parentsInside) {
					Double minDist = parentMinDists.get(parentID);
					if (minDist != null) {
						int index = jumpBinning.getClosestXIndex(minDist);
						for (int i=index; i<jumpBinning.size(); i++)
							sectsWithDistCounts[i]++;
					}
				}
				CSVFile<String> jumpCSV = new CSVFile<>(true);
				jumpCSV.addLine("Bin Start", "Bin End", "Jump Count", "Jump %", "Parents with <=", "% Parents with <=");
				for (int d=0; d<jumpBinning.size(); d++) {
					double middle = jumpBinning.getX(d);
					double start = middle - 0.5;
					double end = middle + 0.5;
					double jumpFract = (double)jumpCounts[d]/jumpsInside.size();
					double parentFract = (double)sectsWithDistCounts[d]/parentsInside.size();
					jumpCSV.addLine((float)start+"", (float)end+"", jumpCounts[d]+"", pDF.format(jumpFract),
							sectsWithDistCounts[d]+"", pDF.format(parentFract));
				}
				jumpCSV.writeToFile(new File(segOutputDir, "jump_dist_stats.csv"));
				
				CPT distCPT = GMT_CPT_Files.BLACK_RED_YELLOW_UNIFORM.instance().rescale(0d, 15d);
				distCPT = distCPT.asDiscrete(1d, true);
				
				mapMaker.plotJumpScalars(jumpDists, distCPT, "Jump Distance (km)");
				
				mapMaker.plot(segOutputDir, solPrefix+"_jump_dists", title, 1400);
			}
			
			mapMaker.clearJumpScalars();
		}
		ratesCSV.writeToFile(new File(particOutputDir, "rates.csv"));
	}

}
