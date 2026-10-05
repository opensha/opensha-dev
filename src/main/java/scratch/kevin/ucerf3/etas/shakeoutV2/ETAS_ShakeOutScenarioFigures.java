package scratch.kevin.ucerf3.etas.shakeoutV2;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.function.DoubleUnaryOperator;

import org.jfree.data.Range;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ColorUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.faultSurface.FaultSection;

import net.mahdilamb.colormap.Colors;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.ETAS_Catalog;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.launcher.ETAS_Config;
import scratch.UCERF3.erf.ETAS.launcher.ETAS_Launcher;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;

public class ETAS_ShakeOutScenarioFigures {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/home/kevin/OpenSHA/UCERF3/etas/shakeout_v2");
		File simsDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations");
		File[] scenarioDirs = {
				new File(simsDir, "2026_05_27-FSS_Rupture_201887_M7p8_Start_2026_10_15_1_yr_kCOV_1p5_MaxPtSrcM_6"),
				new File(simsDir, "2026_06_12-FSS_Rupture_201887_M7p8_Start_2026_10_15_1_yr_kCOV_1p5_MaxPtSrcM_6")
		};
		
		ETAS_Config config1 = ETAS_Config.readJSON(new File(scenarioDirs[0], "config.json"));
		
		ETAS_Launcher launcher = new ETAS_Launcher(config1, false);
		
		FaultSystemSolution sol = launcher.checkOutFSS();
		
		List<? extends FaultSection> sects = sol.getRupSet().getFaultSectionDataList();
		
		String[] parentNames = {"Hollywood", "Raymond"};
		int[] scenarioParents = new int[parentNames.length];
		for (int p=0; p<parentNames.length; p++)
			scenarioParents[p] = FaultSectionUtils.findParentSectionID(sects, parentNames[p]);
		
		int inputScenarioID = 201887;
		int targetScenarioID = 218331;
		
//		Region plotReg = new Region(new Location(32, -121), new Location(36, -115));
		Region plotReg = new Region(new Location(32.5, -121), new Location(35.5, -115.5));
		GriddedRegion gridReg = new GriddedRegion(plotReg, 0.02, GriddedRegion.ANCHOR_0_0);
		
		double[] minMags = {5d, 6d, 7d};
		GriddedGeoDataSet[] nuclXYZ = new GriddedGeoDataSet[minMags.length];
		double[] sectPartics = new double[sects.size()];
		GriddedRegion simReg = launcher.getRegion();
		for (int m=0; m<minMags.length; m++) {
			nuclXYZ[m] = new GriddedGeoDataSet(gridReg);
			for (int i=0; i<nuclXYZ[m].size(); i++) {
				if (simReg.indexForLocation(nuclXYZ[m].getLocation(i)) < 0)
					nuclXYZ[m].set(i, Double.NaN);
			}
		}
		
		EvenlyDiscretizedFunc refWeeks = new EvenlyDiscretizedFunc(0.5, 52, 1d);
		EvenlyDiscretizedFunc refDays = new EvenlyDiscretizedFunc(0.5, 365, 1d);

		EvenlyDiscretizedFunc targetWeeks = refWeeks.deepClone();
		EvenlyDiscretizedFunc targetDays = refDays.deepClone();
		EvenlyDiscretizedFunc corupWeeks = refWeeks.deepClone();
		EvenlyDiscretizedFunc corupDays = refDays.deepClone();
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(gridReg);
		mapMaker.setFaultSections(sects);
		mapMaker.setSectOutlineChar(null);
		mapMaker.setWriteGeoJSON(false);
		mapMaker.setPoliticalBoundaryChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, Color.BLACK));
		
		PlotCurveCharacterstics otherSectChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY);
		
//		CPT magCPT = GMT_CPT_Files.SEQUENTIAL_OSLO_UNIFORM.instance().rescale(5d, 7.5d);
//		DoubleUnaryOperator magSizeFunc = M -> 3d + 3*(M-5d);
		
		CPT magCPT = GMT_CPT_Files.SEQUENTIAL_NAVIA_UNIFORM.instance().rescale(5d, 8d).trim(5d, 7.5d);
		DoubleUnaryOperator magSizeFunc = M -> 3d + 8*(M-5d);
		
		int numCatalogs = 0;
		int numCatalogsWithTarget = 0;
		int[] numCatalogsWithParents = new int[scenarioParents.length];
		int numCatalogsWithParentsCorup = 0;
		int[] catMagCounts = new int[minMags.length];
		for (File dir : scenarioDirs) {
			File binFile = new File(dir, "results_m5_preserve_chain.bin");
			
			for (ETAS_Catalog catalog : ETAS_CatalogIO.getBinaryCatalogsIterable(binFile, 5d)) {
				numCatalogs++;
				boolean hasTarget = false;
				boolean[] hasParents = new boolean[scenarioParents.length];
				boolean hasCorup = false;
				boolean[] partics = new boolean[sects.size()];
				boolean[] hasMags = new boolean[minMags.length];
				for (ETAS_EqkRupture rup : catalog) {
					Location hypo = rup.getHypocenterLocation();
					double mag = rup.getMag();
					for (int m=0; m<minMags.length; m++)
						hasMags[m] |= mag >= minMags[m];
					int gridIndex = gridReg.indexForLocation(hypo);
					if (gridIndex >= 0) {
						for (int m=0; m<minMags.length; m++) {
							if (mag >= minMags[m])
								nuclXYZ[m].add(gridIndex, 1d);
						}
					}
					
					if (rup.getFSSIndex() >= 0) {
						boolean target = rup.getFSSIndex() == targetScenarioID;
						List<FaultSection> rupSects = sol.getRupSet().getFaultSectionDataForRupture(rup.getFSSIndex());
						boolean[] myHasParents = new boolean[scenarioParents.length];
						for (FaultSection sect : rupSects) {
							int parent = sect.getParentSectionId();
							for (int p=0; p<scenarioParents.length; p++)
								myHasParents[p] |= parent == scenarioParents[p];
							partics[sect.getSectionId()] = true;
						}
						boolean all = true;
						for (int p=0; p<scenarioParents.length; p++) {
							all &= myHasParents[p];
							hasParents[p] |= myHasParents[p];
						}
						if (target || all) {
							long timeDelta = rup.getOriginTime() - config1.getSimulationStartTimeMillis();
							double days = (double)timeDelta/(double)ProbabilityModelsCalc.MILLISEC_PER_DAY;
							int daysIndex = refDays.getClosestXIndex(days);
							double weeks = days / 7d;
							int weeksIndex = refWeeks.getClosestXIndex(weeks);
							if (target && !hasTarget) {
								targetWeeks.add(weeksIndex, 1d);
								targetDays.add(daysIndex, 1d);
							}
							if (all && !hasCorup) {
								corupWeeks.add(weeksIndex, 1d);
								corupDays.add(daysIndex, 1d);
							}
						}
						hasCorup |= all;
						hasTarget |= target;
					}
				}
				
				if (hasTarget) {
					LocationList rupLocs = new LocationList();
					List<PlotCurveCharacterstics> chars = new ArrayList<>();
					
					double[] maxMags = new double[sects.size()];
					
					for (ETAS_EqkRupture rup : catalog) {
						rupLocs.add(rup.getHypocenterLocation());
						double mag = rup.getMag();
						chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE,
								(float)magSizeFunc.applyAsDouble(mag), magCPT.getColor(mag)));
						if (rup.getFSSIndex() >= 0) {
							for (int sectIndex : sol.getRupSet().getSectionsIndicesForRup(rup.getFSSIndex()))
								maxMags[sectIndex] = Math.max(maxMags[sectIndex], mag);
						}
					}
					List<PlotCurveCharacterstics> sectChars = new ArrayList<>();
					List<Double> sectSortables = new ArrayList<>();
					for (int s=0; s<maxMags.length; s++) {
						if (maxMags[s] == 0d)
							sectChars.add(otherSectChar);
						else
							sectChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, magCPT.getColor(maxMags[s])));
						sectSortables.add(maxMags[s]);
					}
					mapMaker.plotScatters(rupLocs, chars);
//					mapMaker.plotSectScalars(maxMags, magCPT, "Magnitude");
					mapMaker.plotSectChars(sectChars, magCPT, "Magnitude", sectSortables);
					mapMaker.setSectNaNChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
					mapMaker.plot(outputDir, "scenario_catalog_"+numCatalogsWithTarget, " ");
					
					numCatalogsWithTarget++;
					
					mapMaker.clearScatters();
					mapMaker.clearSectScalars();
				}
				if (hasCorup)
					numCatalogsWithParentsCorup++;
				for (int p=0; p<scenarioParents.length; p++)
					if (hasParents[p])
						numCatalogsWithParents[p]++;
				for (int s=0; s<partics.length; s++)
					if (partics[s])
						sectPartics[s]++;
				if (numCatalogs % 1000 == 0)
					System.out.println("DONE catalog "+numCatalogs);
				for (int m=0; m<minMags.length; m++)
					if (hasMags[m])
						catMagCounts[m]++;
				
			}
		}
		
		DecimalFormat pDF = new DecimalFormat("0.00%");
		System.out.println(numCatalogsWithTarget+"/"+numCatalogs+" ("
				+pDF.format((double)numCatalogsWithTarget/(double)numCatalogs)+") had our exact scenario");
		System.out.println(numCatalogsWithParentsCorup+"/"+numCatalogs+" ("
				+pDF.format((double)numCatalogsWithParentsCorup/(double)numCatalogs)+") had a corupture");
		for (int p=0; p<parentNames.length; p++) {
			System.out.println(numCatalogsWithParents[p]+"/"+numCatalogs+" ("
					+pDF.format((double)numCatalogsWithParents[p]/(double)numCatalogs)+") had a "+parentNames[p]+" rupture");
		}
		CPT nuclCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-6d, -1d);
		nuclCPT.setLog10(true);
		nuclCPT.setNanColor(Color.WHITE);
		
		for (int m=0; m<minMags.length; m++) {
			nuclXYZ[m].scale(1d/(double)numCatalogs);
			
			mapMaker.plotXYZData(nuclXYZ[m], nuclCPT, "Expected num M>"+(int)minMags[m]);
			
			mapMaker.plot(outputDir, "nucl_m"+(int)minMags[m], " ");
		}
		
		mapMaker.clearXYZData();
		for (int s=0; s<sectPartics.length; s++)
			sectPartics[s] /= (double)numCatalogs;
		
		CPT particCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-5d, -1d);
		particCPT.setLog10(true);
		
		mapMaker.plotSectScalars(sectPartics, particCPT, "Fault participation probability");
		
		mapMaker.plot(outputDir, "sect_partic", " ");
		
		mapMaker.clearSectScalars();
		
		PlotCurveCharacterstics scenarioChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 6f, Colors.tab_green);
		PlotCurveCharacterstics[] parentChars = new PlotCurveCharacterstics[parentNames.length];
		for (int p=0; p<parentNames.length; p++)
			parentChars[p] = new PlotCurveCharacterstics(PlotLineType.SOLID, 6f, ColorUtils.TAB_10[p]);
		List<PlotCurveCharacterstics> sectChars = new ArrayList<>();
		HashSet<Integer> scenarioParentIDs = new HashSet<>();
		for (FaultSection sect : sol.getRupSet().getFaultSectionDataForRupture(inputScenarioID))
			scenarioParentIDs.add(sect.getParentSectionId());
		List<Double> sortables = new ArrayList<>();
		for (FaultSection sect : sects) {
			int parentID = sect.getParentSectionId();
			if (scenarioParentIDs.contains(parentID)) {
				sectChars.add(scenarioChar);
				sortables.add(100d);
			} else {
				boolean match = false;
				for (int p=0; p<parentNames.length; p++) {
					if (parentID == scenarioParents[p]) {
						match = true;
						sectChars.add(parentChars[p]);
						sortables.add(100d);
						break;
					}
				}
				if (!match) {
					sectChars.add(otherSectChar);
					sortables.add(0d);
				}
			}
		}
		mapMaker.plotSectChars(sectChars, null, null, sortables);
		
//		List<String> legendNames = new ArrayList<>();
//		List<PlotCurveCharacterstics> legendChars = new ArrayList<>();
//		legendNames.add("Scenario M7.8");
//		legendChars.add(scenarioChar);
//		for (int p=0; p<parentNames.length; p++) {
//			legendNames.add(parentNames[p]);
//			legendChars.add(parentChars[p]);
//		}
//		mapMaker.setCustomLegendItems(legendNames, legendChars);
		
		mapMaker.plot(outputDir, "sect_parents", " ");
		
		System.out.println("Scenario mag is "+(float)sol.getRupSet().getMagForRup(targetScenarioID));
		
		System.out.println("Overall catalog M probs");
		for (int m=0; m<minMags.length; m++)
			System.out.println("\tM>"+(int)minMags[m]+":\t"+pDF.format((double)catMagCounts[m]/(double)numCatalogs));
		
		// write cumulative timing funcs
		for (boolean isTarget : new boolean[] {true,false}) {
			EvenlyDiscretizedFunc func = isTarget ? targetWeeks : corupWeeks;
			EvenlyDiscretizedFunc cmlFunc = new EvenlyDiscretizedFunc(func.getMinX(), func.size(), func.getDelta());
			double sum = 0d;
			for (int i=0; i<func.size(); i++) {
				sum += func.getY(i);
				cmlFunc.set(i, sum);
			}
			
			String title = isTarget ? "Target rupture" : "Hollywood-Raymond corupture";
			
			PlotSpec plot = new PlotSpec(List.of(cmlFunc),
					List.of(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Colors.tab_blue)),
					title, "Weeks after scenario mainshock", "Cumulative number");
			
			HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
			
			gp.drawGraphPanel(plot, false, false, new Range(0d, 52d), null);
			
			String prefix = isTarget ? "timing_target" : "timing_corup";
			PlotUtils.writePlots(outputDir, prefix, gp, 800, 600, true, true, false);

			EvenlyDiscretizedFunc daysFunc = isTarget ? targetDays : corupDays;
			EvenlyDiscretizedFunc cmlDaysFunc = new EvenlyDiscretizedFunc(daysFunc.getMinX(), daysFunc.size(), daysFunc.getDelta());
			sum = 0d;
			for (int i=0; i<daysFunc.size(); i++) {
				sum += daysFunc.getY(i);
				cmlDaysFunc.set(i, sum);
			}
			
			System.out.println(title+" stats");
			System.out.println("\tCount: "+(int)sum);
			System.out.println("\t<1 day: "+(int)cmlDaysFunc.getY(0));
			System.out.println("\t<1 week: "+(int)cmlDaysFunc.getY(6));
			System.out.println("\t<1 month: "+(int)cmlDaysFunc.getY(29));
		}
	}

}
