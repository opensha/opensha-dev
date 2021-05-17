package scratch.kevin.simulators.ruptures.rotation;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipException;

import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.ngaw2.ASK_2014;
import org.opensha.sha.imr.attenRelImpl.ngaw2.BSSA_2014;
import org.opensha.sha.imr.attenRelImpl.ngaw2.CB_2014;
import org.opensha.sha.imr.attenRelImpl.ngaw2.CY_2014;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_Wrappers;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceJBParameter;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceRupParameter;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceX_Parameter;
import org.opensha.sha.simulators.RSQSimEvent;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.io.Files;

import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.rotation.RotatedRupVariabilityConfig.Quantity;
import scratch.kevin.simulators.ruptures.rotation.RotatedRupVariabilityConfig.RotationSpec;

public class RotatedRupVelModelScenarioComparePageGen {

	private Table<VelocityModel, Scenario, RSQSimRotatedRupVariabilityConfig> configTable;
	private Table<VelocityModel, Scenario, BBP_RotatedRupSimLoader> loaderTable;
	private Map<Scenario, Integer> eventCountMap;
	
	private Map<Integer, RSQSimEvent> eventsMap;
	private Map<Integer, EqkRupture> gmpeEventsMap;
	
	public RotatedRupVelModelScenarioComparePageGen(RSQSimCatalog catalog, Map<VelocityModel, File> vmDirsMap)
			throws ZipException, IOException {
		HashSet<Integer> eventIDsSet = new HashSet<>();

		configTable = HashBasedTable.create();
		loaderTable = HashBasedTable.create();
		eventCountMap = new HashMap<>();
		
		for (VelocityModel vm : vmDirsMap.keySet()) {
			File bbpDir = vmDirsMap.get(vm);
			File bbpZipFile = new File(bbpDir, "results_rotD.zip");
			Preconditions.checkState(bbpZipFile.exists());
			
			List<BBP_Site> bbpSites = BBP_Site.readFile(bbpDir);
			
			List<Site> sites = new ArrayList<>();
			for (BBP_Site site : bbpSites)
				sites.add(site.buildGMPE_Site(vm));
			for (Scenario scenario : Scenario.values()) {
				File rotConfFile = new File(bbpDir, "rotation_config_"+scenario.getPrefix()+".csv");
				if (rotConfFile.exists()) {
					RSQSimRotatedRupVariabilityConfig config = RSQSimRotatedRupVariabilityConfig.loadCSV(
							catalog, rotConfFile, null, sites);
					
					BBP_RotatedRupSimLoader bbpLoader = new BBP_RotatedRupSimLoader(bbpZipFile, bbpSites, scenario.getPrefix());
					
					List<Integer> events = config.getValues(Integer.class, Quantity.EVENT_ID);
					eventIDsSet.addAll(events);

					configTable.put(vm, scenario, config);
					loaderTable.put(vm, scenario, bbpLoader);
					
					if (eventCountMap.containsKey(scenario))
						eventCountMap.put(scenario, Integer.min(eventCountMap.get(scenario), events.size()));
					else
						eventCountMap.put(scenario, events.size());
				}
			}
		}
		
		eventsMap = RotatedRupVariabilityPageGen.loadEvents(catalog, eventIDsSet);
		
		System.out.println("Building GMPE ruptures...");
		gmpeEventsMap = new HashMap<>();
		for (RSQSimEvent event : eventsMap.values())
			gmpeEventsMap.put(event.getID(), catalog.getMappedSubSectRupture(event));
		System.out.println("DONE.");
	}
	
	public void generatePage(File outputDir, double[] periods, ScalarIMR[] refGMPEs) throws IOException {
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		List<String> lines = new LinkedList<>();
		
		lines.add("# Scenario-Based BBP Velocity Model Comparisons");
		lines.add("");
		lines.add("This page contains comparisons of BBP simulations for the same RSQSim catalog with different 1-dimensional layered"
			+ " velocity structures. Simulations are done for various scenario ruptures (we search the RSQSim catalog for matching"
			+ " events). This page uses results computed for the rupture rotation variability study.");
		lines.add("");
		lines.add("The following velocity models are shown:");
		lines.add("");
		for (VelocityModel vm : VelocityModel.values())
			if (loaderTable.containsRow(vm))
				lines.add("* "+vm+", Vs30="+(float)vm.getVs30()+" m/s");
		lines.add("");
		if (refGMPEs != null) {
			lines.add("We compare results with the following GMPE(s):");
			lines.add("");
			for (ScalarIMR gmpe : refGMPEs)
				lines.add("* "+gmpe.getName());
			lines.add("");
		}
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		Map<VelocityModel, Color> vmColors = new HashMap<>();
		Color[] possibleColors = { Color.BLUE.darker(), Color.RED.darker(), Color.GREEN.darker(), Color.CYAN.darker() };
		for (int i=0; i<VelocityModel.values().length; i++)
			vmColors.put(VelocityModel.values()[i], possibleColors[i % possibleColors.length]);
		
		String gfDir = System.getenv("BBP_INSTALL_GF");
		if (gfDir != null && gfDir.length() > 0) {
			System.out.println("Plotting velocity profiles! GF Dir: "+gfDir);
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			Map<VelocityModel, TableBuilder> vmTables = new HashMap<>();
			Map<VelocityModel, XY_DataSet> vmVsProfiles = new HashMap<>();
			
			for (VelocityModel vm : VelocityModel.values()) {
				if (!configTable.containsRow(vm))
					continue;
				File vmFile = vm.getFilePath(gfDir);
				XY_DataSet vsProfile = new DefaultXY_DataSet();
				vmVsProfiles.put(vm, vsProfile);
				XY_DataSet vpProfile = new DefaultXY_DataSet();
				
				TableBuilder table = MarkdownUtils.tableBuilder();
				table.addLine("Start Depth (km)", "End Depth (km)", "Thickness (km)", "Vs (km/s)", "Vp (km/s)");
				
				double zStart = 0;
				for (String line : Files.readLines(vmFile, Charset.defaultCharset())) {
					line = line.trim();
					while (line.contains("  "))
						line = line.replaceAll("  ", " ");
					line.replaceAll("\t", " ");
					String[] split = line.split(" ");
					if (split.length > 2) {
						double thickness = Double.parseDouble(split[0]);
						double vp = Double.parseDouble(split[1]);
						double vs = Double.parseDouble(split[2]);
						
						vpProfile.set(vp, zStart);
						vpProfile.set(vp, zStart+thickness);
						vsProfile.set(vs, zStart);
						vsProfile.set(vs, zStart+thickness);
						
						table.addLine((float)zStart, (float)(zStart+thickness), (float)thickness, (float)vs, (float)vp);
						
						zStart += thickness;
					}
				}
				
				vmTables.put(vm, table);
				
				Color color = vmColors.get(vm);
				
				vsProfile.setName(vm+" Vs");
				funcs.add(vsProfile);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, color));
				
				vpProfile.setName(vm+" Vp");
				funcs.add(vpProfile);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1.5f, color));
			}
			
			PlotSpec spec = new PlotSpec(funcs, chars, "Full Velocity Profiles", "Velocity (km/s)", "Depth (km)");
			spec.setLegendVisible(true);
			
			PlotPreferences plotPrefs = PlotPreferences.getDefault();
			plotPrefs.setTickLabelFontSize(18);
			plotPrefs.setAxisLabelFontSize(20);
			plotPrefs.setPlotLabelFontSize(21);
			plotPrefs.setBackgroundColor(Color.WHITE);

			HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
			gp.setyAxisInverted(true);
			
			Range xRange = new Range(0, 8);
			Range yRange = new Range(0, 10);

			gp.drawGraphPanel(spec, false, false, xRange, yRange);
			gp.getChartPanel().setSize(800, 600);
			File pngFile = new File(resourcesDir, "vm_profiles.png");
			gp.saveAsPNG(pngFile.getAbsolutePath());
			
			lines.add("## Velocity Model Profiles");
			lines.add(topLink); lines.add("");

			lines.add("![VM Profiles](resources/"+pngFile.getName()+")");
			lines.add("");
			
			xRange = new Range(0, 4);
			yRange = new Range(0, 1);

			spec.setTitle("Upper Velocity Profiles");
			gp.drawGraphPanel(spec, false, false, xRange, yRange);
			gp.getChartPanel().setSize(800, 600);
			pngFile = new File(resourcesDir, "vm_profiles_upper.png");
			gp.saveAsPNG(pngFile.getAbsolutePath());
			
			lines.add("![Upper VM Profiles](resources/"+pngFile.getName()+")");
			lines.add("");
			
			for (VelocityModel vm : VelocityModel.values()) {
				if (!configTable.containsRow(vm))
					continue;
				
				lines.add("### "+vm+" Profile");
				lines.add(topLink); lines.add("");
				
				lines.addAll(vmTables.get(vm).build());
				lines.add("");
			}
			
			TableBuilder zTable = MarkdownUtils.tableBuilder();
			
			zTable.addLine("VM Name", "VM Z1.0 (km)", "ASK 2014 Z1.0 (km)", "BSSA 2014 Z1.0 (km)", "CY 2014 Z1.0 (km)",
					"VM Z2.5 (km)", "CB 2014 Z2.5 (km)");
			
			for (VelocityModel vm : VelocityModel.values()) {
				if (!configTable.containsRow(vm))
					continue;
				
				zTable.initNewLine();
				zTable.addColumn("**"+vm+"**");
				
				XY_DataSet vsProfile = vmVsProfiles.get(vm);
				double z10 = Double.NaN;
				double z25 = Double.NaN;
				for (Point2D pt : vsProfile) {
					if (Double.isNaN(z10) && pt.getX() >= 1d)
						z10 = pt.getY();
					if (Double.isNaN(z25) && pt.getX() >= 2.5d)
						z25 = pt.getY();
				}
				
				zTable.addColumn("**"+(float)z10+"**");
				double vs30 = vm.getVs30();
				zTable.addColumn((float)ASK_2014.calcZ1ref(vs30));
				zTable.addColumn((float)BSSA_2014.calcZ1ref(vs30));
				zTable.addColumn((float)CY_2014.calcZ1ref(vs30));
				zTable.addColumn("**"+(float)z25+"**");
				zTable.addColumn((float)CB_2014.calcZ25ref(vs30));
				
				zTable.finalizeLine();
			}
			
			lines.add("### Z1.0 and Z2.5 Comparisons");
			lines.add(topLink); lines.add("");
			
			lines.addAll(zTable.build());
			lines.add("");
		}
		
		for (Scenario scenario : Scenario.values()) {
			if (!configTable.containsColumn(scenario))
				continue;
			
			lines.add("## "+scenario.getName());
			lines.add(topLink); lines.add("");
			
			lines.add("Median spectra shown across "+eventCountMap.get(scenario)
				+" RSQSim events which match the following scenario criteria:");
			lines.add("");
			for (String criterion : scenario.getMatchCriteria())
				lines.add("* "+criterion);
			lines.add("");
			
			HashSet<Float> distanceSet = new HashSet<>();
			for (VelocityModel vm : configTable.rowKeySet())
				distanceSet.addAll(configTable.get(vm, scenario).getValues(Float.class, Quantity.DISTANCE));
			
			List<Float> distances = new ArrayList<>(distanceSet);
			Collections.sort(distances);
			
			for (Float distance : distances) {
				lines.add("### "+scenario.getName()+", "+distance+" km");
				lines.add(topLink); lines.add("");
				
				List<DiscretizedFunc> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				for (VelocityModel vm : VelocityModel.values()) {
					if (!configTable.contains(vm, scenario))
						continue;
					
					Color color = vmColors.get(vm);
					
					DiscretizedFunc simFunc = calcSim(scenario, vm, distance);
					if (simFunc == null)
						continue;
					simFunc.setName("Simulated, Vs30="+(float)vm.getVs30());
					funcs.add(simFunc);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, color));
					
					if (refGMPEs != null) {
						List<DiscretizedFunc> gmpeFuncs = new ArrayList<>();
						for (ScalarIMR gmpe : refGMPEs)
							gmpeFuncs.add(calcGMPE(scenario, vm, gmpe, distance, periods));
						DiscretizedFunc gmpeFunc;
						if (gmpeFuncs.size() == 1)
							gmpeFunc = gmpeFuncs.get(0);
						else
							gmpeFunc = meanSpectrum(gmpeFuncs);
						
						gmpeFunc.setName("GMPE");
						funcs.add(gmpeFunc);
						chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1.5f, color));
					}
				}
				
				PlotSpec spec = new PlotSpec(funcs, chars, scenario.getName()+", "+distance+" km",
						"Period (s)", "Median Ground Motion (g)");
				spec.setLegendVisible(true);
				
				PlotPreferences plotPrefs = PlotPreferences.getDefault();
				plotPrefs.setTickLabelFontSize(18);
				plotPrefs.setAxisLabelFontSize(20);
				plotPrefs.setPlotLabelFontSize(21);
				plotPrefs.setBackgroundColor(Color.WHITE);

				HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
				
				double minY = Double.POSITIVE_INFINITY;
				double maxY = Double.NEGATIVE_INFINITY;
				for (DiscretizedFunc func : funcs) {
					minY = Double.min(minY, func.getMinY());
					maxY = Double.max(maxY, func.getMaxY());
				}
				
				Range xRange = new Range(StatUtils.min(periods), StatUtils.max(periods));
				Range yRange = new Range(Math.pow(10, Math.floor(Math.log10(minY))), Math.pow(10, Math.ceil(Math.log10(maxY))));

				gp.drawGraphPanel(spec, false, true, xRange, yRange);
				gp.getChartPanel().setSize(800, 600);
				File pngFile = new File(resourcesDir, scenario.getPrefix()+"_dist"+distance+".png");
				gp.saveAsPNG(pngFile.getAbsolutePath());
				
				lines.add("!["+distance+" km](resources/"+pngFile.getName()+")");
				lines.add("");
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 2));
		lines.add(tocIndex, "## Table Of Contents");
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private DiscretizedFunc calcSim(Scenario scenario, VelocityModel vm, Float distance) throws IOException {
		RSQSimRotatedRupVariabilityConfig config = configTable.get(vm, scenario);
		BBP_RotatedRupSimLoader loader = loaderTable.get(vm, scenario);
		
		List<DiscretizedFunc> spectra = new ArrayList<>();
		for (RotationSpec rotation : config.getRotationsForQuantities(Quantity.DISTANCE, distance))
			spectra.add(loader.getRotD50(rotation));
		if (spectra.isEmpty())
			return null;
		
		return medianSpectrum(spectra);
	}
	
	private static DiscretizedFunc medianSpectrum(List<DiscretizedFunc> spectra) {
		DiscretizedFunc meanSpectrum = new ArbitrarilyDiscretizedFunc();
		DiscretizedFunc refSpectrum = spectra.get(0);
		for (int p=0; p<refSpectrum.size(); p++) {
			double[] vals = new double[spectra.size()];
			for (int i=0; i<spectra.size(); i++) {
				Preconditions.checkState(spectra.get(i).size() == refSpectrum.size());
				vals[i] = spectra.get(i).getY(p);
			}
			meanSpectrum.set(refSpectrum.getX(p), DataUtils.median(vals));
		}
		return meanSpectrum;
	}
	
	private static DiscretizedFunc meanSpectrum(List<DiscretizedFunc> spectra) {
		DiscretizedFunc meanSpectrum = new ArbitrarilyDiscretizedFunc();
		DiscretizedFunc refSpectrum = spectra.get(0);
		for (int p=0; p<refSpectrum.size(); p++) {
			double[] vals = new double[spectra.size()];
			for (int i=0; i<spectra.size(); i++) {
				Preconditions.checkState(spectra.get(i).size() == refSpectrum.size());
				vals[i] = spectra.get(i).getY(p);
			}
			meanSpectrum.set(refSpectrum.getX(p), StatUtils.mean(vals));
		}
		return meanSpectrum;
	}
	
	private DiscretizedFunc calcGMPE(Scenario scenario, VelocityModel vm, ScalarIMR gmpe, Float distance, double[] periods) {
		gmpe.setParamDefaults();
		
		RSQSimRotatedRupVariabilityConfig config = configTable.get(vm, scenario);
		List<Integer> eventIDs = config.getValues(Integer.class, Quantity.EVENT_ID);
		
		List<double[]> results = new ArrayList<>();
		
		gmpe.setSite(config.getValues(Site.class, Quantity.SITE).get(0));
		
		DistanceJBParameter rJBParam = (DistanceJBParameter) gmpe.getParameter(DistanceJBParameter.NAME);
		DistanceRupParameter rRupParam = (DistanceRupParameter) gmpe.getParameter(DistanceRupParameter.NAME);
		DistanceX_Parameter rXParam = (DistanceX_Parameter) gmpe.getParameter(DistanceX_Parameter.NAME);
		
		List<Float> sourceAzimuths = config.getValues(Float.class, Quantity.SOURCE_AZIMUTH);
		int numSourceAz = sourceAzimuths.size();
		
		gmpe.setIntensityMeasure(SA_Param.NAME);
		SA_Param saParam = (SA_Param) gmpe.getIntensityMeasure();
		
		for (int eventID : eventIDs) {
			EqkRupture rup = gmpeEventsMap.get(eventID);
			gmpe.setEqkRupture(rup);
			
			double zTOR = rup.getRuptureSurface().getAveRupTopDepth();
			double rRup, rJB;
			if (BBP_PartBValidationConfig.DIST_JB) {
				rJB = distance;
				rRup = zTOR == 0d ? rJB : Math.sqrt(zTOR*zTOR + rJB*rJB);
			} else {
				rRup = distance;
				rJB = zTOR == 0d ? rRup : Math.sqrt(rRup*rRup - zTOR*zTOR);
			}
			
			// override distances
			rJBParam.setValueIgnoreWarning(rJB);
			rRupParam.setValueIgnoreWarning(rRup);
			
			for (int i=0; i<numSourceAz; i++) {
				double sourceAz = sourceAzimuths.get(i).doubleValue();
				double rX = rJB * Math.sin(Math.toRadians(sourceAz));
				
				rXParam.setValueIgnoreWarning(rX);
				
				double[] vals = new double[periods.length];
				for (int p=0; p<periods.length; p++) {
					SA_Param.setPeriodInSA_Param(saParam, periods[p]);
					vals[p] = gmpe.getMean();
				}
				results.add(vals);
			}
		}
		
		DiscretizedFunc ret = new ArbitrarilyDiscretizedFunc(gmpe.getShortName());
		
		for (int p=0; p<periods.length; p++) {
			double[] vals = new double[results.size()];
			for (int i=0; i<vals.length; i++)
				vals[i] = results.get(i)[p];
			ret.set(periods[p], Math.exp(StatUtils.mean(vals)));
		}
		
		return ret;
	}

	public static void main(String[] args) throws IOException, DocumentException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		File outputDir = new File("/home/kevin/markdown/rsqsim-analysis/catalogs");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");

		Map<VelocityModel, File> vmDirsMap = new HashMap<>();
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
//		vmDirsMap.put(VelocityModel.LA_BASIN_500, new File(bbpParallelDir,
//				"2019_02_27-rundir2585_1myrs-rotatedRups-6scenarios-3dists-18srcAz-1siteSrcAz-400rups-skipYears5000-vmLA_BASIN_500-noHF-1site"));
//		vmDirsMap.put(VelocityModel.LA_BASIN_863, new File(bbpParallelDir,
//				"2019_02_21-rundir2585_1myrs-rotatedRups-6scenarios-3dists-18srcAz-1siteSrcAz-400rups-skipYears5000-vmLA_BASIN_863-noHF-1site"));
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2740.instance(baseDir);
		vmDirsMap.put(VelocityModel.LA_BASIN_500, new File(bbpParallelDir,
				"2019_02_27-rundir2740-rotatedRups-6scenarios-3dists-18srcAz-1siteSrcAz-400rups-skipYears5000-vmLA_BASIN_500-noHF-1site"));
		vmDirsMap.put(VelocityModel.LA_BASIN_863, new File(bbpParallelDir,
				"2019_02_21-rundir2740-rotatedRups-6scenarios-3dists-18srcAz-1siteSrcAz-400rups-skipYears5000-vmLA_BASIN_863-noHF-1site"));
		
		double[] periods = {1d, 2d, 3d, 4d, 5d, 7.5, 10d};
		
		ScalarIMR[] refGMPEs = { new NGAW2_Wrappers.ASK_2014_Wrapper(), new NGAW2_Wrappers.BSSA_2014_Wrapper(),
				new NGAW2_Wrappers.CB_2014_Wrapper(), new NGAW2_Wrappers.CY_2014_Wrapper()};
		
		RotatedRupVelModelScenarioComparePageGen pageGen = new RotatedRupVelModelScenarioComparePageGen(catalog, vmDirsMap);
		
		File catalogOutputDir = new File(outputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		File compareDir = new File(catalogOutputDir, "bbp_vm_rot_rup_compare");
		Preconditions.checkState(compareDir.exists() || compareDir.mkdir());
		
		pageGen.generatePage(compareDir, periods, refGMPEs);
		
		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(outputDir);
	}

}
