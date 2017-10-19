package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.kevin.MarkdownUtils;
import scratch.kevin.MarkdownUtils.TableBuilder;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_SimZipLoader;
import scratch.kevin.bbp.BBP_SimZipLoader.BBP_RupGenSimZipLoader;
import scratch.kevin.bbp.BBP_SimZipLoader.BBP_ShakeMapSimZipLoader;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.SeismogramPlotter;
import scratch.kevin.bbp.ShakemapPlotter;
import scratch.kevin.bbp.SpectraPlotter;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.plots.RuptureVelocityPlot;

class RupSpectraPageGen {
	
	private RSQSimCatalog catalog;
	private File outputDir;
	private VelocityModel vm;
	
	private ScalarIMR[] spectraGMPEs;
	private File refBBPDir;
	private String refName;
	private EqkRupture gmpeRup;
	
	private BBP_ShakeMapSimZipLoader shakemap;
	private GriddedRegion shakemapRegion;
	private List<BBP_Site> shakemapSites;
	private ScalarIMR shakemapGMPE;

	public RupSpectraPageGen(RSQSimCatalog catalog, File outputDir, VelocityModel vm) {
		this.catalog = catalog;
		this.outputDir = outputDir;
		this.vm = vm;
	}
	
	public void setRefBBP(File refBBPDir, String refName, ScalarIMR[] spectraGMPEs, EqkRupture gmpeRup) {
		this.refBBPDir = refBBPDir;
		this.refName = refName;
		this.spectraGMPEs = spectraGMPEs;
		this.gmpeRup = gmpeRup;
	}
	
	public void setShakeMap(BBP_ShakeMapSimZipLoader shakemap, GriddedRegion shakemapRegion,
			List<BBP_Site> shakemapSites, ScalarIMR shakemapGMPE) {
		this.shakemap = shakemap;
		this.shakemapRegion = shakemapRegion;
		this.shakemapSites = shakemapSites;
		this.shakemapGMPE = shakemapGMPE;
	}
	
	private static DecimalFormat twoDigitsDF = new DecimalFormat("0.00");
	
	public void generatePage(RSQSimEvent event, File eventBBPDir, List<BBP_Site> sites,
			boolean rebuildGIF, boolean rebuildMaps) throws IOException, GMT_MapException {
		int eventID = event.getID();
		File eventDir = new File(outputDir, "event_"+eventID);
		Preconditions.checkState(eventDir.exists() || eventDir.mkdir());
		File resourcesDir = new File(eventDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		LinkedList<String> lines = new LinkedList<>();
		
		// header
		lines.add("# Event "+eventID+", M"+twoDigitsDF.format(event.getMagnitude()));
		lines.add("");
		lines.add("[Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
		lines.add("");
		BBP_RupGenSimZipLoader refLoader = null;
		if (refBBPDir != null) {
			File zipFile = new File(refBBPDir, "results.zip");
			refLoader = new BBP_RupGenSimZipLoader(zipFile, sites);
			lines.add("**"+refName+" Simulations: "+refLoader.getNumSims()+"**");
			lines.add("");
		}
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Rupture Plots");
		lines.add(topLink); lines.add("");
		BBP_SourceFile bbpSource = null;
		Location[] bbpSourceRect = null;
		Location bbpSourceHypo = null;
		lines.add("**Legend**");
		lines.add("* Colored, Filled Triangles: RSQSim Elements");
		lines.add("* Red Star: RSQSim Hypocenter");
		if (refBBPDir != null) {
			lines.add("* Dark Green Solid Outline: BBP Equivalent Planar Surface");
			lines.add("* Green Star: BBP Equivalent Hypocenter");
			bbpSource = loadSourceFile(refBBPDir);
			bbpSourceRect = bbpSource.getSurface().getRectangle();
			bbpSourceHypo = bbpSource.getHypoLoc();
		}
		RuptureSurface gmpeSurf = null;
		if (gmpeRup != null) {
			gmpeSurf = gmpeRup.getRuptureSurface();
			lines.add("* Dark Gray Dashed Outline: GMPE Surface");
		}
		lines.add("");
		RSQSimEventSlipTimeFunc func = catalog.getSlipTimeFunc(event);
		// slip time plot
		String rupPlotPrefix = "rupture_plot_"+eventID;
		RupturePlotGenerator.writeSlipPlot(event, func, resourcesDir, rupPlotPrefix, bbpSourceRect, bbpSourceHypo, gmpeSurf);
		File rupPlot = new File(resourcesDir, rupPlotPrefix+".png");
		Preconditions.checkState(rupPlot.exists());
		lines.add("### Slip/Time Plot");
		lines.add(topLink); lines.add("");
		lines.add("![Slip/Time Plot]("+resourcesDir.getName()+"/"+rupPlot.getName()+")");
		
		// animation
		File rupAnim = new File(resourcesDir, rupPlotPrefix+".gif");
		if (rebuildGIF || !rupAnim.exists())
			RupturePlotGenerator.writeSlipAnimation(event, func, rupAnim, 5, bbpSourceRect[0], bbpSourceRect[1]);
		Preconditions.checkState(rupAnim.exists());
		lines.add("### Slip/Vel Animation");
		lines.add(topLink); lines.add("");
		lines.add("[Click here to view Slip/Velocity Animation]("+resourcesDir.getName()+"/"+rupAnim.getName()+")");
		
		// map view plot
		String mapRupPlotPrefix = "rupture_map_plot_"+eventID;
		RupturePlotGenerator.writeMapPlot(catalog.getElements(), event, func, resourcesDir, mapRupPlotPrefix,
				bbpSourceRect, bbpSourceHypo, gmpeSurf);
		File rupMapPlot = new File(resourcesDir, mapRupPlotPrefix+".png");
		Preconditions.checkState(rupMapPlot.exists());
		lines.add("### Map Plot");
		lines.add(topLink); lines.add("");
		lines.add("![Map Plot]("+resourcesDir.getName()+"/"+rupMapPlot.getName()+")");
		
		// rupture velocity plot
		lines.add("### Rupture Velocity Plot");
		RuptureVelocityPlot rupVelPlot = new RuptureVelocityPlot(catalog.getElements(), 0d);
		rupVelPlot.initialize(catalog.getName(), resourcesDir, "rupture_velocity");
		rupVelPlot.processEvent(event);
		rupVelPlot.finalizePlot();
		lines.add(topLink); lines.add("");
		lines.add("![Rupture Velocity Plot]("+resourcesDir.getName()+"/rupture_velocity_scatter_dist.png)");
		
		lines.add("");
		lines.add("## Spectra Plots");
		lines.add(topLink); lines.add("");
		
		for (BBP_Site site : sites) {
			String siteName = RSQSimBBP_Config.siteCleanName(site);
			System.out.println("Site: "+siteName);
			lines.add("### Site "+siteName);
			lines.add(topLink); lines.add("");
			Location loc = site.getLoc();
			lines.add("*Location: "+(float)loc.getLatitude()+", "+(float)loc.getLongitude()+"*");
			
			// calculate distances
			TableBuilder table = MarkdownUtils.tableBuilder();
			List<String> distHeader = new ArrayList<>();
			distHeader.add("Distance");
			distHeader.add("Actual RSQSim Surface");
			if (bbpSource != null)
				distHeader.add("BBP Equivalent Planar Surface");
			if (gmpeSurf != null)
				distHeader.add("GMPE Surface");
			table.addLine(distHeader);
			double horzDist = Double.POSITIVE_INFINITY;
			double linearDist = Double.POSITIVE_INFINITY;
			for (SimulatorElement e : event.getAllElements()) {
				for (Location v : e.getVertices()) {
					horzDist = Math.min(horzDist, LocationUtils.horzDistanceFast(loc, v));
					linearDist = Math.min(linearDist, LocationUtils.linearDistanceFast(loc, v));
				}
			}
			List<String> horzLine = new ArrayList<>();
			horzLine.add("Horizontal");
			horzLine.add(twoDigitsDF.format(horzDist)+" km");
			if (bbpSource != null)
				horzLine.add(twoDigitsDF.format(bbpSource.getSurface().getHorizontalDistance(loc))+" km");
			if (gmpeSurf != null)
				horzLine.add(twoDigitsDF.format(gmpeSurf.getDistanceJB(loc))+" km");
			table.addLine(horzLine);
			List<String> linearLine = new ArrayList<>();
			linearLine.add("3-D");
			linearLine.add(twoDigitsDF.format(linearDist)+" km");
			if (bbpSource != null)
				linearLine.add(twoDigitsDF.format(bbpSource.getSurface().getLinearDistance(loc))+" km");
			if (gmpeSurf != null)
				linearLine.add(twoDigitsDF.format(gmpeSurf.getDistanceRup(loc))+" km");
			table.addLine(linearLine);
			lines.add("");
			lines.addAll(table.build());
			lines.add("");
			lines.add("*NOTE: RSQSim ruptures sometimes have a few co-rupturing elements on faults some distance from "
					+ "the main rupture. This may cause discrepancies in the table above, consult rupture map plot.*");

			
			File fasFile = SpectraPlotter.findFASFile(eventBBPDir, site.getName());
			if (fasFile.exists()) {
				lines.add("#### "+siteName+" Fourier Amplitude Spectra");
				lines.add(topLink); lines.add("");
				String title = "Event "+eventID+" "+siteName+" FAS Spectra";
				String prefix = "fas_spectra_"+site.getName();
				if (refLoader != null) {
					List<DiscretizedFunc> refSpectra = new ArrayList<>();
					int numSims = refLoader.getNumSims();
					for (int i=0; i<numSims; i++)
						refSpectra.add(refLoader.readFAS(site, i));
					DiscretizedFunc dataSpectra = SpectraPlotter.loadFAS(fasFile);
					SpectraPlotter.plotMultiFAS(refSpectra, refName, dataSpectra, "RSQSim", title, resourcesDir, prefix);
				} else {
					SpectraPlotter.plotFAS(fasFile, resourcesDir, prefix);
				}
				lines.add("!["+siteName+" FAS Plot]("+resourcesDir.getName()+"/"+prefix+".png)");
				
			}
			
			File rotD50File = SpectraPlotter.findRotD50File(eventBBPDir, site.getName());
			if (rotD50File.exists()) {
				UncertainArbDiscDataset[] gmpeSpectra = null;
				if (spectraGMPEs != null && spectraGMPEs.length > 0) {
					System.out.print("Calculating GMPEs...");
					gmpeSpectra = new UncertainArbDiscDataset[spectraGMPEs.length];
					for (int i=0; i<gmpeSpectra.length; i++)
						gmpeSpectra[i] = SpectraPlotter.calcGMPE_RotD50(gmpeRup, site, vm, spectraGMPEs[i]);
					System.out.println("DONE.");
				}
				
				lines.add("#### "+siteName+" RotD50 Spectra");
				lines.add(topLink); lines.add("");
				String title = "Event "+eventID+" "+siteName+" RotD50 Spectra";
				String prefix = "rotd50_spectra_"+site.getName();
				if (refLoader != null) {
					List<DiscretizedFunc> refSpectra = new ArrayList<>();
					int numSims = refLoader.getNumSims();
					for (int i=0; i<numSims; i++)
						refSpectra.add(refLoader.readRotD50(site, i));
					DiscretizedFunc dataSpectra = SpectraPlotter.loadRotD50(rotD50File);
					SpectraPlotter.plotMultiRotD50(refSpectra, refName, dataSpectra, "RSQSim", title,
							resourcesDir, prefix, gmpeSpectra);
				} else {
					SpectraPlotter.plotRotD50(rotD50File, resourcesDir, prefix, gmpeSpectra);
				}
				lines.add("!["+siteName+" RotD50 Plot]("+resourcesDir.getName()+"/"+prefix+".png)");
			}
			
			File accelFile = SeismogramPlotter.findBBP_SeisFile(eventBBPDir, site.getName(), true);
			int numComps = 3;
			DiscretizedFunc[] accelSeis = SeismogramPlotter.loadBBP_Seis(accelFile);
			List<DiscretizedFunc[]> accelComps = new ArrayList<>();
			if (refLoader != null)
				for (int i=0; i<numComps; i++)
					accelComps.add(refLoader.readAccelSeis(site, i));
			String accelPrefix = "seis_accel_"+site.getName();
			SeismogramPlotter.plotSeismograms(accelSeis, "Event "+eventID+", "+siteName, true,
					resourcesDir, accelPrefix, false, accelComps);
			lines.add("#### "+siteName+" Acceleration Seismograms");
			lines.add(topLink); lines.add("");
			if (!accelComps.isEmpty()) {
				lines.add("RSQSim ruptures in blue. Gray seismograms are "+refName+" comparisons.");
				lines.add("");
			}
			lines.add("!["+siteName+" Acceleration Seismograms]("+resourcesDir.getName()+"/"+accelPrefix+".png)");
			
			File velFile = SeismogramPlotter.findBBP_SeisFile(eventBBPDir, site.getName(), false);
			DiscretizedFunc[] velSeis = SeismogramPlotter.loadBBP_Seis(velFile);
			List<DiscretizedFunc[]> velComps = new ArrayList<>();
			if (refLoader != null)
				for (int i=0; i<numComps; i++)
					velComps.add(refLoader.readVelSeis(site, i));
			String velPrefix = "seis_vel_"+site.getName();
			SeismogramPlotter.plotSeismograms(velSeis, "Event "+eventID+", "+siteName, false,
					resourcesDir, velPrefix, false, velComps);
			lines.add("#### "+siteName+" Velocity Seismograms");
			lines.add(topLink); lines.add("");
			if (!velComps.isEmpty()) {
				lines.add("RSQSim ruptures in blue. Gray seismograms are "+refName+" comparisons.");
				lines.add("");
			}
			lines.add("!["+siteName+" Velocity Seismograms]("+resourcesDir.getName()+"/"+velPrefix+".png)");
		}
		
		if (shakemap != null) {
			lines.add("");
			lines.add("## ShakeMaps");
			lines.add(topLink); lines.add("");
			
			double[] periods = { 1d, 2d, 3d, 4d, 5d, 7.5, 10d };
			String prefix = "shakemap";
			File[] mapFiles = new File[periods.length];
			File[] gmpeMapFiles = null;
			File[] ratioMapFiles = null;
			if (shakemapGMPE != null) {
				gmpeMapFiles = new File[periods.length];
				ratioMapFiles = new File[periods.length];
			}
			for (int i=0; i<periods.length; i++) {
				double p = periods[i];
				String name;
				if (p == Math.round(p))
					name = (int)p+"s";
				else
					name = (float)p+"s";
				mapFiles[i] = new File(resourcesDir, prefix+"_"+name+".png");
				if (shakemapGMPE != null) {
					gmpeMapFiles[i] = new File(resourcesDir,
							prefix+"_"+name+"_"+shakemapGMPE.getShortName()+".png");
					ratioMapFiles[i] = new File(resourcesDir,
							prefix+"_"+name+"_"+shakemapGMPE.getShortName()+"_ratio.png");
				}
			}
			
			double[] plotPeriods = periods;
			if (!rebuildMaps) {
				// see if we can skip any...
				List<Double> retainedPeriods = new ArrayList<>();
				for (int i=0; i<periods.length; i++) {
					boolean skip = mapFiles[i].exists();
					if (shakemapGMPE != null)
						skip = skip && gmpeMapFiles[i].exists() && ratioMapFiles[i].exists();
					if (!skip)
						retainedPeriods.add(periods[i]);
				}
				plotPeriods = Doubles.toArray(retainedPeriods);
			}
			
			if (periods.length > 0) {
				String label = "Event "+eventID;
				ShakemapPlotter.plotShakemaps(shakemap, shakemapRegion, shakemapSites, label, resourcesDir,
						prefix, true, shakemapGMPE, gmpeRup, vm, plotPeriods);
			} else {
				System.out.println("Skipping all maps!");
			}
			
			TableBuilder table = MarkdownUtils.tableBuilder();
			table.initNewLine().addColumn("SA Period").addColumn("RSQSim");
			if (shakemapGMPE != null)
				table.addColumn(shakemapGMPE.getShortName()).addColumn("Ratio");
			table.finalizeLine();
			for (int i=0; i<periods.length; i++) {
				table.initNewLine().addColumn("**"+(float)periods[i]+" s**");
				table.addColumn("![RSQSim Map]("+resourcesDir.getName()+"/"+mapFiles[i].getName()+")");
				if (shakemapGMPE != null) {
					table.addColumn("!["+shakemapGMPE.getShortName()+" Map]("
							+resourcesDir.getName()+"/"+gmpeMapFiles[i].getName()+")");
					table.addColumn("!["+shakemapGMPE.getShortName()+" Ratio]("
							+resourcesDir.getName()+"/"+ratioMapFiles[i].getName()+")");
				}
				table.finalizeLine();
			}
			
			lines.add("");
			lines.addAll(table.build());
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, eventDir);
	}
	
	private BBP_SourceFile loadSourceFile(File bbpDir) throws IOException {
		for (File file : bbpDir.listFiles())
			if (file.getName().toLowerCase().endsWith(".src"))
				return BBP_SourceFile.readFile(file);
		throw new FileNotFoundException("No .src files found in "+bbpDir.getAbsolutePath());
	}
	
	public static void main(String[] args) throws IOException, DocumentException, GMT_MapException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		File outputDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");
		String refName = "Graves & Pitarka (2015)";
		boolean rebuildGIF = false;
		boolean rebuildMaps = false;
		FaultBasedMapGen.LOCAL_MAPGEN = true;
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2194_LONG.instance(baseDir);
//		int eventID = 136704;
		
//		RSQSimCatalog catalog = Catalogs.JG_UCERF3_millionElement.instance(baseDir);
//		int eventID = 4099020;
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2273.instance(baseDir);
//		int eventID = 412778;
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2310.instance(baseDir);
//		int eventID = 3802809;
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2320.instance(baseDir);
		int eventID = 6195527;
		
		File eventBBPDir = new File(new File(catalog.getCatalogDir(), RSQSimBBP_Config.EVENT_SRF_DIR_NAME),
				"event_"+eventID+"_"+RSQSimBBP_Config.SRF_DT+"s_"+RSQSimBBP_Config.SRF_INTERP_MODE.name()+"_bbp");
		
		// find BBP parallel dir
		File refBBPDir = null;
		File[] refDirs = bbpParallelDir.listFiles();
		Arrays.sort(refDirs, new FileNameComparator());
		for (File dir : refDirs) {
			String name = dir.getName();
			if (dir.isDirectory() && name.contains(catalog.getCatalogDir().getName()) && name.contains(eventID+"")
					&& name.contains("-gp-") && !name.contains("shakemap")) {
				File zipFile = new File(dir, "results.zip");
				if (zipFile.exists())
					refBBPDir = dir;
			}
		}
		if (refBBPDir != null)
			System.out.println("Located ref BBP dir: "+refBBPDir.getAbsolutePath());
		
		File refShakeMapZip = null;
		for (File dir : refDirs) {
			String name = dir.getName();
			if (dir.isDirectory() && name.contains(catalog.getCatalogDir().getName()) && name.contains(eventID+"")
					&& name.contains("shakemap")) {
				File shakemapFile = new File(dir, "results_rotD.zip");
				if (!shakemapFile.exists())
					shakemapFile = new File(dir, "results.zip");
				if (shakemapFile.exists())
					refShakeMapZip = shakemapFile;
			}
		}
		if (refShakeMapZip != null)
			System.out.println("Located ref ShakeMap zip: "+refShakeMapZip.getAbsolutePath());
		
		RSQSimEvent event = catalog.loader().byID(eventID);
		
		List<BBP_Site> sites = null;
		if (refBBPDir != null)
			sites = BBP_Site.readFile(refBBPDir);
		
		boolean runBBP = !eventBBPDir.exists();
		if (!runBBP) {
			if (sites == null) {
				sites = BBP_Site.readFile(eventBBPDir);
				if (sites.isEmpty()) {
					sites = null;
					runBBP = true;
				}
			}
		}
		if (runBBP) {
			System.out.println("Need to run the BBP...");
			RSQSimBBP_Config.runBBP(catalog, event, sites);
		}
		
		ScalarIMR[] gmpes = { AttenRelRef.ASK_2014.instance(null), AttenRelRef.BSSA_2014.instance(null),
				AttenRelRef.CB_2014.instance(null), AttenRelRef.CY_2014.instance(null) };
		for (ScalarIMR gmpe : gmpes)
			gmpe.setParamDefaults();
		ScalarIMR shakemapGMPE = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);
		shakemapGMPE.setParamDefaults();
		VelocityModel vm = VelocityModel.LA_BASIN;
		
		File catalogOutputDir = new File(outputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		RupSpectraPageGen gen = new RupSpectraPageGen(catalog, catalogOutputDir, vm);
		
		EqkRupture gmpeRup = null;
		
		if ((refBBPDir != null && gmpes != null && gmpes.length > 0)
				|| (refShakeMapZip != null && shakemapGMPE != null))
			gmpeRup = catalog.getGMPE_Rupture(event, RSQSimBBP_Config.MIN_SUB_SECT_FRACT);
		
		if (refBBPDir != null)
			gen.setRefBBP(refBBPDir, refName, gmpes, gmpeRup);
		
		if (refShakeMapZip != null) {
			List<BBP_Site> shakemapSites = BBP_Site.readFile(refShakeMapZip.getParentFile());
			File sitesXML = new File(refShakeMapZip.getParentFile(), "sites.xml");
			GriddedRegion reg = ShakemapPlotter.loadGriddedRegion(sitesXML);
			BBP_ShakeMapSimZipLoader loader = new BBP_ShakeMapSimZipLoader(refShakeMapZip, shakemapSites);
			gen.setShakeMap(loader, reg, shakemapSites, shakemapGMPE);
		}
		
		gen.generatePage(event, eventBBPDir, sites, rebuildGIF, rebuildMaps);
		
		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(outputDir);
	}

}
