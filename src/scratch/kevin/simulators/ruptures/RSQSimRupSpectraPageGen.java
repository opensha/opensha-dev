package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.opensha.commons.data.function.UncertainArbDiscDataset;
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

import scratch.kevin.MarkdownUtils;
import scratch.kevin.MarkdownUtils.TableBuilder;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.SpectraPlotter;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class RSQSimRupSpectraPageGen {
	
	private RSQSimCatalog catalog;
	private File outputDir;
	private ScalarIMR[] gmpes;
	private VelocityModel vm;

	public RSQSimRupSpectraPageGen(RSQSimCatalog catalog, File outputDir, ScalarIMR[] gmpes, VelocityModel vm) {
		this.catalog = catalog;
		this.outputDir = outputDir;
		this.gmpes = gmpes;
		this.vm = vm;
	}
	
	private static DecimalFormat twoDigitsDF = new DecimalFormat("0.00");
	
	public void generatePage(RSQSimEvent event, File eventBBPDir, File refBBPDir, String refName,
			EqkRupture gmpeRup, List<BBP_Site> sites) throws IOException {
		int eventID = event.getID();
		File eventDir = new File(outputDir, "event_"+eventID);
		Preconditions.checkState(eventDir.exists() || eventDir.mkdir());
		File resourcesDir = new File(eventDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		LinkedList<String> lines = new LinkedList<>();
		
		// header
		lines.add("# Event "+eventID+", M"+twoDigitsDF.format(event.getMagnitude()));
		lines.add("");
		lines.addAll(catalog.getMarkdownMetadataTable());
		lines.add("");
		
		int tocIndex = lines.size();
		
		lines.add("## Rupture Plots");
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
			bbpSourceRect = bbpSource.buildRectangle();
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
		lines.add("![Slip/Time Plot]("+resourcesDir.getName()+"/"+rupPlot.getName()+")");
		// map view plot
		String mapRupPlotPrefix = "rupture_map_plot_"+eventID;
		RupturePlotGenerator.writeMapPlot(catalog.getElements(), event, func, resourcesDir, mapRupPlotPrefix,
				bbpSourceRect, bbpSourceHypo, gmpeSurf);
		File rupMapPlot = new File(resourcesDir, mapRupPlotPrefix+".png");
		Preconditions.checkState(rupMapPlot.exists());
		lines.add("### Map Plot");
		lines.add("![Map Plot]("+resourcesDir.getName()+"/"+rupMapPlot.getName()+")");
		
		lines.add("");
		lines.add("## Spectra Plots");
		
		for (BBP_Site site : sites) {
			System.out.println("Site: "+site.getName());
			lines.add("## Site "+site.getName());
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
				horzLine.add(twoDigitsDF.format(bbpSource.getHorizontalDistance(loc))+" km");
			if (gmpeSurf != null)
				horzLine.add(twoDigitsDF.format(gmpeSurf.getDistanceJB(loc))+" km");
			table.addLine(horzLine);
			List<String> linearLine = new ArrayList<>();
			linearLine.add("3-D");
			linearLine.add(twoDigitsDF.format(linearDist)+" km");
			if (bbpSource != null)
				linearLine.add(twoDigitsDF.format(bbpSource.getLinearDistance(loc))+" km");
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
				lines.add("### "+site.getName()+" Fourrier Amplitude Spectra");
				String title = "Event "+eventID+" "+site.getName()+" FAS Spectra";
				String prefix = "fas_spectra_"+site.getName();
				if (refBBPDir.exists()) {
					List<File> refFiles = SpectraPlotter.findRunFiles(refBBPDir, site.getName(), false);
					SpectraPlotter.plotMultiFAS(refFiles, refName, fasFile, "RSQSim", title, resourcesDir, prefix);
				} else {
					SpectraPlotter.plotFAS(fasFile, resourcesDir, prefix);
				}
				lines.add("!["+site.getName()+" FAS Plot]("+resourcesDir.getName()+"/"+prefix+".png)");
				
			}
			
			File rotD50File = SpectraPlotter.findRotD50File(eventBBPDir, site.getName());
			if (rotD50File.exists()) {
				UncertainArbDiscDataset[] gmpeSpectra = null;
				if (gmpes != null && gmpes.length > 0) {
					System.out.print("Calculating GMPEs...");
					gmpeSpectra = new UncertainArbDiscDataset[gmpes.length];
					for (int i=0; i<gmpeSpectra.length; i++)
						gmpeSpectra[i] = SpectraPlotter.calcGMPE_RotD50(gmpeRup, loc, vm, gmpes[i]);
					System.out.println("DONE.");
				}
				
				lines.add("### "+site.getName()+" RotD50 Spectra");
				String title = "Event "+eventID+" "+site.getName()+" RotD50 Spectra";
				String prefix = "rotd50_spectra_"+site.getName();
				if (refBBPDir.exists()) {
					List<File> refFiles = SpectraPlotter.findRunFiles(refBBPDir, site.getName(), true);
					SpectraPlotter.plotMultiRotD50(refFiles, refName, rotD50File, "RSQSim", title,
							resourcesDir, prefix, gmpeSpectra);
				} else {
					SpectraPlotter.plotFAS(fasFile, resourcesDir, prefix);
				}
				lines.add("!["+site.getName()+" RotD50 Plot]("+resourcesDir.getName()+"/"+prefix+".png)");
			}
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
	
	private static List<BBP_Site> loadSites(File bbpDir) throws IOException {
		for (File file : bbpDir.listFiles())
			if (file.getName().toLowerCase().endsWith(".stl"))
				return BBP_Site.readFile(file);
		throw new FileNotFoundException("No .stl files found in "+bbpDir.getAbsolutePath());
	}
	
	public static void main(String[] args) throws IOException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		File outputDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");
		String refName = "Graves & Pitarka (2015)";
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2194_LONG.instance(baseDir);
//		int eventID = 136704;
		
		RSQSimCatalog catalog = Catalogs.JG_UCERF3_millionElement.instance(baseDir);
		int eventID = 4099020;
		
		File eventBBPDir = new File(new File(catalog.getCatalogDir(), "event_srfs"),
				"event_"+eventID+"_0.05s_ADJ_VEL_bbp");
		Preconditions.checkState(eventBBPDir.exists());
		
		// find BBP parallel dir
		File refBBPDir = null;
		File[] refDirs = bbpParallelDir.listFiles();
		Arrays.sort(refDirs, new FileNameComparator());
		for (File dir : refDirs) {
			String name = dir.getName();
			if (dir.isDirectory() && name.contains(catalog.getCatalogDir().getName()) && name.contains(eventID+""))
				refBBPDir = dir;
		}
		if (refBBPDir != null)
			System.out.println("Located ref BBP dir: "+refBBPDir.getAbsolutePath());
		
		List<BBP_Site> sites;
		if (refBBPDir != null)
			sites = loadSites(refBBPDir);
		else
			sites = loadSites(new File("/data/kevin/bbp/bbp_data/run"));
		
		ScalarIMR[] gmpes = { AttenRelRef.ASK_2014.instance(null), AttenRelRef.BSSA_2014.instance(null),
				AttenRelRef.CB_2014.instance(null), AttenRelRef.CY_2014.instance(null) };
		for (ScalarIMR gmpe : gmpes)
			gmpe.setParamDefaults();
		VelocityModel vm = VelocityModel.LA_BASIN;
		
		File catalogOutputDir = new File(outputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		catalog.writeMarkdownSummary(catalogOutputDir);
		
		RSQSimRupSpectraPageGen gen = new RSQSimRupSpectraPageGen(catalog, catalogOutputDir, gmpes, vm);
		
		RSQSimEvent event = catalog.loadEventByID(eventID);
		EqkRupture gmpeRup = null;
		if (gmpes != null && gmpes.length > 0)
			gmpeRup = catalog.getGMPE_Rupture(event, 0.2);
		
		gen.generatePage(event, eventBBPDir, refBBPDir, refName, gmpeRup, sites);
	}

}
