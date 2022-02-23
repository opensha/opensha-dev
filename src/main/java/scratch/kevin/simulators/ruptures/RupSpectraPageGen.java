package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.apache.commons.collections4.map.HashedMap;
import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.DocumentException;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.commons.util.FileUtils;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator;
import org.opensha.sha.simulators.srf.SRF_PointData;
import org.opensha.sha.simulators.utils.RSQSimUtils;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_SimZipLoader.BBP_RupGenSimZipLoader;
import scratch.kevin.bbp.BBP_SimZipLoader.BBP_ShakeMapSimZipLoader;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;
import scratch.kevin.bbp.BBP_Wrapper;
import scratch.kevin.bbp.SeismogramPlotter;
import scratch.kevin.bbp.ShakemapPlotter;
import scratch.kevin.bbp.SpectraPlotter;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.plots.RuptureVelocityPlot;

class RupSpectraPageGen {
	
	private RSQSimCatalog catalog;
	private File outputDir;
	private double timeScale;
	private boolean scaleVelocities;
	
	private VelocityModel vm;
	
	private ScalarIMR[] spectraGMPEs;
	private File refBBPDir;
	private String refName;
	private EqkRupture gmpeRup;
	
	private BBP_ShakeMapSimZipLoader shakemap;
	private GriddedRegion shakemapRegion;
	private List<BBP_Site> shakemapSites;
	private ScalarIMR shakemapGMPE;
	
	private List<File> gpShakemapFiles;
	private GriddedRegion gpShakemapRegion;
	private List<BBP_Site> gpShakemapSites;

	public RupSpectraPageGen(RSQSimCatalog catalog, File outputDir, double timeScale, boolean scaleVelocities, VelocityModel vm) {
		this.catalog = catalog;
		this.outputDir = outputDir;
		this.timeScale = timeScale;
		this.scaleVelocities = scaleVelocities;
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
	
	public void setGPShakeMaps(List<File> shakemapFiles, GriddedRegion shakemapRegion,
			List<BBP_Site> shakemapSites) {
		this.gpShakemapFiles = shakemapFiles;
		this.gpShakemapRegion = shakemapRegion;
		this.gpShakemapSites = shakemapSites;
	}
	
	private static DecimalFormat twoDigitsDF = new DecimalFormat("0.00");
	
	public void generatePage(RSQSimEvent event, File eventBBPDir, List<BBP_Site> sites,
			boolean rebuildGIF, boolean rebuildMaps) throws IOException, GMT_MapException {
		int eventID = event.getID();
		String eventDirName = "event_"+eventID;
		if (timeScale != 1d) {
			eventDirName += "_timeScale"+(float)timeScale;
			if (scaleVelocities)
				eventDirName += "_velScale";
		}
		boolean refAdjustDDW = refBBPDir != null && refBBPDir.getName().contains("adjustDDW");
		if (refAdjustDDW)
			eventDirName += "_adjustDDW";
		File eventDir = new File(outputDir, eventDirName);
		Preconditions.checkState(eventDir.exists() || eventDir.mkdir());
		File resourcesDir = new File(eventDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		LinkedList<String> lines = new LinkedList<>();
		
		// header
		if (timeScale != 1d) {
			String title = "Event "+eventID+", M"+twoDigitsDF.format(event.getMagnitude())+", Time Scale Factor: "+(float)+timeScale;
			if (scaleVelocities)
				title += ", Velocities Scaled";
			if (refAdjustDDW)
				title += ", GP DDW Extension";
			lines.add("# "+title);
			lines.add("");
			String str = "**NOTE: RSQSim Slip/Time function has been MODIFIED on this page. "
					+ "Slip initiation for each patch are modified to be "+(float)timeScale+" times sooner.";
			if (scaleVelocities)
				str += " Slip velocities are also scaled by "+(float)+timeScale+"x, resulting in faster (but shorter lasting) slip.";
			else
				str += " Velocities remain the same, the slip time function is just shifted in time.";
			str += "**";
			lines.add(str);
			
			// also scale event first slip times
			for (EventRecord rec : event)
				rec.scaleElementTimeFirstSlips(timeScale, event.getTime());
		} else {
			String title = "Event "+eventID+", M"+twoDigitsDF.format(event.getMagnitude());
			if (refAdjustDDW)
				title += ", GP DDW Extension";
			lines.add("# "+title);
			lines.add("");
		}
		
		// add SRF download
		String srfPrefix = eventBBPDir.getName().replaceAll("_bbp", "");
		File srfFile = new File(eventBBPDir.getParentFile(), srfPrefix+".srf");
		Preconditions.checkState(srfFile.exists(), "SRF file doesn't exist: %s", srfFile.getAbsolutePath());
		File srfDest = new File(resourcesDir, srfFile.getName());
		Files.copy(srfFile, srfDest);
		
		BBP_SourceFile bbpSource = null;
		Location[] bbpSourceRect = null;
		Location bbpSourceHypo = null;
		if (refBBPDir != null) {
			File bbpSourceFile = locateSourceFile(refBBPDir);
			bbpSource = BBP_SourceFile.readFile(bbpSourceFile);
			File srcDest = new File(resourcesDir, bbpSourceFile.getName());
			Files.copy(bbpSourceFile, srcDest);
			lines.add("Download Rupture Files: [SRF]("+resourcesDir.getName()+"/"+srfFile.getName()
				+") [SRC]("+resourcesDir.getName()+"/"+srcDest.getName()+")");
			lines.add("");
			
			bbpSourceRect = bbpSource.getSurface().getRectangle();
			bbpSourceHypo = bbpSource.getHypoLoc();
			if (refAdjustDDW) {
				lines.add("**NOTE: "+refName+" DDW has been modified for the area to match the Somerville (2006) relationship:**");
				lines.add("");
				BBP_PlanarSurface bbpSurface = bbpSource.getSurface();
				double origLen = bbpSurface.getLength();
				double eventArea = event.getArea()*1e-6;
				double origDDW = eventArea/origLen;
				double newArea = origLen * bbpSurface.getWidth();
				lines.add("* Original surface: "+twoDigitsDF.format(origLen)+" km x "+twoDigitsDF.format(origDDW)
					+" km  = "+twoDigitsDF.format(eventArea)+" km^2");
				lines.add("* Modified surface: "+twoDigitsDF.format(origLen)+" km x "+twoDigitsDF.format(bbpSurface.getWidth())
					+" km = "+twoDigitsDF.format(newArea)+" km^2");
				lines.add("");
			}
		} else {
			lines.add("Download Rupture Files: [SRF]("+resourcesDir.getName()+"/"+srfFile.getName()+")");
			lines.add("");
		}
		lines.add("");
		lines.add("[Catalog Details](../../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
		lines.add("");
		lines.add("");
		lines.add("");
		BBP_RupGenSimZipLoader refLoader = null;
		if (refBBPDir != null) {
			File zipFile = new File(refBBPDir, "results.zip");
			refLoader = new BBP_RupGenSimZipLoader(zipFile, sites);
			lines.add("**"+refName+" Simulations: "+refLoader.getNumSims()+" per site**");
			System.out.println("Found "+refLoader.getNumSims()+" sims in "+zipFile.getAbsolutePath());
			lines.add("");
		}
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Fault List");
		lines.add(topLink); lines.add("");
		TableBuilder table = MarkdownUtils.tableBuilder();
		table.addLine("Section Name", "Area Ruptured", "Patches Ruptured", "Moment",
				"Equiv. Mag", "Max Slip", "Time First Slip");
		Map<String, Double> sectAreas = new HashedMap<>();
		Map<String, Integer> sectPatchCounts = new HashMap<>();
		Map<String, Double> sectMoments = new HashMap<>();
		Map<String, Double> sectMaxSlips = new HashMap<>();
		Map<String, Double> sectTimeFirsts = new HashMap<>();
		List<SimulatorElement> elems = event.getAllElements();
		double[] elemSlips = event.getAllElementSlips();
		double[] elemTimes = event.getAllElementTimes();
		double eventTime = event.getTime();
		List<? extends FaultSection> u3SubSects = catalog.getU3SubSects();
		int subSectOffset = RSQSimUtils.getSubSectIndexOffset(catalog.getElements(), u3SubSects);
		for (int i=0; i<elems.size(); i++) {
			SimulatorElement elem = elems.get(i);
			String fault = u3SubSects.get(elem.getSectionID()-subSectOffset).getParentSectionName();
			Double totArea, totMoment, maxSlip, timeFirst;
			Integer patchCount;
			if (sectAreas.containsKey(fault)) {
				totArea = sectAreas.get(fault);
				totMoment = sectMoments.get(fault);
				patchCount = sectPatchCounts.get(fault);
				maxSlip = sectMaxSlips.get(fault);
				timeFirst = sectTimeFirsts.get(fault);
			} else {
				totArea = 0d;
				totMoment = 0d;
				patchCount = 0;
				maxSlip = 0d;
				timeFirst = Double.POSITIVE_INFINITY;
			}
			sectTimeFirsts.put(fault, Math.min(timeFirst, elemTimes[i] - eventTime));
			sectAreas.put(fault, totArea + elem.getArea());
			sectMoments.put(fault, totMoment + FaultMomentCalc.getMoment(elem.getArea(), elemSlips[i]));
			sectPatchCounts.put(fault, patchCount + 1);
			sectMaxSlips.put(fault, Math.max(maxSlip, elemSlips[i]));
		}
		List<String> faultNames = new ArrayList<>();
		List<Double> faultMoments = new ArrayList<>();
		double totalMoment = 0d;
		double totalArea = 0d;
		int totalPatches = 0;
		double maxSlips = 0d;
		for (String fault : sectMoments.keySet()) {
			faultNames.add(fault);
			double moment = sectMoments.get(fault);
			faultMoments.add(moment);
			totalMoment += moment;
			totalArea += sectAreas.get(fault);
			totalPatches += sectPatchCounts.get(fault);
			maxSlips = Math.max(maxSlips, sectMaxSlips.get(fault));
		}
		faultNames = ComparablePairing.getSortedData(faultMoments, faultNames);
		Collections.reverse(faultNames);
		DecimalFormat momentDF = new DecimalFormat("0.00E0");
		if (faultNames.size() > 1) {
			// add total
			String totalName = "*(Total)*";
			faultNames.add(0, totalName);
			sectAreas.put(totalName, totalArea);
			sectPatchCounts.put(totalName, totalPatches);
			sectMoments.put(totalName, totalMoment);
			sectMaxSlips.put(totalName, maxSlips);
			sectTimeFirsts.put(totalName, 0d);
		}
		for (String fault : faultNames) {
			table.initNewLine().addColumn(fault);
			table.addColumn(twoDigitsDF.format(sectAreas.get(fault)*1e-6)+" km^2"); // to km^2
			table.addColumn(sectPatchCounts.get(fault).toString());
			double moment = sectMoments.get(fault);
			table.addColumn(momentDF.format(moment).replace("E", "e")+" N-m");
			double mag = MagUtils.momentToMag(moment);
			table.addColumn("M"+twoDigitsDF.format(mag));
			table.addColumn(twoDigitsDF.format(sectMaxSlips.get(fault))+" m");
			table.addColumn(twoDigitsDF.format(sectTimeFirsts.get(fault))+" s");
			table.finalizeLine();
		}
		lines.addAll(table.build());
		lines.add("");
		
		lines.add("## Rupture Plots");
		lines.add(topLink); lines.add("");
		lines.add("**Legend**");
		lines.add("* Colored, Filled Triangles: RSQSim Elements");
		lines.add("* Red Star: RSQSim Hypocenter");
		if (refBBPDir != null) {
			lines.add("* Dark Green Solid Outline: BBP Equivalent Planar Surface");
			lines.add("* Green Star: BBP Equivalent Hypocenter");
		}
		RuptureSurface gmpeSurf = null;
		if (gmpeRup != null) {
			gmpeSurf = gmpeRup.getRuptureSurface();
			lines.add("* Dark Gray Dashed Outline: GMPE Surface");
		}
		lines.add("");
		RSQSimEventSlipTimeFunc func = catalog.getSlipTimeFunc(event);
		if (timeScale != 1d)
			func = func.getTimeScaledFunc(timeScale, scaleVelocities);
		// slip time plot
		String rupPlotPrefix = "rupture_plot_"+eventID;
		RupturePlotGenerator.writeSlipPlot(event, func, resourcesDir, rupPlotPrefix, bbpSourceRect, bbpSourceHypo, gmpeSurf);
		// write grayscale friendly version
		RupturePlotGenerator.writeSlipPlot(event, func, resourcesDir, rupPlotPrefix+"_pub",
				null, null, null, true, false);
		File rupPlot = new File(resourcesDir, rupPlotPrefix+".png");
		Preconditions.checkState(rupPlot.exists());
		lines.add("### Slip/Time Plot");
		lines.add(topLink); lines.add("");
		lines.add("![Slip/Time Plot]("+resourcesDir.getName()+"/"+rupPlot.getName()+")");
		
		// animation
		if (!refAdjustDDW) {
			File rupAnim = new File(resourcesDir, rupPlotPrefix+".gif");
			if (rebuildGIF || !rupAnim.exists()) {
				if (bbpSourceRect == null)
					RupturePlotGenerator.writeSlipAnimation(event, func, rupAnim, 10);
				else
					RupturePlotGenerator.writeSlipAnimation(event, func, rupAnim, 10, bbpSourceRect[0], bbpSourceRect[1]);
			}
			Preconditions.checkState(rupAnim.exists());
			lines.add("### Slip/Vel Animation");
			lines.add(topLink); lines.add("");
			lines.add("[Click here to view Slip/Velocity Animation]("+resourcesDir.getName()+"/"+rupAnim.getName()+")");
		}
		
		// map view plot
		String mapRupPlotPrefix = "rupture_map_plot_"+eventID;
		double[] elementSlips = event.getAllElementSlips();
		CPT slipCPT = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, StatUtils.max(elementSlips));
		RupturePlotGenerator.writeMapPlot(catalog.getElements(), event, func, resourcesDir, mapRupPlotPrefix,
				bbpSourceRect, bbpSourceHypo, gmpeSurf, elementSlips, slipCPT, "Slip (m)");
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
		lines.add("### RSQSim Patch Slip-Time Histories");
		lines.add(topLink); lines.add("");
		
		lines.add("Slip-time histories on a subset of patches which participate in the event. The top panel "
				+ "shows cumulative slipi, and the bottom shows instantaneous velocity. Raw slip-time functions "
				+ "are in black, and interpolated/discretized for SRF files representation are colored.");
		lines.add("");
		
		int targetPatches = 20;
		int mod = Integer.max(1, elems.size()/targetPatches);
		
		table = MarkdownUtils.tableBuilder();
		
		table.initNewLine();
		for (int i=0; i<elems.size(); i++) {
			if (i % mod > 0)
				continue;
			SimulatorElement elem = elems.get(i);
			String prefix = "rsqsim_patch_"+elem.getID();
			RSQSimSRFGenerator.plotSlip(resourcesDir, prefix, event, func, elem, RSQSimBBP_Config.SRF_DT, false,
					RSQSimBBP_Config.SRF_INTERP_MODE);
			RSQSimSRFGenerator.plotSlip(resourcesDir, prefix+"_pub", event, func, elem, RSQSimBBP_Config.SRF_DT, true,
					RSQSimBBP_Config.SRF_INTERP_MODE);
			table.addColumn("![plot](resources/"+prefix+".png)");
		}
		table.finalizeLine();
		lines.addAll(table.wrap(4, 0).build());
		lines.add("");
		
		if (refBBPDir != null) {
			lines.add("");
			lines.add("### GP Comparison Patch Slip-Time Histories");
			lines.add(topLink); lines.add("");
			
			lines.add("Same as above, but for a Graves & Pitarka rupture.");
			lines.add("");
			
			File refSRF = new File(resourcesDir, "ref_srf.srf");
			if (!refSRF.exists()) {
				System.out.println("Generating ref SRF");
				File srcFile = null;
				for (File file : refBBPDir.listFiles()) {
					if (file.getName().toLowerCase().endsWith(".src")) {
						srcFile = file;
						break;
					}
				}
				Preconditions.checkNotNull("SRC file not found in "+refBBPDir.getAbsolutePath());
				File tempDir = FileUtils.createTempDir();
				BBP_Wrapper wrapper = new BBP_Wrapper(RSQSimBBP_Config.VM, RSQSimBBP_Config.METHOD,
						srcFile, null, null, null, tempDir);
				wrapper.setSRFGenOnly(true);
				wrapper.run();
				File tmpSRF = null;
				for (File file : tempDir.listFiles()) {
					if (file.getName().toLowerCase().endsWith(".srf")) {
						tmpSRF = file;
						break;
					}
				}
				Preconditions.checkNotNull("SRF file not found in "+tempDir.getAbsolutePath());
				Files.copy(tmpSRF, refSRF);
			}
			
			List<SRF_PointData> points = SRF_PointData.readSRF(refSRF);
			double maxTime = 0d;
			for (SRF_PointData point : points)
				maxTime = Math.max(maxTime, point.getEndTime());
			
			mod = Integer.max(1, points.size()/targetPatches);
			
			table = MarkdownUtils.tableBuilder();
			
			table.initNewLine();
			for (int i=0; i<points.size(); i++) {
				if (i % mod > 0)
					continue;
				SRF_PointData point = points.get(i);
				String prefix = "gp_patch_"+i;
				RSQSimSRFGenerator.plotSlip(resourcesDir, prefix, point, maxTime, false);
				table.addColumn("![plot](resources/"+prefix+".png)");
				RSQSimSRFGenerator.plotSlip(resourcesDir, prefix+"_zoom", point, maxTime, true);
			}
			table.finalizeLine();
			lines.addAll(table.wrap(4, 0).build());
			lines.add("");
		}
		
		lines.add("");
		lines.add("## Spectra Plots");
		lines.add(topLink); lines.add("");
		
		File noScaleBBPDir = null;
		if (timeScale != 1d) {
			noScaleBBPDir = RSQSimBBP_Config.getEventBBPDir(catalog, eventID, RSQSimBBP_Config.SRF_INTERP_MODE,
					RSQSimBBP_Config.SRF_DT, 1d, false);
			if (!noScaleBBPDir.exists())
				noScaleBBPDir = null;
			else
				System.out.println("Found unscaled reference BBP dir: "+noScaleBBPDir.getAbsolutePath());
		}
		
		for (BBP_Site site : sites) {
			String siteName = RSQSimBBP_Config.siteCleanName(site);
			System.out.println("Site: "+siteName);
			lines.add("### Site "+siteName);
			lines.add(topLink); lines.add("");
			Location loc = site.getLoc();
			lines.add("*Location: "+(float)loc.getLatitude()+", "+(float)loc.getLongitude()+"*");
			
			// calculate distances
			table = MarkdownUtils.tableBuilder();
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
					
					DiscretizedFunc[] otherSpectra = new DiscretizedFunc[0];
					if (noScaleBBPDir != null) {
						try {
							File noScaleFASFile = SpectraPlotter.findFASFile(noScaleBBPDir, site.getName());
							otherSpectra = new DiscretizedFunc[] { SpectraPlotter.loadFAS(noScaleFASFile) };
							otherSpectra[0].setName("RSQSim (unscaled)");
						} catch (Exception e) {
						}
					}
					
					SpectraPlotter.plotMultiFAS(refSpectra, refName, dataSpectra, "RSQSim", title, resourcesDir, prefix, otherSpectra);
				} else {
					SpectraPlotter.plotFAS(fasFile, resourcesDir, prefix);
				}
				lines.add("!["+siteName+" FAS Plot]("+resourcesDir.getName()+"/"+prefix+".png)");
				
			}
			
			boolean hasRD100;
			File rotDFile;
			try {
				rotDFile = SpectraPlotter.findRotD100File(eventBBPDir, site.getName());
				hasRD100 = true;
			} catch (FileNotFoundException e) {
				rotDFile = SpectraPlotter.findRotD50File(eventBBPDir, site.getName());
				hasRD100 = false;
			}
			UncertainArbDiscFunc[] gmpeSpectra = null;
			if (spectraGMPEs != null && spectraGMPEs.length > 0 && gmpeRup != null) {
				System.out.print("Calculating GMPEs...");
				gmpeSpectra = new UncertainArbDiscFunc[spectraGMPEs.length];
				for (int i=0; i<gmpeSpectra.length; i++)
					gmpeSpectra[i] = SpectraPlotter.calcGMPE_RotD50(gmpeRup, site, spectraGMPEs[i], vm);
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
				DiscretizedFunc dataSpectra = SpectraPlotter.loadRotD50(rotDFile);
				
				DiscretizedFunc[] otherSpectra = new DiscretizedFunc[0];
				if (noScaleBBPDir != null) {
					try {
						File noScaleRotDFile;
						if (hasRD100)
							noScaleRotDFile = SpectraPlotter.findRotD100File(noScaleBBPDir, site.getName());
						else
							noScaleRotDFile = SpectraPlotter.findRotD50File(noScaleBBPDir, site.getName());
						otherSpectra = new DiscretizedFunc[] { SpectraPlotter.loadRotD50(noScaleRotDFile) };
						otherSpectra[0].setName("RSQSim (unscaled)");
					} catch (Exception e) {
					}
				}
				
				SpectraPlotter.plotMultiRotD50(refSpectra, refName, dataSpectra, "RSQSim", title,
						resourcesDir, prefix, gmpeSpectra, otherSpectra);
			} else {
				SpectraPlotter.plotRotD(rotDFile, resourcesDir, prefix, true, hasRD100, gmpeSpectra);
			}
			lines.add("!["+siteName+" RotD50 Plot]("+resourcesDir.getName()+"/"+prefix+".png)");
			
			if (hasRD100) {
				lines.add("#### "+siteName+" RotD Ratio");
				lines.add(topLink); lines.add("");
				title = "Event "+eventID+" "+siteName+" RotD100/RotD50 Ratio";
				prefix = "rotd_ratio_"+site.getName();
				if (refLoader != null) {
					List<DiscretizedFunc[]> refSpectra = new ArrayList<>();
					int numSims = refLoader.getNumSims();
					for (int i=0; i<numSims; i++)
						refSpectra.add(refLoader.readRotD(site, i));
					DiscretizedFunc[] dataSpectra = SpectraPlotter.loadRotD(rotDFile);
					SpectraPlotter.plotRotDRatio(dataSpectra, "RSQSim", refSpectra, refName, title, resourcesDir, prefix);
				} else {
					SpectraPlotter.plotRotDRatio(rotDFile, "RSQSim", title, resourcesDir, prefix);
				}
				lines.add("!["+siteName+" RotD Ratio Plot]("+resourcesDir.getName()+"/"+prefix+".png)");
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
		
		if (shakemap != null && !refAdjustDDW) {
			lines.add("");
			lines.add("## ShakeMaps");
			lines.add(topLink); lines.add("");
			
			double[] periods;
			if (shakemap.hasPGV())
				periods = new double[] { -1, 1d, 2d, 3d, 4d, 5d, 7.5, 10d };
			else
				periods = new double[] { 1d, 2d, 3d, 4d, 5d, 7.5, 10d };
			String prefix = "shakemap";
			File[] mapFiles = new File[periods.length];
			File[] gmpeMapFiles = null;
			File[] ratioMapFiles = null;
			if (shakemapGMPE != null) {
				gmpeMapFiles = new File[periods.length];
				ratioMapFiles = new File[periods.length];
			}
			File[] rd100MapFiles = null;
			File[] rdRatioMapFiles = null;
			if (shakemap.hasRotD100()) {
				rd100MapFiles = new File[periods.length];
				rdRatioMapFiles = new File[periods.length];
			}
			for (int i=0; i<periods.length; i++) {
				double p = periods[i];
				String name;
				if (p == -1d)
					name = "pgv";
				else if (p == Math.round(p))
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
				if (shakemap.hasRotD100()) {
					rd100MapFiles[i] = new File(resourcesDir, prefix+"_"+name+"_rd100.png");
					rdRatioMapFiles[i] = new File(resourcesDir, prefix+"_"+name+"_rd100_ratio.png");
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
					if (shakemap.hasRotD100())
						skip = skip && rd100MapFiles[i].exists() && rdRatioMapFiles[i].exists();
					if (!skip)
						retainedPeriods.add(periods[i]);
				}
				plotPeriods = Doubles.toArray(retainedPeriods);
			}
			
			if (periods.length > 0) {
				String label = "Event "+eventID;
				ShakemapPlotter.plotShakemaps(shakemap, shakemapRegion, shakemapSites, vm, label, resourcesDir,
						prefix, true, shakemapGMPE, gmpeRup, plotPeriods);
			} else {
				System.out.println("Skipping all maps!");
			}
			
			table = MarkdownUtils.tableBuilder();
			table.initNewLine().addColumn("").addColumn("RSQSim");
			if (shakemapGMPE != null)
				table.addColumn(shakemapGMPE.getShortName()).addColumn("Ratio");
			table.finalizeLine();
			for (int i=0; i<periods.length; i++) {
				if (periods[i] == -1d)
					table.initNewLine().addColumn("**PGV**");
				else
					table.initNewLine().addColumn("**"+(float)periods[i]+" s SA**");
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
			
			if (shakemap.hasRotD100()) {
				lines.add("### Region RotD100/RotD50 Ratio");
				lines.add(topLink); lines.add("");
				List<DiscretizedFunc[]> rds = new ArrayList<>();
				for (int i=0; i<shakemap.getNumSims(); i++)
					rds.add(shakemap.readRotD(i));
				String title = "Event "+eventID+" ShakeMap RotD100/RotD50 Ratio";
				prefix = "rotd_ratio_shakemap";
				SpectraPlotter.plotRotDRatio(rds, "RSQSim", title, resourcesDir, prefix);
				lines.add("![ShakeMap RotD Ratio Plot]("+resourcesDir.getName()+"/"+prefix+".png)");
				List<Double> dists = new ArrayList<>();
				for (BBP_Site shakeSite : shakemapSites) {
					Location loc = shakeSite.getLoc();
					if (CatalogGMPE_Compare.DIST_JB)
						dists.add(gmpeSurf.getDistanceJB(loc));
					else
						dists.add(gmpeSurf.getDistanceRup(loc));
				}
				String distLabel = CatalogGMPE_Compare.DIST_JB ? "rJB (km)" : "rRup (km)";
				prefix = "rotd_ratio_shakemap_dist_dependence";
				SpectraPlotter.plotRotDRatioDependence(rds, dists, distLabel, 10, periods,
						"RSQSim", title, resourcesDir, prefix, false);
				lines.add("");
				lines.add("![ShakeMap RotD Dist Dependence]("+resourcesDir.getName()+"/"+prefix+".png)");
				
				table = MarkdownUtils.tableBuilder();
				table.addLine("", "RotD50", "RotD100", "RotD100/RotD50 Ratio");
				for (int i=0; i<periods.length; i++) {
					if (periods[i] == -1d)
						table.initNewLine().addColumn("**PGV**");
					else
						table.initNewLine().addColumn("**"+(float)periods[i]+" s SA**");
					table.addColumn("![RotD50 Map]("+resourcesDir.getName()+"/"+mapFiles[i].getName()+")");
					table.addColumn("![RotD100 Map]("+resourcesDir.getName()+"/"+rd100MapFiles[i].getName()+")");
					table.addColumn("![RotD Ratio Map]("+resourcesDir.getName()+"/"+rdRatioMapFiles[i].getName()+")");
					table.finalizeLine();
				}
				
				lines.add("");
				lines.addAll(table.build());
			}
		}
		
		if (gpShakemapFiles != null && !gpShakemapFiles.isEmpty()) {
			lines.add("");
			lines.add("## G&P ShakeMaps");
			lines.add(topLink); lines.add("");
			
			double[] periods = { 1d, 2d, 3d, 4d, 5d, 7.5, 10d };
			
			List<File[]> mapFilesList = new ArrayList<>();
			List<File[]> rd100MapFilesList = new ArrayList<>();
			List<File[]> rdRatioMapFilesList = new ArrayList<>();
			
			List<String> seeds = new ArrayList<>();
			
			for (File gpShakemapFile : gpShakemapFiles) {
				String parentName = gpShakemapFile.getParentFile().getName();
				System.out.println("Processing GP ShakeMaps for "+parentName);
				Preconditions.checkState(parentName.contains("-seed"));
				String seed = parentName.substring(parentName.indexOf("-seed")+1);
				if (seed.contains("-"))
					seed = seed.substring(0, seed.indexOf("-"));
				String prefix = "gp_shakemap_"+seed;
				
				seed = seed.replaceAll("seed", "");
				System.out.print("Seed: "+seed);
				seeds.add(seed);
				
				BBP_ShakeMapSimZipLoader shakemap =
						new BBP_ShakeMapSimZipLoader(gpShakemapFile, gpShakemapSites);
				
				File[] mapFiles = new File[periods.length];
				mapFilesList.add(mapFiles);
				File[] rd100MapFiles = null;
				File[] rdRatioMapFiles = null;
				if (shakemap.hasRotD100()) {
					rd100MapFiles = new File[periods.length];
					rdRatioMapFiles = new File[periods.length];
					rd100MapFilesList.add(rd100MapFiles);
					rdRatioMapFilesList.add(rdRatioMapFiles);
				}
				for (int i=0; i<periods.length; i++) {
					double p = periods[i];
					String name;
					if (p == Math.round(p))
						name = (int)p+"s";
					else
						name = (float)p+"s";
					mapFiles[i] = new File(resourcesDir, prefix+"_"+name+".png");
					if (shakemap.hasRotD100()) {
						rd100MapFiles[i] = new File(resourcesDir, prefix+"_"+name+"_rd100.png");
						rdRatioMapFiles[i] = new File(resourcesDir, prefix+"_"+name+"_rd100_ratio.png");
					}
				}
				
				double[] plotPeriods = periods;
				if (!rebuildMaps) {
					// see if we can skip any...
					List<Double> retainedPeriods = new ArrayList<>();
					for (int i=0; i<periods.length; i++) {
						boolean skip = mapFiles[i].exists();
						if (shakemap.hasRotD100())
							skip = skip && rd100MapFiles[i].exists() && rdRatioMapFiles[i].exists();
						if (!skip)
							retainedPeriods.add(periods[i]);
					}
					plotPeriods = Doubles.toArray(retainedPeriods);
				}
				
				if (periods.length > 0) {
					String label = "GP Seed "+seed;
					ShakemapPlotter.plotShakemaps(shakemap, gpShakemapRegion, gpShakemapSites, vm, label, resourcesDir,
							prefix, true, null, null, plotPeriods);
				} else {
					System.out.println("Skipping all maps!");
				}
			}
			
			table = MarkdownUtils.tableBuilder();
			table.initNewLine().addColumn("SA Period");
			for (String seed : seeds)
				table.addColumn("G&P Seed "+seed);
			table.finalizeLine();
			for (int i=0; i<periods.length; i++) {
				table.initNewLine().addColumn("**"+(float)periods[i]+" s**");
				for (File[] mapFiles : mapFilesList)
					table.addColumn("![GP Map]("+resourcesDir.getName()+"/"+mapFiles[i].getName()+")");
				table.finalizeLine();
			}
			
			lines.add("");
			lines.addAll(table.build());
			
			for (int s=0; s<seeds.size(); s++) {
				File[] mapFiles = mapFilesList.get(s);
				File[] rd100MapFiles = rd100MapFilesList.get(s);
				File[] rdRatioMapFiles = rdRatioMapFilesList.get(s);
				if (rd100MapFiles != null) {
					lines.add("### GP Region RotD100/RotD50 Ratio, Seed "+seeds.get(s));
					lines.add(topLink); lines.add("");
					table = MarkdownUtils.tableBuilder();
					table.addLine("SA Period", "RotD50", "RotD100", "RotD100/RotD50 Ratio");
					for (int i=0; i<periods.length; i++) {
						table.initNewLine().addColumn("**"+(float)periods[i]+" s**");
						table.addColumn("![RotD50 Map]("+resourcesDir.getName()+"/"+mapFiles[i].getName()+")");
						table.addColumn("![RotD100 Map]("+resourcesDir.getName()+"/"+rd100MapFiles[i].getName()+")");
						table.addColumn("![RotD Ratio Map]("+resourcesDir.getName()+"/"+rdRatioMapFiles[i].getName()+")");
						table.finalizeLine();
					}
					
					lines.add("");
					lines.addAll(table.build());
				}
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, eventDir);
	}
	
	private File locateSourceFile(File bbpDir) throws FileNotFoundException {
		for (File file : bbpDir.listFiles())
			if (file.getName().toLowerCase().endsWith(".src"))
				return file;
		throw new FileNotFoundException("No .src files found in "+bbpDir.getAbsolutePath());
	}
	
	public static void main(String[] args) throws IOException, DocumentException, GMT_MapException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		File outputDir = new File("/home/kevin/markdown/rsqsim-analysis/catalogs");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");
		String refName = "Graves & Pitarka (2016)";
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
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2320.instance(baseDir);
//		int eventID = 6195527;
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2336.instance(baseDir);
//		int eventID = 131670;
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2337.instance(baseDir);
//		int eventID = 203769;
		
//		RSQSimCatalog catalog = Catalogs.JG_2194_K2.instance(baseDir);
//		int eventID = 18840012;
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2194_LONG.instance(baseDir);
//		int eventID = 526885;
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2585.instance(baseDir);
////		int eventID = 81854;
//		int eventID = 1670183;
////		int eventID = 2637969;
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2740.instance(baseDir);
//		int eventID = 385955;
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2829.instance(baseDir);
////		int eventID = 5304;
//		int eventID = 31324;
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
//		int eventID = 9955310;
////		int eventID = 3817386;

//		RSQSimCatalog catalog = Catalogs.BRUCE_4322.instance(baseDir);
//		int eventID = 40636;
////		int eventID = 92236;

//		RSQSimCatalog catalog = Catalogs.BRUCE_4655.instance(baseDir);
////		int eventID = 2106470;
////		int eventID = 1029675;
////		int eventID = 467544;
////		int eventID = 923866;
//		int eventID = 1318657;

//		RSQSimCatalog catalog = Catalogs.BRUCE_4827.instance(baseDir);
////		int eventID = 195167;
//		int eventID = 128149;

//		RSQSimCatalog catalog = Catalogs.BRUCE_4841.instance(baseDir);
//		int eventID = 755070;
////		int eventID = 2441060;

//		RSQSimCatalog catalog = Catalogs.BRUCE_4853.instance(baseDir);
//		int eventID = 516171;

//		RSQSimCatalog catalog = Catalogs.BRUCE_4849.instance(baseDir);
//		int eventID = 80686;
		
//		RSQSimCatalog catalog = Catalogs.TEST_DOUBLE_4860.instance(baseDir);
////		int eventID = 4853;
//		int eventID = 15377;
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_4860_10X.instance(baseDir);
////		int eventID = 51863;
////		int eventID = 39055;
//		int eventID = 12581;
////		int eventID = 77272;
////		int eventID = 41890; // fig 4

//		RSQSimCatalog catalog = Catalogs.BRUCE_4950.instance(baseDir);
////		int eventID = 469461;
//		int eventID = 368122; // event with lowest ave slip ratio from trans to dlist

//		RSQSimCatalog catalog = Catalogs.BRUCE_4983.instance(baseDir);
//		int eventID = 1499589;

		RSQSimCatalog catalog = Catalogs.BRUCE_4983_STITCHED.instance(baseDir);
//		int eventID = 3092204; // this is not an SAF event, multi fault
//		int eventID = 7119753;
//		int eventID = 8242900;
//		int eventID = 1499589;
//		int eventID = 7028377;
//		int eventID = 13383629;
//		int eventID = 6553169;
//		int eventID = 13272163;
//		int eventID = 1651575;
		int eventID = 3015199;
//		int eventID = 1879413;
//		int eventID = 10096082;
//		int eventID = 7992279;
//		int eventID = 7748150;
//		int eventID = 3012841;
//		int eventID = 2614773;
//		int eventID = 1809975;
		
		double timeScale = 1d;
		boolean scaleVelocities = true;
		boolean gpAdjustDDW = false;
		
		File eventBBPDir = RSQSimBBP_Config.getEventBBPDir(catalog, eventID, RSQSimBBP_Config.SRF_INTERP_MODE,
				RSQSimBBP_Config.SRF_DT, timeScale, scaleVelocities);
		System.out.println("Event BBP dir: "+eventBBPDir.getAbsolutePath()+" exists ? "+eventBBPDir.exists());
		
		// find BBP parallel dir
		File refBBPDir = null;
		File[] refDirs = bbpParallelDir.listFiles();
		Arrays.sort(refDirs, new FileNameComparator());
		for (File dir : refDirs) {
			String name = dir.getName();
			if (dir.isDirectory() && name.contains(catalog.getCatalogDir().getName()) && name.contains(eventID+"")
					&& name.contains("-gp-") && !name.contains("shakemap")) {
				if ((gpAdjustDDW && !name.contains("adjustDDW")) || (!gpAdjustDDW && name.contains("adjustDDW")))
					continue;
				File zipFile = new File(dir, "results.zip");
				if (zipFile.exists())
					refBBPDir = dir;
			}
		}
		VelocityModel vm;
		if (refBBPDir != null) {
			System.out.println("Located ref BBP dir: "+refBBPDir.getAbsolutePath());
			vm = RSQSimBBP_Config.detectVM(refBBPDir);
		} else {
			System.out.println("No BBP comparison dir found");
			vm = RSQSimBBP_Config.VM;
		}
		
		File refShakeMapZip = null;
		for (File dir : refDirs) {
			String name = dir.getName();
			if (dir.isDirectory() && name.contains(catalog.getCatalogDir().getName()) && name.contains(eventID+"")
					&& name.contains("-shakemap-")) {
				if (timeScale != 1d && !name.contains("-timeScale"+(float)timeScale) && (!scaleVelocities || !name.contains("-velScale")))
					continue;
				if (name.contains("-gp-"))
					continue; // TODO
				File shakemapFile = new File(dir, "results_rotD.zip");
				if (!shakemapFile.exists())
					shakemapFile = new File(dir, "results.zip");
				if (shakemapFile.exists())
					refShakeMapZip = shakemapFile;
			}
		}
		if (refShakeMapZip != null)
			System.out.println("Located ref ShakeMap zip: "+refShakeMapZip.getAbsolutePath());
		
		List<File> gpShakeMapFiles = new ArrayList<>();
		for (File dir : refDirs) {
			String name = dir.getName();
			if (dir.isDirectory() && name.contains(catalog.getCatalogDir().getName()) && name.contains(eventID+"")
					&& name.contains("-shakemap-") && name.contains("-gp-")) {
				if ((gpAdjustDDW && !name.contains("adjustDDW")) || (!gpAdjustDDW && name.contains("adjustDDW")))
					continue;
				File shakemapFile = new File(dir, "results_rotD.zip");
				if (!shakemapFile.exists())
					shakemapFile = new File(dir, "results.zip");
				if (shakemapFile.exists())
					gpShakeMapFiles.add(shakemapFile);
			}
		}
		if (!gpShakeMapFiles.isEmpty())
			System.out.println("Found "+gpShakeMapFiles.size()+" GP ShakeMaps");
		
		RSQSimEvent event = catalog.loader().byID(eventID);
		
		List<BBP_Site> sites = null;
		if (refBBPDir != null)
			sites = BBP_Site.readFile(refBBPDir);
		else
//			sites = RSQSimBBP_Config.getCyberShakeInitialLASites();
			sites = RSQSimBBP_Config.getCyberShakeVs500LASites();
		
		boolean runBBP = true;
		if (eventBBPDir.exists()) {
			// look for RotD100 files
			for (File file : eventBBPDir.listFiles()) {
				if (file.getName().endsWith(".rd100")) {
					runBBP = false;
					break;
				}
			}
		}
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
			RSQSimBBP_Config.runBBP(catalog, event, sites, timeScale, scaleVelocities);
		}
		
		ScalarIMR[] gmpes = { AttenRelRef.ASK_2014.instance(null), AttenRelRef.BSSA_2014.instance(null),
				AttenRelRef.CB_2014.instance(null), AttenRelRef.CY_2014.instance(null) };
		for (ScalarIMR gmpe : gmpes)
			gmpe.setParamDefaults();
		ScalarIMR shakemapGMPE = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);
		shakemapGMPE.setParamDefaults();
		
		File catalogOutputDir = new File(outputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		File vmDir = new File(catalogOutputDir, "bbp_"+vm.name());
		Preconditions.checkState(vmDir.exists() || vmDir.mkdir());
		
		RupSpectraPageGen gen = new RupSpectraPageGen(catalog, vmDir, timeScale, scaleVelocities, vm);
		
		EqkRupture gmpeRup = null;
		
//		if ((refBBPDir != null && gmpes != null && gmpes.length > 0)
//				|| (refShakeMapZip != null && shakemapGMPE != null))
			gmpeRup = catalog.getMappedSubSectRupture(event);
		
		if (refBBPDir != null)
			gen.setRefBBP(refBBPDir, refName, gmpes, gmpeRup);
		else
			gen.setRefBBP(null, null, gmpes, gmpeRup);
		
		if (refShakeMapZip != null) {
			List<BBP_Site> shakemapSites = BBP_Site.readFile(refShakeMapZip.getParentFile());
			File sitesXML = new File(refShakeMapZip.getParentFile(), "sites.xml");
			GriddedRegion reg = ShakemapPlotter.loadGriddedRegion(sitesXML);
			BBP_ShakeMapSimZipLoader loader = new BBP_ShakeMapSimZipLoader(refShakeMapZip, shakemapSites);
			gen.setShakeMap(loader, reg, shakemapSites, shakemapGMPE);
		}
		
		if (!gpShakeMapFiles.isEmpty()) {
			File dir0 = gpShakeMapFiles.get(0).getParentFile();
			List<BBP_Site> shakemapSites = BBP_Site.readFile(dir0);
			File sitesXML = new File(dir0, "sites.xml");
			GriddedRegion reg = ShakemapPlotter.loadGriddedRegion(sitesXML);
			gen.setGPShakeMaps(gpShakeMapFiles, reg, shakemapSites);
		}
		
		gen.generatePage(event, eventBBPDir, sites, rebuildGIF, rebuildMaps);
		
		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(outputDir);
	}

}
