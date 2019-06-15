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
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.utils.RSQSimUtils;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_SimZipLoader;
import scratch.kevin.bbp.BBP_SimZipLoader.BBP_RupGenSimZipLoader;
import scratch.kevin.bbp.BBP_SimZipLoader.BBP_ShakeMapSimZipLoader;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;
import scratch.kevin.bbp.SeismogramPlotter;
import scratch.kevin.bbp.ShakemapPlotter;
import scratch.kevin.bbp.SpectraPlotter;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.plots.RuptureVelocityPlot;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;

class RupSpectraPageGen {
	
	private RSQSimCatalog catalog;
	private File outputDir;
	private double timeScale;
	private boolean scaleVelocities;
	
	private ScalarIMR[] spectraGMPEs;
	private File refBBPDir;
	private String refName;
	private EqkRupture gmpeRup;
	
	private BBP_ShakeMapSimZipLoader shakemap;
	private GriddedRegion shakemapRegion;
	private List<BBP_Site> shakemapSites;
	private ScalarIMR shakemapGMPE;

	public RupSpectraPageGen(RSQSimCatalog catalog, File outputDir, double timeScale, boolean scaleVelocities) {
		this.catalog = catalog;
		this.outputDir = outputDir;
		this.timeScale = timeScale;
		this.scaleVelocities = scaleVelocities;
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
		BBP_SourceFile bbpSource = null;
		Location[] bbpSourceRect = null;
		Location bbpSourceHypo = null;
		if (refBBPDir != null) {
			bbpSource = loadSourceFile(refBBPDir);
			bbpSourceRect = bbpSource.getSurface().getRectangle();
			bbpSourceHypo = bbpSource.getHypoLoc();
			if (refAdjustDDW) {
				lines.add("**NOTE: "+refName+" DDW with has been modified for the area to match the Somerville (2006) relationship.**");
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
		}
		lines.add("");
		lines.add("[Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
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
		table.addLine("Section Name", "Area Ruptured", "Patches Ruptured", "Moment", "Equiv. Mag", "Max Slip");
		Map<String, Double> sectAreas = new HashedMap<>();
		Map<String, Integer> sectPatchCounts = new HashMap<>();
		Map<String, Double> sectMoments = new HashMap<>();
		Map<String, Double> sectMaxSlips = new HashMap<>();
		List<SimulatorElement> elems = event.getAllElements();
		double[] elemSlips = event.getAllElementSlips();
		List<FaultSectionPrefData> u3SubSects = catalog.getU3SubSects();
		int subSectOffset = RSQSimUtils.getSubSectIndexOffset(catalog.getElements(), u3SubSects);
		for (int i=0; i<elems.size(); i++) {
			SimulatorElement elem = elems.get(i);
			String fault = u3SubSects.get(elem.getSectionID()-subSectOffset).getParentSectionName();
			Double totArea, totMoment, maxSlip;
			Integer patchCount;
			if (sectAreas.containsKey(fault)) {
				totArea = sectAreas.get(fault);
				totMoment = sectMoments.get(fault);
				patchCount = sectPatchCounts.get(fault);
				maxSlip = sectMaxSlips.get(fault);
			} else {
				totArea = 0d;
				totMoment = 0d;
				patchCount = 0;
				maxSlip = 0d;
			}
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
					RupturePlotGenerator.writeSlipAnimation(event, func, rupAnim, 5);
				else
					RupturePlotGenerator.writeSlipAnimation(event, func, rupAnim, 5, bbpSourceRect[0], bbpSourceRect[1]);
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
			UncertainArbDiscDataset[] gmpeSpectra = null;
			if (spectraGMPEs != null && spectraGMPEs.length > 0 && gmpeRup != null) {
				System.out.print("Calculating GMPEs...");
				gmpeSpectra = new UncertainArbDiscDataset[spectraGMPEs.length];
				for (int i=0; i<gmpeSpectra.length; i++)
					gmpeSpectra[i] = SpectraPlotter.calcGMPE_RotD50(gmpeRup, site, spectraGMPEs[i]);
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
			
			double[] periods = { 1d, 2d, 3d, 4d, 5d, 7.5, 10d };
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
				ShakemapPlotter.plotShakemaps(shakemap, shakemapRegion, shakemapSites, label, resourcesDir,
						prefix, true, shakemapGMPE, gmpeRup, plotPeriods);
			} else {
				System.out.println("Skipping all maps!");
			}
			
			table = MarkdownUtils.tableBuilder();
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
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2829.instance(baseDir);
//		int eventID = 5304;
		int eventID = 31324;
		
		double timeScale = 1d;
		boolean scaleVelocities = false;
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
			vm = RSQSimBBP_Config.VM;
		}
		
		File refShakeMapZip = null;
		for (File dir : refDirs) {
			String name = dir.getName();
			if (dir.isDirectory() && name.contains(catalog.getCatalogDir().getName()) && name.contains(eventID+"")
					&& name.contains("shakemap")) {
				if (timeScale != 1d && !name.contains("-timeScale"+(float)timeScale) && (!scaleVelocities || !name.contains("-velScale")))
					continue;
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
		else
			sites = RSQSimBBP_Config.getCyberShakeInitialLASites();
		
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
		
		RupSpectraPageGen gen = new RupSpectraPageGen(catalog, vmDir, timeScale, scaleVelocities);
		
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
		
		gen.generatePage(event, eventBBPDir, sites, rebuildGIF, rebuildMaps);
		
		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(outputDir);
	}

}
