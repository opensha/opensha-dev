package scratch.kevin.nshm27.figures;

import static scratch.kevin.nshm27.figures.NSHM27_PaperPaths.*;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationUtils.LocationAverager;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.FaultUtils;
import org.opensha.commons.util.FaultUtils.AngleAverager;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceCouplingDepthModels;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceDeformationModels;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceDeformationModels.DeformationFront;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceFaultModels;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_LogicTree;
import net.mahdilamb.colormap.Colors;

public class SlipProjectionFigures {

	public static void main(String[] args) throws IOException {
//		File outputDir = new File(FIGURES_DIR, "slip_projection");
		File outputDir = new File("/tmp/nshm27_slip_projection");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());

		NSHM27_InterfaceFaultModels fm = NSHM27_InterfaceFaultModels.AMSAM_V1;
		double maxSlip = 120;
//		NSHM27_InterfaceFaultModels fm = NSHM27_InterfaceFaultModels.GNMI_V1;
//		double maxSlip = 20d;
		
		List<? extends FaultSection> sects = fm.buildSubSects(fm);
		
		CPT dipCPT = GMT_CPT_Files.SEQUENTIAL_LAJOLLA_UNIFORM.instance().reverse().rescale(0d, 60d);
		double maxRatio = 1.5;
		CPT ratioCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(0d, 2d).trim(1d, 2d).rescale(1d, maxRatio);
		
		NSHM27_InterfaceDeformationModels.Aggregated dm = fm.getDefaultDeformationModel();
		
		LogicTreeBranch<LogicTreeNode> branch = NSHM27_LogicTree.buildDefault(fm.getSeismicityRegion(), TectonicRegionType.SUBDUCTION_INTERFACE, false);
		branch.setValue(NSHM27_InterfaceCouplingDepthModels.NONE);
		
//		double maxSlip = maxSlip(dm.apply(fm, branch, sects));
		CPT slipCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, maxSlip);
		
		DeformationFront df = NSHM27_InterfaceDeformationModels.getDeformationFront(fm);
		double moveOffset = 15d; // km
		
		LocationList dfTrace = df.trace();
		List<FaultSection> regSects = new ArrayList<>(sects);
		if (moveOffset > 0d) {
			LocationList movedTrace = new LocationList(dfTrace.size());
			for (int i=0; i<dfTrace.size(); i++) {
				Location loc = dfTrace.get(i);
				double az;
				if (i == 0) {
					az = LocationUtils.azimuth(loc, dfTrace.get(i+1));
				} else if (i == dfTrace.size()-1) {
					az = LocationUtils.azimuth(dfTrace.get(i-1), loc);
				} else {
					Location before = dfTrace.get(i-1);
					Location after = dfTrace.get(i+1);
					double az1 = LocationUtils.azimuth(before, loc);
					double az2 = LocationUtils.azimuth(loc, after);
					AngleAverager avg = new AngleAverager();
					avg.add(az1, LocationUtils.horzDistanceFast(before, loc));
					avg.add(az2, LocationUtils.horzDistanceFast(loc, after));
					az = avg.getAverage();
				}
				movedTrace.add(LocationUtils.location(loc, Math.toRadians(az-90d), moveOffset));
			}
			dfTrace = movedTrace;
			FaultTrace offsetTrace = new FaultTrace();
			offsetTrace.addAll(movedTrace);
			FaultSection offsetSect = new GeoJSONFaultSection.Builder(999, "Offset", offsetTrace).dip(90d).upperDepth(0d).lowerDepth(10d).build();
			regSects.add(offsetSect);
		}
		
		GeographicMapMaker.buildBufferedRegion(regSects);
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(GeographicMapMaker.buildBufferedRegion(regSects));
		mapMaker.setFaultSections(sects);
		mapMaker.setWriteGeoJSON(false);
		mapMaker.setFillSurfaces(true);
		
		mapMaker.plotSectScalars(s->s.getAveDip(), dipCPT, "Dip (degrees)");
		mapMaker.plot(outputDir, fm.name()+"_dip", " ");
		
		// cos(dip) = horizontal / on-plane
		// on-plane = horizontal / cos(dip)
		mapMaker.plotSectScalars(s->1d/Math.cos(Math.toRadians(s.getAveDip())),
				ratioCPT, "Projected / Horizontal Slip Rate Ratio");
		mapMaker.plot(outputDir, fm.name()+"_slip_proj_ratio", " ");
		
		for (NSHM27_InterfaceCouplingDepthModels depthCoupling : NSHM27_InterfaceCouplingDepthModels.values()) {
			branch.setValue(depthCoupling);
			List<? extends FaultSection> dmSects = dm.apply(fm, branch, sects);
			mapMaker.setFaultSections(dmSects);
			mapMaker.plotSectScalars(s->s.getReducedAveSlipRate(),
					slipCPT, "Slip deficit rate (mm/yr)");
			
			List<LocationList> traces = new ArrayList<>();
			List<PlotCurveCharacterstics> traceChars = new ArrayList<>();
			
			traces.add(dfTrace);
			traceChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 15f, Color.BLACK));
//			traces.add(dfTrace);
//			traceChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 10f, Color.LIGHT_GRAY));
			int num = dfTrace.size();
			double[] slips = NSHM27_InterfaceDeformationModels.getCoupledSlipRates(df, dm.getSampler());
			for (int i=0; i<num; i++) {
				Location loc = dfTrace.get(i);
				Location l1, l2;
				if (i == 0)
					l1 = loc;
				else
					l1 = middle(loc, dfTrace.get(i-1));
				if (i == num-1)
					l2 = loc;
				else
					l2 = middle(loc, dfTrace.get(i+1));
				traces.add(LocationList.of(l1, l2));
				traceChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 10f, slipCPT.getColor(slips[i])));
			}
			mapMaker.plotLines(traces, traceChars);
			
			mapMaker.plot(outputDir, fm.name()+"_"+dm.getFilePrefix()+"_"+depthCoupling.name()+"_slip_deficit_rate",
					dm.getShortName()+" DM, "+depthCoupling.getShortName()+" Taper");
			
			if (depthCoupling == NSHM27_InterfaceCouplingDepthModels.NONE) {
				// add length annotations
				
				FaultSection sect= fm.getFaultSections().get(0);
				FaultTrace upper = sect.getFaultTrace();
				FaultTrace lower = sect.getLowerFaultTrace();
				
				int numResample = 5000;
				
				upper = FaultUtils.resampleTrace(upper, numResample);
				lower = FaultUtils.resampleTrace(lower, numResample);
				
				FaultTrace middle = new FaultTrace();
				for (int i=0; i<numResample; i++) {
					Location upperLoc = upper.get(i);
					Location lowerLoc = lower.get(i);
					double dist = LocationUtils.horzDistance(lowerLoc, upperLoc);
					double az = LocationUtils.azimuthRad(lowerLoc, upperLoc);
					Location middleLoc = LocationUtils.location(lowerLoc, az, 0.5*dist);
					middle.add(middleLoc);
				}
				
				double markerDelta = 50d;
				double nextMarker = 0d;
				
				LocationList lengthLocs = new LocationList();
				List<Double> lengthVals = new ArrayList<>();
				
				double curLen = 0d;
				for (int i=0; i<middle.size(); i++) {
					if (i > 0)
						curLen += LocationUtils.horzDistance(middle.get(i-1), middle.get(i));
					if (curLen >= nextMarker || i == numResample-1) {
						Location loc;
						double markerLen;
						if (i > 0 && i < numResample-1) {
							double overshoot = curLen - nextMarker;
							double backAz = LocationUtils.azimuthRad(middle.get(i), middle.get(i-1));
							loc = LocationUtils.location(middle.get(i), backAz, overshoot);
							markerLen = nextMarker;
						} else {
							loc = middle.get(i);
							markerLen = curLen;
						}
						lengthLocs.add(loc);
						lengthVals.add(markerLen);
						
						nextMarker += markerDelta;
					}
				}
				
				for (int i=0; i<lengthLocs.size(); i++) {
					Location loc = lengthLocs.get(i);
					double markerLen = lengthVals.get(i);
					
					XYTextAnnotation labelAnn = new XYTextAnnotation(" "+(int)Math.round(markerLen), loc.lon, loc.lat);
					labelAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 24));
					labelAnn.setTextAnchor(TextAnchor.CENTER_LEFT);
					labelAnn.setRotationAnchor(TextAnchor.CENTER_LEFT);
					double azimuth;
					if (i == 0) {
						azimuth = LocationUtils.azimuthRad(lengthLocs.get(0), lengthLocs.get(1));
					} else if (i == lengthLocs.size()-1) {
						azimuth = LocationUtils.azimuthRad(lengthLocs.get(lengthLocs.size()-2), lengthLocs.get(lengthLocs.size()-1));
					} else {
						AngleAverager avg = new AngleAverager();
						avg.add(LocationUtils.azimuth(lengthLocs.get(i-1), lengthLocs.get(i)), 1d);
						avg.add(LocationUtils.azimuth(lengthLocs.get(i), lengthLocs.get(i+1)), 1d);
						azimuth = Math.toRadians(avg.getAverage());
					}
					double rotAngle = azimuth;
					// correct for aspect ratio
					rotAngle = mapMaker.getRotationAngleCorrectedForAspectRatio(rotAngle);
					labelAnn.setRotationAngle(rotAngle);
//					labelAnn.setPaint(Colors.tab_lightbrown);
//					labelAnn.setPaint(Colors.tab_brown);
					labelAnn.setPaint(Colors.tab_brown.darker().darker().darker());
					mapMaker.addAnnotation(labelAnn);
				}
				
				mapMaker.setScatterSymbol(PlotSymbol.FILLED_CIRCLE, 6f, PlotSymbol.CIRCLE, Color.BLACK);
				mapMaker.plotScatters(lengthLocs, Color.WHITE);
				
				mapMaker.plot(outputDir, fm.name()+"_"+dm.getFilePrefix()+"_slip_and_length",
						dm.getShortName()+" DM & Kilometer Length Markers");
				
				mapMaker.clearScatters();
				mapMaker.clearAnnotations();
			}
		}
		
		mapMaker.clearLines();
		
		maxSlip = maxSlip(NSHM27_InterfaceDeformationModels.Aggregated.HIGH_COUPLING.apply(fm, branch, sects));
		slipCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, maxSlip);
		
		branch.setValue(NSHM27_InterfaceCouplingDepthModels.AVERAGE);
		for (NSHM27_InterfaceDeformationModels.Aggregated odm : NSHM27_InterfaceDeformationModels.Aggregated.values()) {
			
			branch.setValue(odm);
			List<? extends FaultSection> dmSects = odm.apply(fm, branch, sects);
			
			mapMaker.setFaultSections(dmSects);
			mapMaker.plotSectScalars(s->s.getReducedAveSlipRate(),
					slipCPT, "Slip deficit rate (mm/yr)");
			
			mapMaker.plot(outputDir, fm.name()+"_"+odm.name()+"_slip_deficit_rate",
					odm.getShortName()+" DM, Average Taper");
		}
	}
	
	private static double maxSlip(List<? extends FaultSection> sects) {
		double maxSlip = sects.stream().mapToDouble(S->S.getOrigAveSlipRate()).max().getAsDouble();
		double div10 = maxSlip / 10d;
		if (div10 % 1 < 0.4)
			return Math.floor(div10)*10d;
		return Math.ceil(div10)*10d;
	}
	
	private static Location middle(Location l1, Location l2) {
		return new LocationAverager().add(l1, 1d).add(l2, 1d).getAverage();
	}

}
