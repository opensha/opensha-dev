package scratch.kevin.nshm27.figures;

import static scratch.kevin.nshm27.figures.NSHM27_PaperPaths.*;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils.LocationAverager;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.GeoJSONFaultReader;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceCouplingDepthModels;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceDeformationModels;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceDeformationModels.DeformationFront;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceFaultModels;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_LogicTree;

public class SlipProjectionFigures {

	public static void main(String[] args) throws IOException {
		File outputDir = new File(FIGURES_DIR, "slip_projection");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());

		NSHM27_InterfaceFaultModels fm = NSHM27_InterfaceFaultModels.AMSAM_V1;
//		NSHM27_InterfaceFaultModels fm = NSHM27_InterfaceFaultModels.GNMI_V1;
		
		List<? extends FaultSection> sects = fm.buildSubSects(fm);
		
		CPT dipCPT = GMT_CPT_Files.SEQUENTIAL_LAJOLLA_UNIFORM.instance().reverse().rescale(0d, 60d);
		double maxRatio = 1.5;
		CPT ratioCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(0d, 2d).trim(1d, 2d).rescale(1d, maxRatio);
		
		NSHM27_InterfaceDeformationModels dm = fm.getDefaultDeformationModel();
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(sects);
		mapMaker.setWriteGeoJSON(false);
		mapMaker.setFillSurfaces(true);
		
		mapMaker.plotSectScalars(s->s.getAveDip(), dipCPT, "Dip (degrees)");
		mapMaker.plot(outputDir, fm.name()+"_dip", " ");
		
		// cos(dip) = horizontal / on-plane
		// on-plane = horizontal / cos(dip)
		mapMaker.plotSectScalars(s->1d/Math.cos(Math.toRadians(s.getAveDip())),
				ratioCPT, "Projected / Horizontal Slip Rate Ratio");
		mapMaker.plot(outputDir, fm.name()+"_slip_proj_ratio", " ");
		
		LogicTreeBranch<LogicTreeNode> branch = NSHM27_LogicTree.buildDefault(fm.getSeisReg(), TectonicRegionType.SUBDUCTION_INTERFACE, false);
		branch.setValue(NSHM27_InterfaceCouplingDepthModels.NONE);
		
		double maxSlip = maxSlip(dm.apply(fm, branch, sects));
		CPT slipCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, maxSlip);
		
		DeformationFront df = dm.getDeformationFront(fm);
		
		for (NSHM27_InterfaceCouplingDepthModels depthCoupling : NSHM27_InterfaceCouplingDepthModels.values()) {
			branch.setValue(depthCoupling);
			List<? extends FaultSection> dmSects = dm.apply(fm, branch, sects);
			mapMaker.setFaultSections(dmSects);
			mapMaker.plotSectScalars(s->s.getReducedAveSlipRate(),
					slipCPT, "Slip deficit rate (mm/yr)");
			
			List<LocationList> traces = new ArrayList<>();
			List<PlotCurveCharacterstics> traceChars = new ArrayList<>();
			traces.add(df.trace());
			traceChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 10f, Color.BLACK));
			traces.add(df.trace());
			traceChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 10f, Color.LIGHT_GRAY));
			for (int i=0; i<df.slips().length; i++) {
				Location loc = df.trace().get(i);
				Location l1, l2;
				if (i == 0)
					l1 = loc;
				else
					l1 = middle(loc, df.trace().get(i-1));
				if (i == df.slips().length-1)
					l2 = loc;
				else
					l2 = middle(loc, df.trace().get(i+1));
				traces.add(LocationList.of(l1, l2));
				traceChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 5f, slipCPT.getColor(df.slips()[i])));
			}
			mapMaker.plotLines(traces, traceChars);
			
			mapMaker.plot(outputDir, fm.name()+"_"+dm.name()+"_"+depthCoupling.name()+"_slip_deficit_rate",
					dm.getShortName()+" DM, "+depthCoupling.getShortName()+" Taper");
		}
		
		mapMaker.clearLines();
		
		maxSlip = maxSlip(NSHM27_InterfaceDeformationModels.HIGH_COUPLING.apply(fm, branch, sects));
		slipCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, maxSlip);
		
		branch.setValue(NSHM27_InterfaceCouplingDepthModels.AVERAGE);
		for (NSHM27_InterfaceDeformationModels odm : NSHM27_InterfaceDeformationModels.values()) {
			if (odm.getNodeWeight() == 0d)
				continue;
			
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
