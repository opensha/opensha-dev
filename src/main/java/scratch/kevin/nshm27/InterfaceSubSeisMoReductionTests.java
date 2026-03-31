package scratch.kevin.nshm27;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationUtils.LocationAverager;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM26_DeclusteringAlgorithms;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM26_InterfaceDeformationModels;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM26_InterfaceFaultModels;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM26_LogicTree;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM26_SeisRateModelBranch;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM26_SeisSmoothingAlgorithms;
import gov.usgs.earthquake.nshmp.erf.nshm27.util.InterfaceGridAssociations;
import gov.usgs.earthquake.nshmp.erf.nshm27.util.NSHM26_SeisPDF_Loader;
import gov.usgs.earthquake.nshmp.erf.nshm27.util.NSHM26_RegionLoader.NSHM26_SeismicityRegions;

public class InterfaceSubSeisMoReductionTests {

	public static void main(String[] args) throws IOException {
		LogicTreeBranch<LogicTreeNode> branch = NSHM26_LogicTree.buildDefault(
//				NSHM26_SeismicityRegions.AMSAM, TectonicRegionType.SUBDUCTION_INTERFACE, false);
				NSHM26_SeismicityRegions.GNMI, TectonicRegionType.SUBDUCTION_INTERFACE, false);
		NSHM26_InterfaceFaultModels fm = branch.requireValue(NSHM26_InterfaceFaultModels.class);
		NSHM26_InterfaceDeformationModels dm = branch.requireValue(NSHM26_InterfaceDeformationModels.class);

		NSHM26_SeismicityRegions reg = fm.getSeisReg();
		File pdfBaseDir = new File("/home/kevin/OpenSHA/nshm26/data/spatial_seis_pdfs/"+reg.name().toLowerCase()+"/2026_03_09-v1_2D/INTERFACE");
		
		System.out.println("Branch: "+branch+"; reg="+reg);
		
		NSHM26_SeisRateModelBranch rateModel = NSHM26_SeisRateModelBranch.PREFFERRED;
		NSHM26_DeclusteringAlgorithms decluster = NSHM26_DeclusteringAlgorithms.AVERAGE;
		NSHM26_SeisSmoothingAlgorithms smooth = NSHM26_SeisSmoothingAlgorithms.AVERAGE;
		double cutoffHorzDist = 50d;
		
		List<? extends FaultSection> sects = dm.build(branch);
		Location[] middles = new Location[sects.size()];
		for (int i=0; i<middles.length; i++) {
			FaultSection sect = sects.get(i);
			FaultTrace upper = sect.getFaultTrace();
			FaultTrace lower = sect.getLowerFaultTrace();
			middles[i] = new LocationAverager().add(upper.first(), 1d).add(upper.last(), 1d).add(lower.first(), 1d).add(lower.last(), 1d).getAverage();
		}
		
		RuptureSurface[] sectSurfs = new RuptureSurface[sects.size()];
		Region[] sectOutlines = new Region[sects.size()];
		for (int s=0; s<sects.size(); s++) {
			sectSurfs[s] = sects.get(s).getFaultSurface(2d);
			sectOutlines[s] = new Region(sectSurfs[s].getPerimeter(), BorderType.MERCATOR_LINEAR);
		}
		
		GriddedGeoDataSet pdf = NSHM26_SeisPDF_Loader.load2D(pdfBaseDir, reg, decluster, smooth);
		
		InterfaceGridAssociations assoc = new InterfaceGridAssociations(sects, pdf.getRegion());
		double sumMapped = 0d;
		for (int i=0; i<pdf.size(); i++)
			sumMapped += pdf.get(i)*assoc.getNodeFraction(i);
		System.out.println("PDF fraction mapped: "+(float)sumMapped);
		
		double[] moments = new double[sects.size()];
		for (int s=0; s<moments.length; s++)
			moments[s] = sects.get(s).calcMomentRate(true);
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(5.01, 9.51);
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(sects);
		CPT cpt = GMT_CPT_Files.SEQUENTIAL_LAJOLLA_UNIFORM.instance().reverse().rescale(0d, 1d);
		
//		CPT pdfCPT = new CPT(0d, pdf.getMaxZ(), new Color(255, 255, 255, 0), new Color(0, 0, 0, 127));
		CPT pdfCPT = GMT_CPT_Files.SEQUENTIAL_OSLO_UNIFORM.instance().reverse().rescale(0d, pdf.getMaxZ());
		
		mapMaker.plotXYZData(pdf, pdfCPT, "Average PDF");
		Color transColor = new Color(255, 255, 255, 0);
		List<Color> transColors = new ArrayList<>(sects.size());
		for (int i=0; i<sects.size(); i++)
			transColors.add(transColor);
		mapMaker.plotSectColors(transColors);
		
		mapMaker.plot(new File("/tmp"), "mo_red_"+reg.name()+"_pdf", " ");
		
		double[] mMins = {6.55, 7.05, 7.55};
		for (double mMin : mMins) {
			IncrementalMagFreqDist seisMFD = rateModel.build(reg, TectonicRegionType.SUBDUCTION_INTERFACE, refMFD,
					refMFD.getX(refMFD.getClosestXIndex(mMin-0.1)));
//			System.out.println("Seis MFD for interface Mmin="+mMin+":\n"+seisMFD);
			double[] impliedMoments = new double[sects.size()];
			for (int i=0; i<pdf.size(); i++) {
				double moSum = 0d;
				for (int j=0; j<seisMFD.size(); j++)
					moSum += MagUtils.magToMoment(seisMFD.getX(j))*seisMFD.getY(j)*pdf.get(i);
				Map<Integer, Double> mappings = assoc.getSectionFracsOnNode(i);
				for (int s : mappings.keySet())
					impliedMoments[s] += moSum*mappings.get(s);
			}
			
			double[] moFracts = new double[sects.size()];
			MinMaxAveTracker track = new MinMaxAveTracker();
			for (int i=0; i<sects.size(); i++) {
				double fract = impliedMoments[i] / moments[i];
				track.addValue(fract);
				moFracts[i] = fract;
			}
			System.out.println("Mmin="+mMin+" sub-seis moment fracts:\t"+track);
			mapMaker.plotSectScalars(moFracts, cpt, "Implied Sub-Seis Mo Red, Mmin="+(float)mMin);
			mapMaker.plotXYZData(pdf, pdfCPT, null);
			mapMaker.plot(new File("/tmp"), "mo_red_"+reg.name()+"_mMin"+(float)mMin, " ");
		}
	}

}
