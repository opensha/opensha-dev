package scratch.kevin;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.simulators.stiffness.AggregatedStiffnessCalculator;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator;
import org.opensha.sha.simulators.stiffness.AggregatedStiffnessCalculator.AggregationMethod;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator.PatchAlignment;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator.StiffnessType;

import com.google.gson.annotations.Expose;

import net.mahdilamb.colormap.Colors;
import scratch.UCERF3.analysis.GMT_CA_Maps;

public class ScenarioStiffnessCalc {

	public static void main(String[] args) throws IOException {
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(new File("/home/kevin/git/ucerf3-etas-launcher/inputs/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip"));
		int rupIndex = 201887;
		float cffRatioThresh = 0.5f;
		double stiffGridSpacing = 2d;
		double coeffOfFriction = 0.5;
		SubSectStiffnessCalculator stiffnessCalc = new SubSectStiffnessCalculator(
				rupSet.getFaultSectionDataList(), stiffGridSpacing, 3e4, 3e4, coeffOfFriction, PatchAlignment.FILL_OVERLAP, 1d);
		AggregatedStiffnessCalculator sumAgg = new AggregatedStiffnessCalculator(StiffnessType.CFF, stiffnessCalc, true,
				AggregationMethod.FLATTEN, AggregationMethod.SUM, AggregationMethod.SUM, AggregationMethod.SUM);
		
		List<FaultSection> sources = rupSet.getFaultSectionDataForRupture(rupIndex);
		
//		Region region = new CaliforniaRegions.RELM_SOCAL();
		Region region = new Region(new Location(32, -120), new Location(36, -114));
		
		CPT cffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-1d, 1d);
		
		List<Color> colors = new ArrayList<>(rupSet.getNumSections());
		double[] sectsInRegion = rupSet.getFractSectsInsideRegion(region, false);
		
		List<Double> scalars = new ArrayList<>();
		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			Color color;
			if (sources.contains(sect)) {
				color = Colors.tab_green;
				scalars.add(Double.POSITIVE_INFINITY);
			} else if (sectsInRegion[sect.getSectionId()] > 0d) {
				double cff = sumAgg.calc(sources, sect);
				color = cffCPT.getColor(cff);
				System.out.println(sect.getSectionId()+". "+sect.getSectionName()+":\t"+cff);
				scalars.add(cff);
			} else {
				color = Colors.LIGHT_GRAY;
				scalars.add(Double.NEGATIVE_INFINITY);
			}
			colors.add(color);
		}
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(region);
		mapMaker.setFaultSections(rupSet.getFaultSectionDataList());
		mapMaker.plotSectColors(colors, cffCPT, "ΔCFF (MPa)", scalars);
		
		mapMaker.plot(new File("/tmp"), "u3_cff_"+rupIndex, "Rupture "+rupIndex+" Coulomb Stress Change");
		
		System.out.println("DONE");
		System.exit(0);
	}

}
