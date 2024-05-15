package scratch.kevin.nshm23.prvi;

import java.io.IOException;
import java.util.List;

import org.opensha.commons.geo.json.Feature;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets.CoulombRupSetConfig;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.util.SubSectionBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalFaultModels;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;

import com.google.common.base.Preconditions;

public class BunceMainRidgeConnectionTests {

	public static void main(String[] args) throws IOException {
		PRVI25_CrustalFaultModels fm = PRVI25_CrustalFaultModels.PRVI_FM_INITIAL;
		List<? extends FaultSection> allSects = fm.getFaultSections();
		
		FaultSection bunce5 = null;
		FaultSection mainRidge1 = null;
		
		for (FaultSection sect : allSects) {
			if (sect.getSectionName().equals("Bunce 5"))
				bunce5 = sect;
			if (sect.getSectionName().equals("Main Ridge 1"))
				mainRidge1 = sect;
		}
		
		Preconditions.checkNotNull(bunce5);
		Preconditions.checkNotNull(mainRidge1);
		
		// this doesn't
		bunce5.setAveRake(0d);
		
		// this does the trick
//		bunce5.setAveRake(90d);
//		bunce5 = setDip(bunce5, 60);
		
		// or this does the trick
		mainRidge1.setAveRake(-45d);
		
		List<FaultSection> subSects = SubSectionBuilder.buildSubSects(List.of(bunce5, mainRidge1));
		
		CoulombRupSetConfig config = new RuptureSets.CoulombRupSetConfig(subSects, "coulomb", NSHM23_ScalingRelationships.AVERAGE);
		FaultSystemRupSet rupSet = config.build(1);
		ClusterRuptures cRups = rupSet.requireModule(ClusterRuptures.class);
		
		boolean hasJump = false;
		double largestMag = 0d;
		double largestMultiMag = 0d;
		int mostSubSects = 0;
		for (int r=0; r<cRups.size(); r++) {
			double mag = rupSet.getMagForRup(r);
			mostSubSects = Integer.max(mostSubSects, rupSet.getSectionsIndicesForRup(r).size());
			largestMag = Double.max(mag, largestMag);
			ClusterRupture cRup = cRups.get(r);
			if (cRup.getTotalNumJumps() > 0) {
				hasJump = true;
				largestMultiMag = Double.max(mag, largestMultiMag);
			}
		}

		System.out.println("Bunce5: rake="+(int)bunce5.getAveRake()+", dip="+(int)bunce5.getAveDip()+", dipDir="+(int)bunce5.getDipDirection());
		System.out.println("MainRidge1: rake="+(int)mainRidge1.getAveRake()+", dip="+(int)mainRidge1.getAveDip()+", dipDir="+(int)mainRidge1.getDipDirection());
		System.out.println("Largest mag: "+largestMag);
		System.out.println("Largest subsect count: "+mostSubSects+"/"+subSects.size());
		System.out.println("Has jump? "+hasJump);
		if (hasJump)
			System.out.println("Largest jumping mag: "+largestMultiMag);
	}
	
	private static FaultSection setDip(FaultSection sect, double dip) {
		Feature feature = new GeoJSONFaultSection(sect).toFeature();
		feature.properties.set(GeoJSONFaultSection.DIP, dip);
		return GeoJSONFaultSection.fromFeature(feature);
	}

}
