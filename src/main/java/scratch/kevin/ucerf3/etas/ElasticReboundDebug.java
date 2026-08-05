package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.Date;

import org.opensha.commons.util.modules.ModuleArchive;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.faultSurface.FaultSection;

import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.ETAS.launcher.ETAS_Config;
import scratch.UCERF3.erf.ETAS.launcher.ETAS_Launcher;
import scratch.UCERF3.erf.ETAS.launcher.TriggerRupture;

public class ElasticReboundDebug {

	public static void main(String[] args) throws IOException {
		ModuleArchive.VERBOSE_DEFAULT = false;
		ETAS_Config config = ETAS_Config.readJSON(new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2026_05_27-FSS_Rupture_201887_M7p8_Start_2026_05_12_1_yr_kCOV_1p5_MaxPtSrcM_6/config.json"));
		
		ETAS_Launcher launcher = new ETAS_Launcher(config, false);
		
		FaultSystemSolution sol = launcher.checkOutFSS();
		FaultSystemRupSet rupSet = sol.getRupSet();
		int mojaveParent = FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Mojave", "S)");
		System.out.println("Mojave parent: "+mojaveParent);
		System.out.println("From Solution");
		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			if (sect.getParentSectionId() == mojaveParent) {
				Date date = new Date(sect.getDateOfLastEvent());
				System.out.println(sect.getSectionId()+". "+sect.getSectionName()+" DOLE: "+date);
			}
		}
		
		FaultSystemSolutionERF erf = (FaultSystemSolutionERF)launcher.checkOutERF();
//		erf = (FaultSystemSolutionERF)launcher.checkOutERF();
		sol = erf.getSolution();
		rupSet = sol.getRupSet();
		System.out.println("From ERF");
		double[] normTS = erf.getNormTimeSinceLastForSections();
		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			if (sect.getParentSectionId() == mojaveParent) {
				Date date = new Date(sect.getDateOfLastEvent());
				System.out.println(sect.getSectionId()+". "+sect.getSectionName()+" DOLE: "+date+"\tnormTS: "+(float)normTS[sect.getSectionId()]);
			}
		}
		
		System.out.println("Start date: "+config.getSimulationStartTimeMillis());
		System.out.println("Triggers:");
		for (TriggerRupture rup : config.getTriggerRuptures()) {
			System.out.println("M"+rup.getMag(rupSet).floatValue()+" at "+rup.getOccurrenceTime(config.getSimulationStartTimeMillis()));
		}
	}

}
