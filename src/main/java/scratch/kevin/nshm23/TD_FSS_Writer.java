package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.util.TimeUtils;
import org.opensha.commons.util.modules.ModuleArchive;
import org.opensha.nshmp.shaded.model.NshmpHazardModel;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.util.MergedSolutionCreator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper.HistoricalRupture;
import org.opensha.sha.faultSurface.FaultSection;

import scratch.ned.nshm23.CEUS_FSS_creator;
import scratch.ned.nshm23.Cascadia_FSS_creator;
import scratch.ned.nshm23.FSS_Fetcher2023;

public class TD_FSS_Writer {

	public static void main(String[] args) throws IOException {
		ModuleArchive.VERBOSE_DEFAULT = false;
		File outputDir = new File("/data/kevin/fss_inversions/2026_10-nshm23-td_erf-solutions");
		FaultSystemSolution wus = getWUS();
		System.out.println(getDOLEstats(wus));
		
		File wusOutput = new File(outputDir, "nshm23-wus-ba-hist_dole.zip");
		System.out.println("Writing WUS to "+wusOutput.getAbsolutePath());
		wus.write(wusOutput);
		
		NshmpHazardModel conusModel = NshmpHazardModel.load(CONUS_DIR.toPath());
		
		FaultSystemSolution cascadia = getCascadia(conusModel);
		System.out.println(getDOLEstats(cascadia));
		File cascadiaOutput = new File(outputDir, "nshm23-cascadia-middle-hist_dole.zip");
		System.out.println("Writing Cascadia to "+cascadiaOutput.getAbsolutePath());
		cascadia.write(cascadiaOutput);
		
		FaultSystemSolution wusPlusCascadia = MergedSolutionCreator.merge(wus, cascadia);
		System.out.println(getDOLEstats(wusPlusCascadia));
		File wusPlusCascadiaOutput = new File(outputDir, "nshm23-wus-cascadia-middle-hist_dole.zip");
		System.out.println("Writing WUS+Cascadia to "+wusPlusCascadiaOutput.getAbsolutePath());
		wusPlusCascadia.write(wusPlusCascadiaOutput);
		
		FaultSystemSolution ceus = getCEUS(conusModel);
		System.out.println(getDOLEstats(ceus));
		File ceusOutput = new File(outputDir, "nshm23-ceus-preferred-hist_dole.zip");
		System.out.println("Writing CEUS to "+ceusOutput.getAbsolutePath());
		ceus.write(ceusOutput);
		
		FaultSystemSolution conusFull = MergedSolutionCreator.merge(wus, ceus, cascadia);
		System.out.println(getDOLEstats(conusFull));
		File conusFullOutput = new File(outputDir, "nshm23-conus-hist_dole.zip");
		System.out.println("Writing WUS+Cascadia to "+conusFullOutput.getAbsolutePath());
		conusFull.write(conusFullOutput);
	}
	
	public static FaultSystemSolution getWUS() throws IOException {
		System.out.println("Loading WUS");
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/data/kevin/fss_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged_gridded_simplified_revised2026.zip"));
		List<HistoricalRupture> histRups = DOLE_SubsectionMapper.loadHistRups();
		System.out.println(DOLE_SubsectionMapper.mapDOLE(sol.getRupSet().getFaultSectionDataList(), histRups, null, null, true));
		return sol;
	}
	
	private static final File MODELS_DIR = new File("/data/kevin/nshm23/nshmp-haz-models");
	private static final File CONUS_DIR = new File(MODELS_DIR, "nshm-conus-6.1.4");
//	private static final File CONUS_DIR = new File(MODELS_DIR, "nshm-conus-6.2.0");
	
	public static FaultSystemSolution getCascadia(NshmpHazardModel model) throws IOException {
		System.out.println("Loading Cascadia");
		FaultSystemSolution sol = Cascadia_FSS_creator.getFaultSystemSolution(CONUS_DIR, model, Cascadia_FSS_creator.FaultModelEnum.MIDDLE);
		List<HistoricalRupture> histRups = DOLE_SubsectionMapper.loadHistRups( );
		System.out.println(DOLE_SubsectionMapper.mapDOLE(sol.getRupSet().getFaultSectionDataList(), histRups, null, null, true));
		return sol;
	}
	
	public static FaultSystemSolution getCEUS(NshmpHazardModel model) throws IOException {
		System.out.println("Loading CEUS");
		ArrayList<FaultSystemSolution> sols = CEUS_FSS_creator.getFaultSystemSolutionList(CONUS_DIR, model, CEUS_FSS_creator.FaultModelEnum.PREFERRED);
		FaultSystemSolution sol = sols.size() == 1 ? sols.get(0) : MergedSolutionCreator.merge(sols);
		List<HistoricalRupture> histRups = DOLE_SubsectionMapper.loadHistRups();
		System.out.println(DOLE_SubsectionMapper.mapDOLE(sol.getRupSet().getFaultSectionDataList(), histRups, null, null, true));
		return sol;
	}
	
	private static String getDOLEstats(FaultSystemSolution sol) {
		int numWith = 0;
		int totNum = 0;
		int earliestYear = Integer.MAX_VALUE;
		int latestYear = Integer.MIN_VALUE;
		for (FaultSection sect : sol.getRupSet().getFaultSectionDataList()) {
			totNum++;
			long dole = sect.getDateOfLastEvent();
			if (dole > Long.MIN_VALUE) {
				numWith++;
				int year = TimeUtils.epochMillisToYear(dole);
				earliestYear = Integer.min(earliestYear, year);
				latestYear = Integer.max(latestYear, year);
			}
		}
		
		return numWith+"/"+totNum+" have DOLE; years in ["+earliestYear+", "+latestYear+"]";
	}

}
