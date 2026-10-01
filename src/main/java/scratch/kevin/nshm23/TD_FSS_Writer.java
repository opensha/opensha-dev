package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.util.TimeUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.timeDependence.DOLE_SubsectionMapper.HistoricalRupture;
import org.opensha.sha.faultSurface.FaultSection;

public class TD_FSS_Writer {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/data/kevin/fss_inversions/2026_10-nshm23-td_erf-solutions");
		FaultSystemSolution wus = getWUS();
		System.out.println(getDOLEstats(wus));
		
		wus.write(new File(outputDir, "nshm23-wus-ba-hist_dole.zip"));
	}
	
	public static FaultSystemSolution getWUS() throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/data/kevin/fss_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged_gridded_simplified_revised2026.zip"));
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
