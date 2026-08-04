package scratch.kevin.nshm23.hazardValidation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.util.Precision;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.nshmp.model.HazardModel;
import gov.usgs.earthquake.nshmp.model.NshmErf;
import gov.usgs.earthquake.nshmp.model.NshmSource;
import gov.usgs.earthquake.nshmp.model.SystemRuptureSet.SystemRupture;

public class WrapperRupSetRupMapper {

	public static void main(String[] args) throws IOException {
		File inputSolFile = new File("/home/kevin/OpenSHA/fss_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/"
//				+ "results_WUS_FM_v3_branch_averaged_gridded_simplified_revised2026.zip");
				+ "results_WUS_FM_v3_branch_averaged_gridded_simplified_matchPrev.zip");
		File modelDir = new File("/data/kevin/nshm23/nshmp-haz-models/nshm-conus-6.1.3");
//		File modelDir = new File("/data/kevin/nshm23/nshmp-haz-models/nshm-conus-6.2.0");
		
		FaultSystemSolution sol = FaultSystemSolution.load(inputSolFile);
		
		HazardModel model = HazardModel.load(modelDir.toPath());
		
		NshmErf wrapper = new NshmErf(model, Set.of(TectonicRegionType.ACTIVE_SHALLOW), IncludeBackgroundOption.EXCLUDE);
		wrapper.getTimeSpan().setDuration(1d);
		wrapper.updateForecast();
		
		map(sol, wrapper);
	}
	
	private static final double MAG_TOL = 0.01;
	private static final double RATE_ABS_TOL = 1e-6;
	private static final double RATE_REL_TOL = 1e-4;
	
	public static WrapperMatch[] map(FaultSystemSolution sol, NshmErf erf) {
		Preconditions.checkState(erf.getTimeSpan().getDuration() == 1d, "Duration must be 1 year");
		List<WrapperMatch> rups = new ArrayList<>();
		for (int sourceID=0; sourceID<erf.getNumSources(); sourceID++) {
			ProbEqkSource source = erf.getSource(sourceID);
			if (source instanceof NshmSource.System) {
				NshmSource.System systemSource = (NshmSource.System)source;
				Preconditions.checkState(systemSource.getNumRuptures() == 1);
				ProbEqkRupture rup = systemSource.getRupture(0);
				rups.add(new WrapperMatch(-1, sourceID, systemSource, systemSource.delegate(), rup));
			}
		}
		
		SystemRupture del0 = rups.get(0).wrapperSource.delegate();
		System.out.println("Wrapper rup0:\t"+del0.magnitude()+"\t"+del0.rate()+"\t"+del0.rake());
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		System.out.println("RupSet rup0:\t"+(float)rupSet.getMagForRup(0)+"\t"+(float)sol.getRateForRup(0)+"\t"+(float)rupSet.getAveRakeForRup(0));
		System.out.println("Building mappings between "+rupSet.getNumRuptures()+" FSS ruptures and "+rups.size()+" wrapper ruptures");
		Preconditions.checkState(rups.size() <= rupSet.getNumRuptures());
		WrapperMatch[] mappings = new WrapperMatch[rupSet.getNumRuptures()];
		
		int curWrapperIndex = 0;
		for (int rupIndex=0; rupIndex<mappings.length; rupIndex++) {
			double rate = sol.getRateForRup(rupIndex);
			if (rate == 0d)
				continue;
			String rupStr = rupString(sol, rupIndex);
			double mag = rupSet.getMagForRup(rupIndex);
			WrapperMatch wrapperRup = rups.get(curWrapperIndex);
			SystemRupture delagate = wrapperRup.wrapperSource.delegate();
			String delStr = rupString(delagate);
			Preconditions.checkState(Precision.equals(mag, delagate.magnitude(), MAG_TOL),
					"Magnitude mismatch:\n\tFSS %s:\t%s\n\tWrapper %s:\t%s",
					rupIndex, rupStr, curWrapperIndex, delStr);
			Preconditions.checkState(Precision.equals(rate, delagate.rate(), RATE_ABS_TOL),
					"Abs rate mismatch:\n\tFSS %s:\t%s\n\tWrapper %s:\t%s",
					rupIndex, rupStr, curWrapperIndex, delStr);
			Preconditions.checkState(Math.abs(rate - delagate.rate())/rate < RATE_REL_TOL,
					"Relative rate mismatch:\n\tFSS %s:\t%s\n\tWrapper %s:\t%s",
					rupIndex, rupStr, curWrapperIndex, delStr);
			mappings[rupIndex] = new WrapperMatch(rupIndex, wrapperRup.wrapperSourceIndex, wrapperRup.wrapperSource,
					wrapperRup.delagate, wrapperRup.wrapperRup);
			
			curWrapperIndex++;
		}
		Preconditions.checkState(curWrapperIndex == rups.size());
		System.out.println("DONE; matched "+curWrapperIndex+" ruptures.");
		
		return mappings;
	}
	
	private static String rupString(SystemRupture rup) {
		return "M="+(float)rup.magnitude()+"\trake="+(float)rup.rake()+"\trate="+(float)rup.rate();
	}
	
	private static String rupString(FaultSystemSolution sol, int rupIndex) {
		return "M="+(float)sol.getRupSet().getMagForRup(rupIndex)
				+"\trake="+(float)sol.getRupSet().getAveRakeForRup(rupIndex)+"\trate="+(float)sol.getRateForRup(rupIndex);
	}
	
	public record WrapperMatch(int fssRupIndex, int wrapperSourceIndex, NshmSource.System wrapperSource,
			SystemRupture delagate, ProbEqkRupture wrapperRup) {};

}
