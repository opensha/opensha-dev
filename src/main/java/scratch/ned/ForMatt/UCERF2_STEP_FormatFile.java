package scratch.ned.ForMatt;

import java.util.ArrayList;

import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.earthquake.calc.ERF_Calculator;

public class UCERF2_STEP_FormatFile {

	public static void main(String[] args) {
		MeanUCERF2 meanUCERF2 = new MeanUCERF2();
		meanUCERF2.setParameter(UCERF2.BACK_SEIS_NAME, UCERF2.BACK_SEIS_INCLUDE);
		meanUCERF2.setParameter(UCERF2.BACK_SEIS_RUP_NAME, UCERF2.BACK_SEIS_RUP_POINT);
		meanUCERF2.setParameter(UCERF2.PROB_MODEL_PARAM_NAME, UCERF2.PROB_MODEL_POISSON);
		meanUCERF2.updateForecast();

		CaliforniaRegions.RELM_GRIDDED gridRELM_Reg  = new CaliforniaRegions.RELM_GRIDDED();
		
		ArrayList<SummedMagFreqDist> mfds = ERF_Calculator.getMagFreqDistsAtLocsInRegion(meanUCERF2, gridRELM_Reg,5.05,36,0.1, true);

		String fileName = "dev/scratch/ned/ForMatt/UCERF2_STEP_FormatFile.txt";
		ERF_Calculator.writeSTEP_FormatFile(mfds, gridRELM_Reg, fileName);
	}

}
