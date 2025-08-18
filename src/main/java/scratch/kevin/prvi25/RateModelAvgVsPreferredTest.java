package scratch.kevin.prvi25;

import java.io.IOException;

import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateModel;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeismicityRateEpoch;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionCaribbeanSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionMuertosSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;

import com.google.common.base.Preconditions;

public class RateModelAvgVsPreferredTest {

	public static void main(String[] args) throws IOException {
		double prefWeight = PRVI25_CrustalSeismicityRate.PREFFERRED.getNodeWeight(null);
		double lowerWeight = PRVI25_CrustalSeismicityRate.LOW.getNodeWeight(null);
		double upperWeight = PRVI25_CrustalSeismicityRate.HIGH.getNodeWeight(null);
		Preconditions.checkState((float)(prefWeight + lowerWeight + upperWeight) == 1f);
		for (PRVI25_SeismicityRegions seisReg : PRVI25_SeismicityRegions.values()) {
			System.out.println(seisReg.getName());
			double sumAvg = 0d;
			double sumPref = 0d;
			for (PRVI25_SeismicityRateEpoch epoch : PRVI25_SeismicityRateEpoch.values()) {
				System.out.println("\t"+epoch.getName());
				SeismicityRateModel model;
				if (seisReg == PRVI25_SeismicityRegions.CRUSTAL)
					model = PRVI25_CrustalSeismicityRate.loadRateModel(epoch);
				else if (seisReg == PRVI25_SeismicityRegions.CAR_INTERFACE)
					model = PRVI25_SubductionCaribbeanSeismicityRate.loadRateModel(epoch, false);
				else if (seisReg == PRVI25_SeismicityRegions.CAR_INTRASLAB)
					model = PRVI25_SubductionCaribbeanSeismicityRate.loadRateModel(epoch, true);
				else if (seisReg == PRVI25_SeismicityRegions.MUE_INTERFACE)
					model = PRVI25_SubductionMuertosSeismicityRate.loadRateModel(epoch, false);
				else if (seisReg == PRVI25_SeismicityRegions.MUE_INTRASLAB)
					model = PRVI25_SubductionMuertosSeismicityRate.loadRateModel(epoch, true);
				else
					throw new IllegalStateException();
				
				double meanRate = model.getMeanRecord().rateAboveM1;
				double lowerRate = model.getLowerRecord().rateAboveM1;
				double upperRate = model.getUpperRecord().rateAboveM1;
				double avgRate = prefWeight*meanRate + lowerWeight*lowerRate + upperWeight*upperRate;
				
				sumPref += epoch.getNodeWeight(null)*meanRate;
				sumAvg += epoch.getNodeWeight(null)*avgRate;
				
				System.out.println("\t\tPREF: "+(float)meanRate);
				System.out.println("\t\tAVG: "+(float)avgRate+"\t"+(avgRate > meanRate ? "ABOVE" : "BELOW"));
			}
			System.out.println("\tPREF: "+(float)sumPref);
			System.out.println("\tAVG: "+(float)sumAvg+"\t"+(sumAvg > sumPref ? "ABOVE" : "BELOW"));
		}
	}

}
