package scratch.aftershockStatistics.cmu;

import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import scratch.aftershockStatistics.RJ_AftershockModel_SequenceSpecific;

/**
 * Created by clark on 11/21/2016.
 */
public class ModelResult_RJSpecific extends RJ_AftershockModel_SequenceSpecific {
    double maxLikelihood_a, maxLikelihood_b, maxLikelihood_c;

    public ModelResult_RJSpecific(){
        super();
    }
    public ModelResult_RJSpecific(ObsEqkRupture mainShock, ObsEqkRupList aftershockList,
                                  double magCat, double capG, double capH,
                                  double b, double dataStartTimeDays, double dataEndTimeDays,
                                  double min_a, double max_a, int num_a,
                                  double min_p, double max_p, int num_p,
                                  double min_c, double max_c, int num_c){
        super(mainShock, aftershockList,
        magCat, capG, capH,
        b, dataStartTimeDays, dataEndTimeDays,
        min_a, max_a, num_a,
        min_p, max_p, num_p,
        min_c, max_c, num_c);
    }
}
