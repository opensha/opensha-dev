package scratch.aftershockStatistics.cmu;

import scratch.aftershockStatistics.GenericRJ_Parameters;
import scratch.aftershockStatistics.RJ_AftershockModel_Generic;

/**
 * Created by clark on 11/21/2016.
 */
public class ModelResult_RJGeneric extends RJ_AftershockModel_Generic {
    double maxLikelihood_a, maxLikelihood_b, maxLikelihood_c;

    public ModelResult_RJGeneric(){ super();}
    public ModelResult_RJGeneric(double mainMag, GenericRJ_Parameters b){
        super(mainMag, b);
    }
}
