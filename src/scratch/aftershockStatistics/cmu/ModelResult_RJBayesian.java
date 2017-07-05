package scratch.aftershockStatistics.cmu;

import scratch.aftershockStatistics.RJ_AftershockModel;
import scratch.aftershockStatistics.RJ_AftershockModel_Bayesian;

/**
 * Created by clark on 11/21/2016.
 */
public class ModelResult_RJBayesian extends RJ_AftershockModel_Bayesian {
    double maxLikelihood_a, maxLikelihood_b, maxLikelihood_c;

    public ModelResult_RJBayesian(){ super(); }
    public ModelResult_RJBayesian(RJ_AftershockModel model1, RJ_AftershockModel model2){
        super(model1, model2);
    }
}
