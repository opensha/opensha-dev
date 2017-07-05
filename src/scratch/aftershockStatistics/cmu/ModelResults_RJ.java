package scratch.aftershockStatistics.cmu;

import org.bson.types.ObjectId;
import org.json.simple.JSONObject;
import org.mongodb.morphia.annotations.Entity;
import org.mongodb.morphia.annotations.Id;
import org.opensha.commons.data.siteData.impl.TectonicRegime;
import scratch.aftershockStatistics.*;



@Entity("ModelResults_RJ")
/**
 * <p>Title: ModelResults_RJ </p>
 *
 * <p>Description: Modified implementation of models results</p>
 *
 * <p>Copyright: Copyright (c) 2016</p>
 *
 * <p>Company: CMU</p>
 *
 * @author rewritten by Clark Jeria
 * @version 1.0
 */
public class ModelResults_RJ implements java.io.Serializable{

    private TectonicRegime regime;
    private GenericRJ_Parameters genericParams;
    private ModelResult_RJGeneric genericModel;
    private ModelResult_RJSpecific specificModel;
    private ModelResult_RJBayesian bayesianModel;
    private JSONObject forecast_generic;
    private JSONObject forecast_specific;
    private JSONObject forecast_bayesian;
    public ObjectId eventDataId;
    public String eventId;

    private long sampleTime;
    @Id
    private ObjectId id;

    public ObjectId getId() {
        return id;
    }

    public void setId(ObjectId id) {
        this.id = id;
    }

    private ObjectId dataSourceId;

    public ObjectId getDataSourceId() {
        return dataSourceId;
    }

    public void setDataSourceId(ObjectId dataSourceId) {
        this.dataSourceId = dataSourceId;
    }

    public ModelResults_RJ(){}


    public TectonicRegime getRegime() {
        return regime;
    }

    public void setRegime(TectonicRegime regime) {
        this.regime = regime;
    }

    public GenericRJ_Parameters getGenericParams() {
        return genericParams;
    }

    public void setGenericParams(GenericRJ_Parameters genericParams) {
        this.genericParams = genericParams;
    }

    public ModelResult_RJGeneric getGenericModel() {
        return genericModel;
    }

    public void setGenericModel(ModelResult_RJGeneric genericModel) {
        this.genericModel = genericModel;
    }

    public ModelResult_RJSpecific getSpecificModel() {
        return specificModel;
    }

    public void setSpecificModel(ModelResult_RJSpecific specificModel) {
        this.specificModel = specificModel;
    }

    public ModelResult_RJBayesian getBayesianModel() {
        return bayesianModel;
    }

    public void setBayesianModel(ModelResult_RJBayesian bayesianModel) {
        this.bayesianModel = bayesianModel;
    }

    public long getSampleTime() {
        return sampleTime;
    }

    public void setSampleTime(long sampleTime) {
        this.sampleTime = sampleTime;
    }

    public void setForecast_generic(JSONObject forecast_generic) {
        this.forecast_generic = forecast_generic;
    }

    public JSONObject getForecast_specific() {
        return forecast_specific;
    }

    public void setForecast_specific(JSONObject forecast_specific) {
        this.forecast_specific = forecast_specific;
    }

    public JSONObject getForecast_bayesian() {
        return forecast_bayesian;
    }

    public void setForecast_bayesian(JSONObject forecast_bayesian) {
        this.forecast_bayesian = forecast_bayesian;
    }
}
