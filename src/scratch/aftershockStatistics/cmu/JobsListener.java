package scratch.aftershockStatistics.cmu;

import javax.jms.*;
import javax.jms.MessageListener;

import com.fasterxml.jackson.core.JsonParseException;
import com.fasterxml.jackson.databind.JsonMappingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.mongodb.morphia.Datastore;
import org.mongodb.morphia.query.Query;
import org.opensha.commons.data.siteData.impl.TectonicRegime;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.aftershockStatistics.*;

import java.io.IOException;
import java.util.GregorianCalendar;

/**
 * Created by clark on 11/11/2016.
 */
public class JobsListener implements MessageListener {
    private String consumerName;

    public JobsListener(String consumerName) {
        this.consumerName = consumerName;
    }

    public void onMessage(Message message) {
        try {
            String objectMessage;
            if (message instanceof TextMessage) {
                TextMessage textMessage = (TextMessage) message;

                objectMessage = textMessage.getText();
            }else {
                BytesMessage bytesMessage = (BytesMessage) message;
                byte data[] = new byte[(int) bytesMessage.getBodyLength()];
                bytesMessage.readBytes(data);

                objectMessage = new String(data);
            }

            System.out.println(consumerName + " received " + objectMessage);

            ObjectMapper mapper = new ObjectMapper();
            JobMessage job = mapper.readValue(objectMessage, JobMessage.class);


            /* Get Data from Shadow Store */
            Datastore datastore = MongoDBHelper.INSTANCE.getDatastore();

            /* Run Job */

            // this is the date range for which we fetching aftershock data
            double dataMinDays = job.dataMinDays;
            double dataMaxDays = job.dataMaxDays;



            /*// this is the date range for which we are forecasting
            double forecastMinDays = 0;
            double forecastMaxDays = 7;

            // depth range. event web service (or at least java wrapper) doesn't allow infinite depth, so set arbitrarily high
            double minDepth = 0; //job.minDepth
            double maxDepth = 1000; //job.maxDepth*/



            Query<ObsEqkSequence> query = datastore.createQuery(ObsEqkSequence.class);
            query.field("id").equal(job.id);
            ObsEqkSequence shockSeq = query.get();
            if(shockSeq == null) {
                System.out.println("No AfterShock data found for ID: " + job.id);
                return;
            }

            ModelResults_RJ model = new ModelResults_RJ();

            /*
             * Fetch generic aftershock parameters
             */
            GenericRJ_ParametersFetch genericFetch = new GenericRJ_ParametersFetch();
            TectonicRegime regime = genericFetch.getRegion(shockSeq.getMainShock().getHypocenterLocation());
            GenericRJ_Parameters genericParams = genericFetch.get(regime);

            model.setRegime(regime);
            model.setGenericParams(genericParams);

            // Make generic model
            ModelResult_RJGeneric genericModel = new ModelResult_RJGeneric(shockSeq.getMainShock().getMag(), genericParams);
            genericModel.maxLikelihood_a = genericModel.getMaxLikelihood_a();
            genericModel.maxLikelihood_b = genericModel.getMaxLikelihood_c();
            genericModel.maxLikelihood_c = genericModel.getMaxLikelihood_p();

            model.setGenericModel(genericModel);

            /*
             * Calculate sequence specific model
             */
            //double g = 0.25;
            double g = job.g;
            //double h = 1.0;
            double h = job.h;
            //double mCat = 4.5;
            double mCat = job.magCat;
            //double b = genericParams.get_bValue();
            double b = job.b;
            if (b < 0){
                b = genericParams.get_bValue();
            }
            double p = genericParams.get_pValue();
            double c = genericParams.get_cValue();
            double aMin = job.minA;
            double aMax = job.maxA;
            int aNum = 101;
            ModelResult_RJSpecific seqSpecificModel = new ModelResult_RJSpecific(
                    shockSeq.getMainShock(), shockSeq.getAfterShocks(), mCat, g, h, b,
                    dataMinDays, dataMaxDays,
                    aMin, aMax, aNum, p, p, 1, c, c, 1);
            seqSpecificModel.maxLikelihood_a = seqSpecificModel.getMaxLikelihood_a();
            seqSpecificModel.maxLikelihood_b = seqSpecificModel.getMaxLikelihood_c();
            seqSpecificModel.maxLikelihood_c = seqSpecificModel.getMaxLikelihood_p();

            model.setSpecificModel(seqSpecificModel);

            /*
             * Calculate Bayesian combination model
             */
            ModelResult_RJBayesian bayesianModel = new ModelResult_RJBayesian(seqSpecificModel, genericModel);
            bayesianModel.maxLikelihood_a = bayesianModel.getMaxLikelihood_a();
            bayesianModel.maxLikelihood_b = bayesianModel.getMaxLikelihood_c();
            bayesianModel.maxLikelihood_c = bayesianModel.getMaxLikelihood_p();

            model.setBayesianModel(bayesianModel);

            /**
             *  Forecasting
             */

            GregorianCalendar eventDate = shockSeq.getMainShock().getOriginTimeCal();
            GregorianCalendar startDate = new GregorianCalendar();

            Double minDays = dataMaxDays;

            double startTime = eventDate.getTime().getTime() + minDays * ProbabilityModelsCalc.MILLISEC_PER_DAY;
            startDate.setTimeInMillis((long)startTime);

            USGS_AftershockForecast forecast_generic = new USGS_AftershockForecast(genericModel, eventDate, startDate);
            USGS_AftershockForecast forecast_specific = new USGS_AftershockForecast(seqSpecificModel, eventDate, startDate);
            USGS_AftershockForecast forecast_bayesian = new USGS_AftershockForecast(bayesianModel, eventDate, startDate);

            model.setForecast_generic(forecast_generic.buildJSON());
            model.setForecast_specific(forecast_specific.buildJSON());
            model.setForecast_bayesian(forecast_bayesian.buildJSON());

            model.eventDataId = job.id;
            model.eventId = shockSeq.getMainShock().getEventId();
            model.setSampleTime(startDate.getTimeInMillis());

            datastore.save(model);

        } catch (JMSException e) {
            e.printStackTrace();
        } catch (JsonParseException e){
            e.printStackTrace();
        } catch (JsonMappingException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e){
            e.printStackTrace();
        } catch (Throwable e){
            e.printStackTrace();
        }
    }
}
