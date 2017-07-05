package scratch.aftershockStatistics.cmu;

import javax.jms.BytesMessage;
import javax.jms.Destination;
import javax.jms.JMSException;
import javax.jms.TextMessage;
import javax.jms.Message;
import javax.jms.MessageListener;
import javax.jms.MessageProducer;
import javax.jms.Session;


import com.fasterxml.jackson.core.JsonParseException;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.JsonMappingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.google.common.base.Preconditions;

import org.apache.activemq.ScheduledMessage;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import org.mongodb.morphia.Datastore;
import org.mongodb.morphia.query.Query;
import org.mongodb.morphia.query.UpdateOperations;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.geo.Region;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import scratch.aftershockStatistics.ComcatAccessor;
import org.mongodb.morphia.aggregation.Sort;

import java.io.IOException;
import java.net.UnknownHostException;
import java.time.Instant;
import java.util.Iterator;
import java.util.List;

/**
 *
 */
public class ComcatListener implements MessageListener {

	long dayToMilli = 86400000;
	long monthToMillis = 2629746000L;
	Session session = null;
	public ComcatListener(Session session) {
		this.session=session;
	}
	/**
	 * callback for messages from subscribed queues
	 * @param message
	 * @return void
	 */

	public void onMessage(Message message) {
		String comcatMessage;
		try {
			String dest = message.getJMSDestination().toString();
			comcatMessage = getMsgTxt(message);
			ObjectMapper mapper = new ObjectMapper();
			ComcatMessage comcat = mapper.readValue(comcatMessage, ComcatMessage.class);
			System.out.println(mapper.writeValueAsString(comcat));
			
			/*if event notified by PDL queue is already in DB ignore it*/
			
			if(dest.equals("queue://PDL") && eventExists(comcat.eventId)==true)
			{
				System.out.println("Event Exists");
				return;
			}
			
			/*fetch event, its aftershocks and queue for forecast generation*/
			processEvent(comcat);
		}
		catch (JMSException e) {
			e.printStackTrace();
		} catch (UnknownHostException e) {
			e.printStackTrace();
		} catch (JsonParseException e) {
			e.printStackTrace();
		} catch (JsonMappingException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch(Exception e){
			e.printStackTrace();
		}

	}

	/**
	 * process the message received to get the message in string format
	 * @param message
	 * @return String
	 */

	public String getMsgTxt(Message message) throws JMSException
	{
		/*JMS messages from comcat queue*/
		if(message instanceof TextMessage)  
		{
			TextMessage msg = (TextMessage) message;
			return msg.getText();
		}
		else
		{
			/*byte messages from PDL queue as sent using stomp */
			BytesMessage byteMessage = (BytesMessage) message; 
			byte[] byteArr = new byte[(int)byteMessage.getBodyLength()];
			byteMessage.readBytes(byteArr);
			return new String(byteArr);
		}
	}
	/**
	 * Determines whether an event already exists in the Aftershock database or is within an active region
	 * @param eventId
	 * @return boolean
	 */

	public boolean eventExists(String eventId) throws UnknownHostException
	{
		boolean result=false;
		Datastore datastore = MongoDBHelper.INSTANCE.getDatastore();
		Query<ObsEqkSequence> q= datastore.createQuery(ObsEqkSequence.class);
		/*counts the records for the eventId*/
		long count = q.field("eqkList.eventId").equal(eventId).countAll(); 
		/*if event exists return true*/
		if (count!=0) 
			return true;
		/*
		 * check whether event is present in Comcat db if not it need not be processed
		 */
		ComcatAccessor accessor = new ComcatAccessor();
		ObsEqkRupture eqk = accessor.fetchEvent(eventId);
		if(eqk==null)
		{
			System.out.println("Event does not exist in Comcat");
			return true;
		}
	    /* If event does not exist check if it comes under any actively processed regions*/
		q= datastore.createQuery(ObsEqkSequence.class);
	    /*obtain all regions with isActive flag set to true from database*/
		List<ObsEqkSequence> isActiveList = q.field("isActive").equal(true).retrievedFields(true, "mainShockRegion").asList();
		isActiveList.addAll(q.field("isActive").equal(true).retrievedFields(true, "foreShockRegion").asList());
		if(isActiveList!=null)
		{
			/*
			 * traverse through list of active regions to see if any of the regions contain the event 
			 */
			for(ObsEqkSequence eqkList:isActiveList)
			{
				RegionParams regionParams = eqkList.getMainShockRegion();
				Region region = regionParams.buildRegion();
				if(region.contains(eqk.getHypocenterLocation()))
				{
					result= true;
					break;
				}
				if (eqkList.getForeShockRegion()!=null)
				{
					region = regionParams.buildRegion();
					if(region.contains(eqk.getHypocenterLocation()))
					{
						result= true;
						break;
					}
				}
			}
		}

		return result;
	}
	/**
	 * Processes the event and enqueues it for forecasts
	 * @param message
	 * @return void
	 */
	public  void processEvent(ComcatMessage message) throws UnknownHostException, JMSException
	{
      /* Get Data from Shadow Store */
		String eventID=message.eventId;
		int probabilityEthqkFlag=1;
		double dataMinDays = 0; //days after event to start fetching aftershocks
		double dataMaxDays; //days after event until which aftershocks are to be fetched
		double minDepth = message.minLocation.depth;
		double maxDepth = message.maxLocation.depth; 
		ObsEqkRupList afs_final=null;
		Datastore datastore = MongoDBHelper.INSTANCE.getDatastore();
		/*change the active state for catalogs in datastore conditionally*/
		updateActiveStatus( datastore, message); 
        
		/*
		 * fetch the event from comcat
		 */
		ComcatAccessor accessor = new ComcatAccessor();
		ObsEqkRupture initialMainShock = accessor.fetchEvent(eventID);
		if(initialMainShock==null)
			System.out.println("Mainshock is null");
		Preconditions.checkNotNull(initialMainShock, "Error fetching mainshock '%s'", eventID);
		long now = Instant.now().toEpochMilli();
	    /*if dataMaxDays is 0 get all aftershocks until now*/
		if(message.dataMaxDays==0){
			dataMaxDays =  ((double)(now-initialMainShock.getOriginTime()) / (double)dayToMilli);
			System.out.println(dataMaxDays);
			message.dataMaxDays=dataMaxDays;
		}
		else
		{
			dataMaxDays=message.dataMaxDays;
		}
		dataMinDays=message.dataMinDays;
		ObsEqkSequence shockSeq = new ObsEqkSequence(initialMainShock.getEventId());
		shockSeq.getEqkList().add(initialMainShock);
		shockSeq.setSampleTime(now);  //time of fetching data from comcat
		/**
		 *do not predict for more than a year from today otherwise run forescasts for at least a month
		 */
		if(now-dataMaxDays*dayToMilli<0)
			probabilityEthqkFlag=0;
		else if(dataMaxDays*dayToMilli > monthToMillis)
			probabilityEthqkFlag=getProbabilityFlag(message.eventId); //checks the previous forecast to see if the event is worth processing
		WC1994_MagLengthRelationship wcMagLen = new WC1994_MagLengthRelationship();
		RegionParams region = new RegionParams(message);
		if(message.regionType.equals(RegionParams.WCCIRCULAR))
		{
			region.setRadius(wcMagLen.getMedianLength(shockSeq.getMainShock().getMag()));
		}
		if(region.isCircular()==true && message.centerType.equals(RegionParams.CUSTOM)==false)
		{
			region.setCicleCenter(initialMainShock.getHypocenterLocation());
		}

		Region finalRegion=region.buildRegion();
	    /*Loop to get the greatest aftershock larger than mainshock*/
		while(true)
		{
			if(message.regionType.equals(RegionParams.WCCIRCULAR))
			{
				region.setRadius(wcMagLen.getMedianLength(shockSeq.getMainShock().getMag()));
			}
			if(region.isCircular()==true && message.centerType.equals(RegionParams.CENTROID))
				finalRegion=region.buildRegion(shockSeq.getMainShock(),dataMinDays, dataMaxDays);
			shockSeq.setMainShockRegion(region);
			afs_final= accessor.fetchAftershocks(shockSeq.getMainShock(), dataMinDays, dataMaxDays, minDepth, maxDepth, finalRegion);
			if(afs_final.isEmpty())
			{
				break;
			}
			shockSeq.getEqkList().addAll(afs_final);
			afs_final.sortByMag();
			ObsEqkRupture largestAfs=afs_final.get(afs_final.size()-1);
			if(largestAfs.getMag() > shockSeq.getMainShock().getMag())
			{
				shockSeq.setMainShockId(largestAfs.getEventId());
			}
			else
			{
				break;
			}
		}
		if(!shockSeq.getMainShock().getEventId().equals(initialMainShock.getEventId())){
			shockSeq.setForeShockId(initialMainShock.getEventId());
			RegionParams foreShock=new RegionParams(message);
			shockSeq.setForeShockRegion(foreShock);
		}
		if(probabilityEthqkFlag==1 && message.persistent==true && message.emulate==false){
			shockSeq.setisActive(true);
		}
		datastore.save(shockSeq);
	     /*schedule the event for forecasting now and after a certain interval*/
		try{
			ObjectMapper mapper = new ObjectMapper();
			
			JobMessage job = new JobMessage(shockSeq.getId(),message.dataMaxDays,message.dataMinDays,message.minA,message.maxA,message.minP,message.maxP,message.minC, message.maxC, message.b, message.g, message.h, message.magCat);
			String jobString = mapper.writeValueAsString(job);

			schedEvent("Jobs.RJ",jobString,0);
			/*schedule event for later processing only if the mag 3 probablity threshold crosses 10%*/
			if(probabilityEthqkFlag==1 && (message.persistent==true || message.emulate==true) ){
				long diff = now-shockSeq.getMainShock().getOriginTime();
				long delay=getDelay(diff);
				/*
				 * if emulating an old event dataMaxDays and delay have to be set acoordingly 
				 */
				if(message.emulate==true)
				{
					double interval = getDelay((long)(dataMaxDays*dayToMilli));
					interval=interval/dayToMilli;
					interval= interval+dataMaxDays;
					message.dataMaxDays = interval; /*cumulative maxDays for emulation*/
					delay=15000; //queue for processing after 15 seconds to lessen race conditions b/w jobs and Comcat queues
				}
				String comcatString = mapper.writeValueAsString(message);
				schedEvent("Comcat",comcatString,delay);

			}
			System.out.println("Processing ended");

		} catch (JsonProcessingException e) {
			e.printStackTrace();
		}

	}

	/**
	 * Gets the probability of eqk with mag 3 happening
	 * @param eventId of catalog on which forecast is based
	 * @return int
	 */
	public int getProbabilityFlag(String eventId)
	{
		int a=0;
		System.out.println("enters probability checking");
		System.out.println(eventId);
		if(eventId!=null)
		{	System.out.println("is not null");
			Datastore datastore = MongoDBHelper.INSTANCE.getDatastore();
			/*
			 *get latest forecast for the mentioned eventid 
			 */
			Query<ModelResults_RJ> q = datastore.createQuery(ModelResults_RJ.class).field("eventId").equal(eventId);
			Iterator<ModelResults_RJ> itr= datastore.createAggregation(ModelResults_RJ.class).
													match(q).
													sort(Sort.descending("sampleTime")).
													limit(1).aggregate(ModelResults_RJ.class);


			while(itr.hasNext())
			{
				ModelResults_RJ forecast = itr.next();
				String  array = forecast.getForecast_bayesian().toJSONString();
				try {
					JSONObject jsonObj = new JSONObject(array);
					JSONArray forecastArray = jsonObj.getJSONArray("forecast");
					JSONObject obj = forecastArray.getJSONObject(0);
					JSONArray binsArray=obj.getJSONArray("bins");
					JSONObject f=binsArray.getJSONObject(0);
					double probability=f.getDouble("probability");
					System.out.println("Probability = "+(probability*100));
					if(probability*100>10){
						a=1;
						System.out.println("returning "+ a);
						return 1;
					}
					else{
						System.out.println("returning "+ a);
						return 0;
					}
				} catch (JSONException e) {
					e.printStackTrace();
				}
			}
		}
		return 1;
	}

	/**
	 * Updates the active flag for events if they are persistent and not emulated to false
	 * @param datastore and message
	 * @return void
	 */
	public void updateActiveStatus(Datastore datastore, ComcatMessage message)
	{
		if(message.persistent==true && message.emulate==false){
			Query<ObsEqkSequence> q = datastore.createQuery(ObsEqkSequence.class);
			q.and(
					q.criteria("eqkList.eventId").equal(message.eventId),
					q.criteria("mainShockRegion.regionType").equal(message.regionType));
			UpdateOperations<ObsEqkSequence> update = datastore.createUpdateOperations(ObsEqkSequence.class).set("isActive",false);
			datastore.update(q,update);
		}
	}

	/**
	 * Schedules a message in a queue with the given delay
	 * @param queuename, msg, delay
	 * @return void
	 */
	public void schedEvent(String queue, String str, long delay) throws JMSException
	{
		final Destination dest = this.session.createQueue(queue);
		MessageProducer producer = this.session.createProducer(dest);
		Message msg = this.session.createTextMessage(str);
		msg.setLongProperty(ScheduledMessage.AMQ_SCHEDULED_DELAY, delay);
		producer.send(msg);
	}

	/**
	 * gets delay of schedule in queue based on diff between event and processing time
	 * @param long
	 * @return long 
	 */
	public long getDelay(long diff)
	{
		long delay;
		if(diff >= 31556952000L) // a year
			delay = 31556952000L;
		else if(diff >= 2629746000L) // a month
			delay = 2629746000L;
		else if(diff >=604800000) // a week
			delay=604800000;
		else if(diff>=86400000) // a day
			delay=86400000;
		else
			delay=14400000;   //4 hrs
		return delay;
	}
}