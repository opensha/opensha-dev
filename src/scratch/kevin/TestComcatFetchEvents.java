package scratch.kevin;

import java.math.BigDecimal;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.List;

import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;

import com.google.common.collect.Lists;

import gov.usgs.earthquake.event.EventQuery;
import gov.usgs.earthquake.event.EventWebService;
import gov.usgs.earthquake.event.Format;
import gov.usgs.earthquake.event.JsonEvent;

public class TestComcatFetchEvents {
	
	public static void main(String[] args) throws Exception {
		EventWebService service = new EventWebService(new URL("http://earthquake.usgs.gov/fdsnws/event/1/"));
		EventQuery query = new EventQuery();
//		query.setStartTime(new GregorianCalendar(2015, 0, 1).getTime());
//		query.setEndTime(new Date());
//		query.setMinMagnitude(new BigDecimal(2.5d));
		Region reg = new CaliforniaRegions.RELM_TESTING();
//		query.setMinLatitude(new BigDecimal(reg.getMinLat()));
//		query.setMaxLatitude(new BigDecimal(reg.getMaxLat()));
//		query.setMinLongitude(new BigDecimal(reg.getMinLon()));
//		query.setMaxLongitude(new BigDecimal(reg.getMaxLon()));
		query.setEventId("ci37166079");
		System.out.println(service.getUrl(query, Format.GEOJSON));
		List<JsonEvent> events = service.getEvents(query);
		
		System.out.println("Fetched "+events.size()+" events");
		
		ObsEqkRupList list = new ObsEqkRupList();
		
		for (JsonEvent event : events) {
			Location hypoLoc = new Location(event.getLatitude().doubleValue(),
					event.getLongitude().doubleValue(), event.getDepth().doubleValue());
			if (!reg.contains(hypoLoc))
				continue;
			ObsEqkRupture rup = new ObsEqkRupture(event.getEventId().toString(), event.getTime().getTime(),
					hypoLoc, event.getMag().doubleValue());
			list.add(rup);
			System.out.println(rup.getEventId()+" "+new Date(rup.getOriginTime())+" "+hypoLoc+" "+rup.getMag());
		}
		
		System.out.println(list.size()+" of which were inside region");
		
		
	}

}
