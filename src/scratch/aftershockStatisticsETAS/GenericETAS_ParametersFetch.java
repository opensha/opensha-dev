package scratch.aftershockStatisticsETAS;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.data.siteData.impl.STREC_DataWrapper;
import org.opensha.commons.data.siteData.impl.TectonicRegime;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.util.ExceptionUtils;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class GenericETAS_ParametersFetch {
	
	private final static boolean D = false; //debug
	
	private OrderedSiteDataProviderList regimeProviders;
	
	private Map<TectonicRegime, GenericETAS_Parameters> dataMap;
	
	public GenericETAS_ParametersFetch() {
		ArrayList<SiteData<?>> providers = Lists.newArrayList();
		providers.add(new STREC_DataWrapper());
		regimeProviders = new OrderedSiteDataProviderList(providers);
		
		dataMap = Maps.newHashMap();

		URL paramsURL = GenericETAS_ParametersFetch.class.getResource("resources/vdEGenericETASParams_080518.csv"); //updated CALIFORNIA added Mref, Mmax, columns 
		
		if(D) System.out.println(paramsURL);
		
		try {
			CSVFile<String> csv = CSVFile.readURL(paramsURL, true);
			Preconditions.checkState(csv.getNumCols() == 16, "unexpected number of columns: %s", csv.getNumCols());
			
			for (int row=1; row<csv.getNumRows(); row++) {
				String regimeName = csv.get(row, 0).trim();
				TectonicRegime regime = TectonicRegime.forName(regimeName);
				
				
				// columns: regionName,aValue_MaxLike,pValue_MaxLike,aValue_Mean,aValue_Sigma,
				//						aValue_Sigma1,aValue_Sigma0,bValue,cValue
				
				double aValue_mean = Double.parseDouble(csv.get(row, 1));
				double aValue_sigma = Double.parseDouble(csv.get(row, 2));
//				double log_cValue = Double.parseDouble(csv.get(row, 3));			
				double pValue = Double.parseDouble(csv.get(row, 4));
				double covaa = Double.parseDouble(csv.get(row, 5));
				double covpp = Double.parseDouble(csv.get(row, 6));
				double covcc = Double.parseDouble(csv.get(row, 7));
				double covap = Double.parseDouble(csv.get(row, 8));
				double covac = Double.parseDouble(csv.get(row, 9));
				double covcp = Double.parseDouble(csv.get(row, 10));
				int numSequences = Integer.parseInt(csv.get(row, 11));
				double alpha = Double.parseDouble(csv.get(row, 12));
				double bValue = Double.parseDouble(csv.get(row, 13));
				double refMag = Double.parseDouble(csv.get(row, 14));
				double maxMag = Double.parseDouble(csv.get(row, 15));

//				double cValue = Math.pow(10, log_cValue);
				// due to memory concerns, we're replacing the generic c-value with the global average c-value
				double cValue = Math.pow(10, -2.5);
				
				double pValue_sigma = Math.sqrt(covpp*numSequences);	//this needs to be replaced with a real estimate
				// due to memory concerns, we're replacing logcValue_sigma with a value that makes 3*sigma value equal to -5, 0
//				double logcValue_sigma = Math.sqrt(covcc*numSequences);		//this needs to be replaced with a real estimate
				double logcValue_sigma = 0.8333;		//this needs to be replaced with a real estimate
				
				double bValue_sigma = 0.1;//this needs to be replaced with a real estimate
//				double maxMag = 9.5;//this needs to be replaced with a real estimate
				
				dataMap.put(regime, new GenericETAS_Parameters(
						aValue_mean, aValue_sigma, pValue, pValue_sigma, cValue, logcValue_sigma, covaa, covpp, covcc, covap, covac, covcp, numSequences, alpha,  bValue, bValue_sigma, refMag, maxMag));
				
			}
		} catch (IOException e) {
			ExceptionUtils.throwAsRuntimeException(e);
		}
	}
	
	/**
	 * Fetches a Garcia region for the given location and returns default omori parameters
	 * @param region
	 * @return array of parameters: {a, p, c};
	 */
	public GenericETAS_Parameters get(Location loc) {
		TectonicRegime region = getRegion(loc);
		return get(region);
	}
	
	/**
	 * 
	 * @param region
	 * @return array of parameters: {a, p, c};
	 */
	public GenericETAS_Parameters get(TectonicRegime region) {
		return dataMap.get(region);
	}
	
	public TectonicRegime getRegion(Location loc) {
		try {
			ArrayList<SiteDataValue<?>> datas = regimeProviders.getBestAvailableData(loc);
			if (datas.isEmpty())
				return null;
			return (TectonicRegime)datas.get(0).getValue();
			
		} catch (Exception e){
			System.out.println("Couldn't retrieve tectonic regime for location " + loc + ". Using global average.");
			return TectonicRegime.GLOBAL_AVERAGE;
		}
	}
	
	public static void main(String[] args) {
		LocationList locs = new LocationList();
		
		locs.add(new Location(28.2305, 84.7314, 8.22));
		locs.add(new Location(35, -118, 7d));
		locs.add(new Location(35, -50, 7d));
		for (int i=0; i<5; i++) {
			double lat = 180d*Math.random()-90d;
			double lon = 360d*Math.random()-180d;
			double depth = 20d*Math.random();
			locs.add(new Location(lat, lon, depth));
		}
		
		GenericETAS_ParametersFetch fetch = new GenericETAS_ParametersFetch();
		for (Location loc : locs) {
			TectonicRegime regime = fetch.getRegion(loc);
			System.out.println(loc+", "+regime+": "+fetch.get(regime));
		}
	}

}
