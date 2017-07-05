package scratch.aftershockStatistics;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.siteData.AbstractSiteData;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.data.siteData.impl.STREC_DataWrapper;
import org.opensha.commons.data.siteData.impl.TectonicRegime;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.ExceptionUtils;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class GenericRJ_ParametersFetch {
	
	private OrderedSiteDataProviderList regimeProviders;
	
	private Map<TectonicRegime, GenericRJ_Parameters> dataMap;
	
	public GenericRJ_ParametersFetch() {
		// from Morgan via e-mail 10/6/2015
		dataMap = Maps.newHashMap();
		
		// array contents: [ a, p, c ];
		
		URL paramsURL = GenericRJ_ParametersFetch.class.getResource("PageEtAlGenericParams_032116.csv");
		try {
			CSVFile<String> csv = CSVFile.readURL(paramsURL, true);
			Preconditions.checkState(csv.getNumCols() == 9, "unexpected number of columns: %s", csv.getNumCols());
			
			for (int row=1; row<csv.getNumRows(); row++) {
				String regimeName = csv.get(row, 0).trim();
				TectonicRegime regime = TectonicRegime.forName(regimeName);
				Preconditions.checkNotNull(regime, "Unknown regime: %s", regimeName);
				
				// columns: regionName,aValue_MaxLike,pValue_MaxLike,aValue_Mean,aValue_Sigma,
				//						aValue_Sigma1,aValue_Sigma0,bValue,cValue
				
				double pValue = Double.parseDouble(csv.get(row, 2));
				double aValue_mean = Double.parseDouble(csv.get(row, 3));
				double aValue_sigma = Double.parseDouble(csv.get(row, 4));
				double aValue_sigma1 = Double.parseDouble(csv.get(row, 5));
				double aValue_sigma0 = Double.parseDouble(csv.get(row, 6));
				double bValue = Double.parseDouble(csv.get(row, 7));
				double cValue = Double.parseDouble(csv.get(row, 8));
				
				dataMap.put(regime, new GenericRJ_Parameters(
						aValue_mean, aValue_sigma, aValue_sigma0, aValue_sigma1, bValue, pValue, cValue));
			}
		} catch (IOException e) {
			ExceptionUtils.throwAsRuntimeException(e);
		}
		
		ArrayList<SiteData<?>> providers = Lists.newArrayList();
		
		// now hardcoded ones
		// via e-mail from Jeanne 2/8/17 and 2/9/17, subject "OAF To-Do List; Jan 31, 2017":
		/*
		 * I dug up the relevant email thread on the California regional parameters.  The summary is to use the modified
		 * R&J parameters a=−1.85, b = 0.91, p = 1.08, and c = 0.05.  R&J don't define a-sigma, but I think we were okay
		 * with using the global values of σ0=0.49 and σ1=750.  Use the same equation for the completeness, except with
		 * Mcat= 2.5.  The California spatial region is given in the attached files that Andy provided (the union of
		 * region.ncsn and region.scsn).
		 */
		// aValue_sigma value of 1.76 is from the second e-mail, "Integrating equation 8 over M5 to M8, weighted by a MFD
		// with b=1, the sigma value is 1.76. I know this is kind of large, but this is the extrapolation to lower
		// magnitude of the observed relation between magnitude and sigma."
		providers.add(new HardcodedRegimeFetch(TectonicRegime.CALIFORNIA, buildANSS_CA_Region()));
		dataMap.put(TectonicRegime.CALIFORNIA, new GenericRJ_Parameters(-1.85, 1.76, 0.49, 750, 0.91, 1.08, 0.05));
		
		// add in STREC as the fallback
		providers.add(new STREC_DataWrapper());
		regimeProviders = new OrderedSiteDataProviderList(providers);
	}
	
	private static Region buildANSS_CA_Region() {
		LocationList locs = new LocationList();
		// region.scsn locations
		locs.add(new Location(36.6847,   -117.793));
		locs.add(new Location(35.8000,   -116.400));
		locs.add(new Location(34.0815,   -114.472));
		locs.add(new Location(32.0000,   -114.333));
		locs.add(new Location(32.0000,   -120.500));
		locs.add(new Location(34.6945,   -121.380));
		// locs.add(new Location(36.6847,   -117.793)); // this was to close it, instead skip and continue on to NCSN
		// region.ncsn locations, in reverse to match the ordering and continue on in SCSN
		//locs.add(new Location(34.6945,   -121.380)); // ignore, duplicate with  SCSN
		locs.add(new Location(40.0000,   -125.500));
		locs.add(new Location(43.0200,   -125.000));
		locs.add(new Location(42.0000,   -122.700));
		locs.add(new Location(42.0000,   -121.417));
		locs.add(new Location(39.5000,   -120.750));
		locs.add(new Location(37.7500,   -119.500));
		locs.add(new Location(37.7500,   -118.250));
		locs.add(new Location(36.6847,   -117.793)); // this is also the first one from SCSN, and closes the region
		//locs.add(new Location(34.6945,   -121.380)); // this was to close it, skip it (and we're already closed with SCSN)
		
		return new Region(locs, null);
	}
	
	/**
	 * Fetches a Garcia region for the given location and returns default omori parameters
	 * @param region
	 * @return array of parameters: {a, p, c};
	 */
	public GenericRJ_Parameters get(Location loc) {
		TectonicRegime region = getRegion(loc);
		return get(region);
	}
	
	/**
	 * 
	 * @param region
	 * @return array of parameters: {a, p, c};
	 */
	public GenericRJ_Parameters get(TectonicRegime region) {
		return dataMap.get(region);
	}
	
	public TectonicRegime getRegion(Location loc) {
		ArrayList<SiteDataValue<?>> datas = regimeProviders.getBestAvailableData(loc);
		if (datas.isEmpty())
			return null;
		return (TectonicRegime)datas.get(0).getValue();
	}
	
	private class HardcodedRegimeFetch extends AbstractSiteData<TectonicRegime> {
		private TectonicRegime regime;
		private Region region;
		
		public HardcodedRegimeFetch( TectonicRegime regime, Region region) {
			Preconditions.checkArgument(regime != null);
			Preconditions.checkArgument(region != null);
			this.region = region;
			this.regime = regime;
		}
		
		@Override
		public Region getApplicableRegion() {
			return region;
		}
		@Override
		public double getResolution() {
			return 0;
		}
		@Override
		public String getName() {
			return "Harcoded for Regime '"+regime.name()+"', Region '"+region.getName()+"'";
		}
		@Override
		public String getShortName() {
			return "Hardcoded";
		}
		@Override
		public String getDataType() {
			return TYPE_TECTONIC_REGIME;
		}
		@Override
		public String getDataMeasurementType() {
			return TYPE_FLAG_MEASURED;
		}
		@Override
		public Location getClosestDataLocation(Location loc) throws IOException {
			return loc;
		}
		@Override
		public TectonicRegime getValue(Location loc) throws IOException {
			if (region.contains(loc))
				return regime;
			return null;
		}
		@Override
		public boolean isValueValid(TectonicRegime el) {
			return regime.equals(el);
		}
		@Override
		public String getMetadata() {
			return getName();
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
		
		GenericRJ_ParametersFetch fetch = new GenericRJ_ParametersFetch();
		for (Location loc : locs) {
			TectonicRegime regime = fetch.getRegion(loc);
			System.out.println(loc+", "+regime+": "+fetch.get(regime));
		}
	}

}
