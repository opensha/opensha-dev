package scratch.aftershockStatistics;

import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import scratch.aftershockStatistics.OAFTectonicRegime;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;

import com.google.common.base.Preconditions;

public class HardCodedExample {

	public static void main(String[] args) {
		/*
		 * Inputs
		 */
		// main shock event ID
		String eventID = "us20002926";
		
		// this is the date range for which we fetching aftershock data
		double dataMinDays = 0;
		double dataMaxDays = 7;
		
		// this is the date range for which we are forecasting
		double forecastMinDays = 0;
		double forecastMaxDays = 7;
		
		// depth range. event web service (or at least java wrapper) doesn't allow infinite depth, so set arbitrarily high
		double minDepth = 0;
		double maxDepth = 1000;
		
		/*
		 * Fetch mainshock
		 */
		// this is our interface to Comcat
		ComcatAccessor accessor = new ComcatAccessor();
		
		ObsEqkRupture mainshock = accessor.fetchEvent(eventID);
		Preconditions.checkNotNull(mainshock, "Error fetching mainshock '%s'", eventID);
		
		/* 
		 * Determine aftershock region and fetch aftershocks
		 */
		// determine search radius
		WC1994_MagLengthRelationship wcMagLen = new WC1994_MagLengthRelationship();
		double radius = wcMagLen.getMedianLength(mainshock.getMag());
		
		// circular region centered around main shock location
		Region tempRegion = new Region(mainshock.getHypocenterLocation(), radius);
		
		// fetch aftershocks
		ObsEqkRupList aftershocks = accessor.fetchAftershocks(mainshock, dataMinDays, dataMaxDays, minDepth, maxDepth, tempRegion);
		
		// now find centroid and use that to rebuild the aftershock list
		if (aftershocks.isEmpty()) {
			System.out.println("No aftershocks found, skipping centroid");
		} else {
			Location centroid = AftershockStatsCalc.getCentroid(mainshock, aftershocks);
			Region finalRegion = new Region(centroid, radius);

			aftershocks = accessor.fetchAftershocks(mainshock, dataMinDays, dataMaxDays, minDepth, maxDepth, finalRegion);
		}
		
		/*
		 * Fetch generic aftershock parameters
		 */
		GenericRJ_ParametersFetch genericFetch = new GenericRJ_ParametersFetch();
		OAFTectonicRegime regime = genericFetch.getRegion(mainshock.getHypocenterLocation());
		GenericRJ_Parameters genericParams = genericFetch.get(regime);
		System.out.println("Generic params for "+regime+": "+genericParams);
		
		// Make generic model
		RJ_AftershockModel_Generic genericModel = new RJ_AftershockModel_Generic(mainshock.getMag(), genericParams);
		
		/*
		 * Calculate sequence specific model
		 */
		double g = 0.25;
		double h = 1.0;
		double mCat = 4.5;
		double b = genericParams.get_bValue();
		double p = genericParams.get_pValue();
		double c = genericParams.get_cValue();
		double aMin = -4.5;
		double aMax = -0.5;
		int aNum = 101;
		RJ_AftershockModel_SequenceSpecific seqSpecificModel = new RJ_AftershockModel_SequenceSpecific(
				mainshock, aftershocks, mCat, g, h, b,
				dataMinDays, dataMaxDays,
				aMin, aMax, aNum, p, p, 1, c, c, 1);
		
		/*
		 * Calculate Bayesian combination model
		 */
		RJ_AftershockModel_Bayesian bayesianModel = new RJ_AftershockModel_Bayesian(seqSpecificModel, genericModel);
		
		/*
		 * Print results
		 */
		double minMag = 3d;
		double maxMag = 9d;
		double deltaMag = 0.1;
		int numMag = (int)((maxMag - minMag)/deltaMag + 1.5);
		
		EvenlyDiscretizedFunc mode = genericModel.getModalCumNumMFD(minMag, maxMag, numMag, forecastMinDays, forecastMaxDays);
		mode.setName("Generic Forecast Mode");
		System.out.println(mode);
		
		mode = seqSpecificModel.getModalCumNumMFD(minMag, maxMag, numMag, forecastMinDays, forecastMaxDays);
		mode.setName("Sequence Specific Forecast Mode");
		System.out.println(mode);
		
		mode = bayesianModel.getModalCumNumMFD(minMag, maxMag, numMag, forecastMinDays, forecastMaxDays);
		mode.setName("Bayesian Forecast Mode");
		System.out.println(mode);
	}

}
