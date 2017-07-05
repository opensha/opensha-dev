package scratch.aftershockStatisticsETAS;

import java.io.*;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.siteData.impl.TectonicRegime;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;

import com.google.common.base.Preconditions;

public class HardCodedETASExample {

	public static void main(String[] args) {
		/*
		 * Inputs
		 */
		// main shock event ID
		String eventID = "us20002926";
		boolean loadFromDisk = false; 
		String tempSaveFilename = "/Users/nvanderelst/Documents/LastRetrievedCatalog.sav";
		
		// this is the date range for which we fetching aftershock data
		double dataMinDays = 0;
		double dataMaxDays = 7;
		
		// this is the date range for which we are forecasting
		double forecastMinDays = 0;
		double forecastMaxDays = 31;
		int nSims = 1000;
		
		// depth range. event web service (or at least java wrapper) doesn't allow infinite depth, so set arbitrarily high
		double minDepth = 0;
		double maxDepth = 1000;
		
		/*
		 * Fetch mainshock
		 */
		ObsEqkRupture mainshock = new ObsEqkRupture();
		ObsEqkRupList aftershocks = new ObsEqkRupList();
		GenericETAS_Parameters genericParams = new GenericETAS_Parameters();
//		GenericRJ_Parameters genericParamsRJ = new GenericRJ_Parameters();
		
		if(loadFromDisk){
			// if the flag is set to load from disk, try to load 
			// the most recently downloaded catalog and fetched parameters
			try{
				FileInputStream f_in = new FileInputStream(tempSaveFilename);
				ObjectInputStream obj_in = new ObjectInputStream (f_in);

				Object obj = obj_in.readObject();
				// Cast object to a ObsEqkRupture
				mainshock = (ObsEqkRupture) obj;

				obj = obj_in.readObject();
				// Cast object to a ObsEqkRupList
				aftershocks = (ObsEqkRupList) obj;

				obj = obj_in.readObject();
				// Cast object to a GenericETAS_Parameters
				genericParams = (GenericETAS_Parameters) obj;
				
				obj = obj_in.readObject();
				// Cast object to a GenericRJ_Parameters
//				genericParamsRJ = (GenericRJ_Parameters) obj;

				obj_in.close();
				
			} catch (IOException e) {
				e.printStackTrace();
				System.out.println("Problem loading last catalog from disk. Trying Comcat instead.");
				loadFromDisk = false;
			} catch (ClassNotFoundException e) {
				e.printStackTrace();
				System.out.println("Problem loading last catalog from disk. Trying Comcat instead.");
				loadFromDisk = false;
			}
		}
		
		if(!loadFromDisk){
			//if the flag is not set to load from disk, or if that failed, go get it from Comcat.
			
			// this is our interface to Comcat
			ComcatAccessor accessor = new ComcatAccessor();
			mainshock = accessor.fetchEvent(eventID);
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
			aftershocks = accessor.fetchAftershocks(mainshock, dataMinDays, dataMaxDays, minDepth, maxDepth, tempRegion);

			// now find centroid and use that to rebuild the aftershock list
			if (aftershocks.isEmpty()) {
				System.out.println("No aftershocks found, skipping centroid");
			} else {
				Location centroid = ETAS_StatsCalc.getCentroid(mainshock, aftershocks);
				Region finalRegion = new Region(centroid, radius);
				aftershocks = accessor.fetchAftershocks(mainshock, dataMinDays, dataMaxDays, minDepth, maxDepth, finalRegion);
			}

			/*
			 * Fetch generic aftershock parameters
			 */
			GenericETAS_ParametersFetch genericFetch = new GenericETAS_ParametersFetch();
//			GenericRJ_ParametersFetch genericFetchRJ = new GenericRJ_ParametersFetch();
			TectonicRegime regime = genericFetch.getRegion(mainshock.getHypocenterLocation());
			Preconditions.checkNotNull(regime, "Error fetching tectonic regime");
			genericParams = genericFetch.get(regime);
//			genericParamsRJ = genericFetchRJ.get(regime);
			
			System.out.println("Generic params for "+regime+": "+genericParams);
		
			// write the retrieved data to a file in case the Internet goes down.
			try {
				FileOutputStream saveFile = new FileOutputStream(tempSaveFilename);
				ObjectOutputStream save = new ObjectOutputStream(saveFile);
				save.writeObject(mainshock);
				save.writeObject(aftershocks);
				save.writeObject(genericParams);
//				save.writeObject(genericParamsRJ);
				save.close();
				
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		System.out.println(genericParams);
		double refMag = 4.5;
		
		// generate forecast object using generic parameters (catalogs produced inside model)
		ETAS_AftershockModel_Generic genericETASmodel = new ETAS_AftershockModel_Generic(
				mainshock, aftershocks, genericParams,
				dataMinDays, dataMaxDays, forecastMinDays, forecastMaxDays, refMag, 9.5, 100, nSims); //maxMag = 9.5, maxGeneratons = 100; nsims = 0);

		// make one-month forecast
		ArbDiscrEmpiricalDistFunc M5distFunc = genericETASmodel.computeNum_DistributionFunc(dataMaxDays+forecastMinDays, 
				dataMaxDays+31, refMag);
		
		System.out.println("dataEndTime = " + genericETASmodel.dataEndTimeDays);
		
		
		EvenlyDiscretizedFunc cumNumWithTimeGen = genericETASmodel.getModalCumNumEventsWithTime( 5.0, dataMaxDays+forecastMinDays,  dataMaxDays+forecastMaxDays,  1);
		System.out.println(cumNumWithTimeGen);
		
		
		/*
		 * Calculate sequence specific model
		 */
		double g = 0.25;
		double h = 1.0;
		double mCat = 4.5;
		double b = genericParams.get_bValue();
		double alpha = 1.0;
//		double refMag = 4.5;
		double aMin = -3.56;
		double aMax = -1.56;
		double pMin = 1.08;
		double pMax = 1.08;
		double cMin = 0.0087096; 	// c as exponent
		double cMax = 0.0087096;	// log10 of c (likelihood is distributed log normal)
		int aNum = 41;
		int pNum = 1;
		int cNum = 1;
		
		// pass as vectors
		double[] aVec = ETAS_StatsCalc.linspace(aMin,aMax,aNum);
		double[] pVec = ETAS_StatsCalc.linspace(pMin,pMax,pNum);
		double[] cVec = ETAS_StatsCalc.logspace(cMin,cMax,cNum);
				
		System.out.println(mainshock.getMag());
//		Make sequence specific model -- fits and computes.
		ETAS_AftershockModel_SequenceSpecific seqSpecificETASmodel = new ETAS_AftershockModel_SequenceSpecific(
				mainshock, aftershocks,
				mCat, aVec, pVec, cVec, alpha, b, refMag, 
				dataMinDays, dataMaxDays, forecastMinDays, forecastMaxDays, 9.5, 100, nSims); //maxMag = 9.5, maxGeneratons = 100;););
		
		// make one-month forecast
		ArbDiscrEmpiricalDistFunc M5distFuncSS = seqSpecificETASmodel.computeNum_DistributionFunc(dataMaxDays+forecastMinDays, 
				dataMaxDays+31, refMag);
		
		System.out.println("dataEndTime = " + seqSpecificETASmodel.dataEndTimeDays);
		
		
		EvenlyDiscretizedFunc cumNumWithTime = seqSpecificETASmodel.getModalCumNumEventsWithTime( 5.0, dataMaxDays+forecastMinDays,  dataMaxDays+forecastMaxDays,  1);
		System.out.println(cumNumWithTime);
		 
//		
//		/*
//		 * Calculate Bayesian combination model
//		 */
//		RJ_AftershockModel_Bayesian bayesianModel = new RJ_AftershockModel_Bayesian(seqSpecificModel, genericModel);
//		

		
//		/*
//		 * Print results
//		 */
//		double minMag = 3d;
//		double maxMag = 9d;
//		double deltaMag = 0.1;
//		int numMag = (int)((maxMag - minMag)/deltaMag + 1.5);
//		
//		EvenlyDiscretizedFunc mode = genericModel.getModalCumNumMFD(minMag, maxMag, numMag, forecastMinDays, forecastMaxDays);
//		mode.setName("Generic Forecast Mode");
//		System.out.println(mode);
		
//		mode = seqSpecificModel.getModalCumNumMFD(minMag, maxMag, numMag, forecastMinDays, forecastMaxDays);
//		mode.setName("Sequence Specific Forecast Mode");
//		System.out.println(mode);
//		
//		mode = bayesianModel.getModalCumNumMFD(minMag, maxMag, numMag, forecastMinDays, forecastMaxDays);
//		mode.setName("Bayesian Forecast Mode");
//		System.out.println(mode);
	}

}
