package scratch.kevin.risk;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.TimeSpan;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.ParameterList;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sra.gui.portfolioeal.Asset;
import org.opensha.sra.gui.portfolioeal.PortfolioParser;
import org.opensha.sra.vulnerability.Vulnerability;
import org.opensha.sra.vulnerability.VulnerabilityFetcher;
import org.opensha.sra.vulnerability.models.SimpleVulnerability;

import com.google.common.base.Preconditions;

import scratch.UCERF3.erf.mean.MeanUCERF3;
import scratch.UCERF3.erf.mean.MeanUCERF3.Presets;

public class U3_EAL_ExampleCalc {

	public static void main(String[] args) throws IOException {
		// this is a pretty weird setup; you instantiate the asset by first giving a comma separated list of parameters
		String assetParamNames = "AssetGroup,AssetID,AssetName,Lat,Lon,ValHi,ValLo,Value,VulnModel,Vs30";
		// then will pass the values in to the asset; here's a sample from one of Keith's files
		// the vulnerability model is in the 2nd to last column, and will be looked up
		String assetValues = "CEA 2019 100pct proxy portfolio,347,06001400100,37.859419,-122.230389,8106,1054,1054,RM1L-m-RES1-DF,627.7137";
		Asset asset = new Asset(assetParamNames);
		
		// you'll notice that the 2nd to last column contained the vulnerability model name
		// that name is matched against a list of vulnerability functions; there is currently no way to just set your
		// own custom vulnerability. This is a pretty bad design that really only makes sense in a portfolio context
		// where you're parsing an input file.
		
		// Vulnerability models can be loaded in Keith's 'VUL06' file format. Here's how you could do that
//		VulnerabilityFetcher.getVulnerabilities(new File("/path/to/vul06.csv"));
		// or you can put in your own custom implementation like this. we need to make this better...
//		Vulnerability vuln = new SimpleVulnerability(name, shortName, imType, imLevels, mfdVals, covVals);
//		VulnerabilityFetcher.getVulnerabilities().put("MyVulnName", vuln);
		
		// output file
		File outputCSV = new File("/tmp/eal_demo.csv");
		
		List<String> assetArrayList = PortfolioParser.parseLine(assetValues);
		String[] assetValueList = assetArrayList.toArray(new String[0]);
		ParameterList assetParams = asset.getParameterList();
		Preconditions.checkState(assetParams.size() == assetValueList.length,
				"Have %s params but %s values", assetParams.size(), assetValueList.length);
		System.out.println("Asset parameters (before parsing)");
		for (int p=0; p<assetValueList.length; p++)
			System.out.println(assetParams.getParameter(p).getName()+":\t"+assetValueList[p]);
		asset.setAssetParameters(assetValueList);
		
		// now we need an ERF, lets use UCERF3 (FM3.1, excluding background seismicity)
		MeanUCERF3 erf = new MeanUCERF3();
		erf.setPreset(Presets.FM3_1_BRANCH_AVG);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		// set to poisson, 1 year
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		erf.getTimeSpan().setDuration(1d);
		// build the forecast; this will download a data file the first time
		erf.updateForecast();
		
		// max source-site dist in KM
		double maxDistance = 200d;
		
		// also need a GMM, let's use ASK 2014
		ScalarIMR gmm = AttenRelRef.ASK_2014.get();
		gmm.setParamDefaults();
		
		// this exists in order to supply all required GMM paremeters, e.g. basin depth, that aren't included in the asset
		// well set those as default. the location and Vs30 value will be overridden automatically when we do the calculation
		Site site = new Site();
		for (Parameter<?> param : gmm.getSiteParams())
			site.addParameter((Parameter<?>) param.clone());
		
		// conditional losses for each rupture
		double[][] rupConditionalLosses = asset.calculateExpectedLossPerRup(gmm, maxDistance, null, site, erf, null);
		double eal = asset.getAssetEAL();
		
		System.out.println("EAL: "+eal);
		
		// lets write out the individual losses to a CSV File
		CSVFile<String> csv = new CSVFile<>(false);
		
		csv.addLine("Source ID", "Rupture ID", "Magnitude", "Rate", "Conditional Loss");
		for (int sourceID=0; sourceID<erf.getNumSources(); sourceID++) {
			if (rupConditionalLosses[sourceID] == null)
				// too far away
				continue;
			ProbEqkSource source = erf.getSource(sourceID);
			for (int rupID=0; rupID<source.getNumRuptures(); rupID++) {
				ProbEqkRupture rup = source.getRupture(rupID);
				double rate = rup.getMeanAnnualRate(erf.getTimeSpan().getDuration());
				csv.addLine(sourceID+"", rupID+"", (float)rup.getMag()+"", rate+"", rupConditionalLosses[sourceID][rupID]+"");
			}
		}
		
		// add EAL at the bottom
		csv.addLine("");
		csv.addLine("EAL:", eal+"");
		
		csv.writeToFile(outputCSV);
	}

}
