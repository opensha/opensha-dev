package scratch.kevin.risk;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.TimeSpan;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.ParameterList;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.PointSource;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.util.FocalMech;
import org.opensha.sra.gui.portfolioeal.Asset;
import org.opensha.sra.gui.portfolioeal.PortfolioParser;
import org.opensha.sra.vulnerability.VulnerabilityFetcher;

import com.google.common.base.Preconditions;

import scratch.UCERF3.erf.mean.MeanUCERF3;
import scratch.UCERF3.erf.mean.MeanUCERF3.Presets;

public class U3_EAL_ExampleCalc {

	public static void main(String[] args) throws IOException, ClassNotFoundException, InstantiationException, IllegalAccessException {
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
		// here's how I call it to preload the vulnerabilities and skip the servlet; uncomment and change to your path
		VulnerabilityFetcher.getVulnerabilities(new File("/home/kevin/OpenSHA/portfolio_lec/2011_11_07_VUL06.txt"));
		// or you can put in your own custom implementation like this. we need to make this better...
//		Vulnerability vuln = new SimpleVulnerability(name, shortName, imType, imLevels, mfdVals, covVals);
//		VulnerabilityFetcher.getVulnerabilities().put("MyVulnName", vuln);
		
		// if true, will also write conditional probability of exceedance for each IML for each rupture
		boolean writeIndividualRupExceedances = true;
		
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
//		MeanUCERF3 erf = new MeanUCERF3();
//		erf.setPreset(Presets.FM3_1_BRANCH_AVG);
//		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
//		// set to poisson, 1 year
//		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
//		erf.getTimeSpan().setDuration(1d);
		// or use a simple test ERF with a GR point source at a single location
		// M5-8 range, shows some discrepancy at present
		double minMag = 5d;
		double deltaMag = 0.1d;
		int numMag = 30; // up to M8
		// single M7 test, matches nearly perfectly
//		double minMag = 7d;
//		double deltaMag = 0.1d; // not used
//		int numMag = 1;
		
		GutenbergRichterMagFreqDist grMFD = new GutenbergRichterMagFreqDist(minMag, numMag, deltaMag);
		Location srcLoc = new Location(37, -122); // make sure this is nearish to to your asset
		double totRate = 1d; // total cumulative event rate
		double b = 1d;
		double duration = 1d;
		grMFD.setAllButTotMoRate(grMFD.getMinX(), grMFD.getMaxX(), totRate, b);
		PointSource ptSrc = PointSource.poissonBuilder(srcLoc).truePointSources()
				.forMFDAndFocalMech(grMFD, FocalMech.STRIKE_SLIP.mechanism).duration(duration).build();
		System.out.println("Simple ERF ruptures:");
		for (ProbEqkRupture rup : ptSrc)
			System.out.println("\tM"+(float)rup.getMag()+" with P="+(float)rup.getProbability()+" and rate="+(float)rup.getMeanAnnualRate(duration));
		AbstractERF erf = new AbstractERF() {
			
			@Override
			public String getName() {
				return null;
			}
			
			@Override
			public void updateForecast() {
				if (this.timeSpan == null)
					this.timeSpan = new TimeSpan(TimeSpan.NONE, TimeSpan.YEARS);
			}
			
			@Override
			public ProbEqkSource getSource(int idx) {
				Preconditions.checkState(idx == 0, "only 1 source");
				return ptSrc;
			}
			
			@Override
			public int getNumSources() {
				return 1;
			}
		};
		erf.updateForecast();
		
		// max source-site dist in KM
		double maxDistance = 200d;
		
		// also need a GMM, let's use ASK 2014
		ScalarIMR gmm = AttenRelRef.ASK_2014.get();
		// or could use the full NSHM23-WUS model:
//		ScalarIMR gmm = AttenRelRef.USGS_NSHM23_ACTIVE.get();
		
		gmm.setParamDefaults();
		
		// this exists in order to supply all required GMM paremeters, e.g. basin depth, that aren't included in the asset
		// well set those as default. the location and Vs30 value will be overridden automatically when we do the calculation
		Site site = new Site();
		for (Parameter<?> param : gmm.getSiteParams())
			site.addParameter((Parameter<?>) param.clone());
		
		// conditional losses for each rupture
		double[][] rupConditionalLosses = asset.calculateExpectedLossPerRup(gmm, maxDistance, null, site, erf, null);
		double eal = asset.getAssetEAL();
		
		System.out.println("EAL from conditional losses: "+eal);
		
		// compute it the other way
		double calcEAL = asset.calculateEAL(gmm, maxDistance, site, erf, null);
		System.out.println("EAL from the curve: "+calcEAL);
		
		System.out.println("Difference: "+Math.abs(eal-calcEAL)+" ("+Math.abs(100d*(eal-calcEAL)/calcEAL)+" %)");
		
		// lets write out the individual losses to a CSV File
		CSVFile<String> csv = new CSVFile<>(false);
		
		List<String> header = new ArrayList<>();
		header.add("Source ID");
		header.add("Rupture ID");
		header.add("Magnitude");
		header.add("Rate");
		header.add("Conditional Loss");
		ArbitrarilyDiscretizedFunc logXVals = null;
		if (writeIndividualRupExceedances) {
			logXVals = new ArbitrarilyDiscretizedFunc();
			double[] imls = asset.getVulnModel().getIMLValues();
			header.add(gmm.getShortName()+" Conditional Exceedance Probs");
			for (double iml : imls) {
				header.add((float)iml+"");
				logXVals.set(Math.log(iml), 0d);
			}
		}
		csv.addLine(header);
		for (int sourceID=0; sourceID<erf.getNumSources(); sourceID++) {
			if (rupConditionalLosses[sourceID] == null)
				// too far away
				continue;
			ProbEqkSource source = erf.getSource(sourceID);
			for (int rupID=0; rupID<source.getNumRuptures(); rupID++) {
				ProbEqkRupture rup = source.getRupture(rupID);
				double rate = rup.getMeanAnnualRate(erf.getTimeSpan().getDuration());
				List<String> line = new ArrayList<>(header.size());
				line.add(sourceID+"");
				line.add(rupID+"");
				line.add((float)rup.getMag()+"");
				line.add(rate+"");
				line.add(rupConditionalLosses[sourceID][rupID]+"");
				if (writeIndividualRupExceedances) {
					line.add(""); // intentionally blank
					gmm.setEqkRupture(rup);
					gmm.getExceedProbabilities(logXVals);
					for (int i=0; i<logXVals.size(); i++)
						line.add((float)logXVals.getY(i)+"");
				}
				csv.addLine(line);
			}
		}
		
		// add EAL at the bottom
		csv.addLine("");
		csv.addLine("EAL:", eal+"");
		
		csv.writeToFile(outputCSV);
	}

}
