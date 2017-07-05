package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.disaggregation.DisaggregationCalculator;
import org.opensha.sha.calc.disaggregation.DisaggregationPlotData;
import org.opensha.sha.calc.hazardMap.HazardDataSetLoader;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.gcim.ui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.SiteParams.DepthTo1pt0kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_TypeParam;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.griddedSeismicity.GridSourceFileReader;
import scratch.UCERF3.griddedSeismicity.GridSourceProvider;
import scratch.UCERF3.utils.FaultSystemIO;

public class MehrdadDemo {

	public static void main(String[] args) throws IOException, DocumentException {
		boolean ucerf2 = false;
		
		int startYear = 2017;
		double duration = 50;
		AbstractERF erf;
		if (ucerf2) {
			erf = new MeanUCERF2();
			erf.setParameter(UCERF2.PROB_MODEL_PARAM_NAME, MeanUCERF2.PROB_MODEL_WGCEP_PREF_BLEND);
		} else {
			// UCERF3
			FaultSystemSolution sol = FaultSystemIO.loadSol(new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
					+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
			FaultSystemRupSet rupSet = sol.getRupSet();
			double[] rates = sol.getRateForAllRups();
			
			
			// modify some rupture rates
			rates[0] = 0; // disable rupture 0
			// create new solution with my modified rates
			sol = new FaultSystemSolution(rupSet, rates);
			// write new solution to disk
			FaultSystemIO.writeSol(sol, new File("/tmp/test_sol.zip"));
			List<FaultSectionPrefData> subSects = rupSet.getFaultSectionDataList();
			FaultSectionPrefData sect0 = subSects.get(0);
			for (FaultSectionPrefData sect : subSects)
				System.out.println(sect.getName());
			int parentSectID = 301; // Mopjave S
			// list of all ruptures which include Mojave S
			List<Integer> mojaveRups = rupSet.getRupturesForParentSection(parentSectID);
			double totMojaveRate = 0d;
			for (int rupIndex : mojaveRups)
				totMojaveRate += sol.getRateForRup(rupIndex);
			System.out.println("Total Mojave supra-seismogenic rate: "+totMojaveRate);
			
			// find all subsections that rupture with the given parent sect ID
			Map<Integer, Double> multiFaultParticRates = Maps.newHashMap();
			for (int rupIndex : rupSet.getRupturesForParentSection(parentSectID)) {
				// for each rupture which Mojave participates
				
				double rate = sol.getRateForRup(rupIndex);
				for (int sectIndex : rupSet.getSectionsIndicesForRup(rupIndex)) {
					// for each subsection of this rupture
					
					Double particRate = multiFaultParticRates.get(sectIndex);
					if (particRate == null)
						particRate = rate; // first rupture that we have encountered for this section
					else
						particRate += rate; // add it
					multiFaultParticRates.put(sectIndex, particRate);
				}
			}
			for (int sectIndex : multiFaultParticRates.keySet()) {
				FaultSectionPrefData sect = rupSet.getFaultSectionData(sectIndex);
				System.out.println("Section "+sectIndex+", "+sect.getSectionName()
					+" participates with parent section "+parentSectID+" with rate "
						+multiFaultParticRates.get(sectIndex));
			}
			
			// create ERF from sol
			erf = new FaultSystemSolutionERF(sol);
			erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
			erf.setParameter(HistoricOpenIntervalParam.NAME, startYear-1875d);
			erf.getTimeSpan().setStartTime(startYear);
			erf.getTimeSpan().setDuration(duration);
		}
		
		erf.updateForecast();
		ScalarIMR imr = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);
		imr.setParamDefaults();
		imr.setIntensityMeasure(PGV_Param.NAME);
		// create site
		Site site = new Site(new Location(34, -118));
		// add site params required by NGA-W2
		Vs30_Param vs30 = new Vs30_Param();
		vs30.setValue(760);
		site.addParameter(vs30);
		Vs30_TypeParam vs30type = new Vs30_TypeParam();
		vs30type.setValue(Vs30_TypeParam.VS30_TYPE_INFERRED);
		site.addParameter(vs30type);
		DepthTo2pt5kmPerSecParam z25 = new DepthTo2pt5kmPerSecParam();
		z25.setValue(null); // disable, use Vs30
		site.addParameter(z25);
		DepthTo1pt0kmPerSecParam z10 = new DepthTo1pt0kmPerSecParam();
		z10.setValue(null); // disable, use Vs30
		site.addParameter(z10);
		// create curve calculator
		HazardCurveCalculator calc = new HazardCurveCalculator();
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(imr.getIntensityMeasure().getName());
		// xvals here are in linear space, we need them in log space for calculation
		DiscretizedFunc logHazFunc = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<xVals.size(); i++)
			logHazFunc.set(Math.log(xVals.getX(i)), 0);
		// this calculates and stores y vals in logHazFunc
		calc.getHazardCurve(logHazFunc, site, imr, erf);
		System.out.println(logHazFunc);
		
		double iml = HazardDataSetLoader.getCurveVal(logHazFunc, true, 0.02); // iml at 2% prob
		
		DisaggregationCalculator disaggCalc = new DisaggregationCalculator();
		disaggCalc.disaggregate(iml, site, imr, erf, calc.getAdjustableParams());
		DisaggregationPlotData disaggData = disaggCalc.getDisaggPlotData();
	}
	
	private static void customRupSet() throws IOException, DocumentException {
		FaultSystemSolution sol = FaultSystemIO.loadSol(new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		FaultSystemRupSet origRupSet = sol.getRupSet();
		List<FaultSectionPrefData> faultSectionData = origRupSet.getFaultSectionDataList();
		// create new list that contains contents of old one, old one probably unmodifiable
		faultSectionData = new ArrayList<FaultSectionPrefData>(faultSectionData);
		// add a new fault
		FaultSectionPrefData newFault = new FaultSectionPrefData();
		newFault.setSectionId(faultSectionData.size());
		//newFault.setAveDip();
//		newFault.set
		faultSectionData.add(newFault);
		// create new rupture that's only on that fault
		List<List<Integer>> sectionForRups = origRupSet.getSectionIndicesForAllRups();
		sectionForRups = Lists.newArrayList(sectionForRups);
		// sections for new rup on that fault
		List<Integer> sectionsForRup = Lists.newArrayList();
		sectionsForRup.add(newFault.getSectionId());
		sectionForRups.add(sectionsForRup);
		// need to add to mag/rake/area/rate arrays
		double[] origMags = origRupSet.getMagForAllRups();
		double[] newMags = new double[origMags.length+1];
		System.arraycopy(origMags, 0, newMags, 0, origMags.length); // copy mags over
		newMags[newMags.length-1] = 7d; // mag of new rupture
		// TODO repeat for rake, area, rate
		double[] newRakes = null;
		double[] newAreas = null;
		double[] newRates = null; // get from sol not rupSet
		
		// grid sources
		GridSourceFileReader gridProv = (GridSourceFileReader)sol.getGridSourceProvider();
		gridProv.writeGriddedSeisFile(new File("/tmp/grid_sources.xml"));
		// modify externally
		gridProv = GridSourceFileReader.fromFile(new File("/tmp/grid_sources.xml"));
		
		FaultSystemRupSet rupSet = new FaultSystemRupSet(
				faultSectionData, null, null, null, sectionForRups, newMags, newRakes, newAreas, null, "metadata");
		FaultSystemSolution newSol = new FaultSystemSolution(rupSet, newRates);
		newSol.setGridSourceProvider(sol.getGridSourceProvider());
		// write to file
		FaultSystemIO.writeSol(newSol, new File("/tmp/test_sol.zip"));
	}

}
