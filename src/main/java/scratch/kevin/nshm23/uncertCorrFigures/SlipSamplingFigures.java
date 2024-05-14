package scratch.kevin.nshm23.uncertCorrFigures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets.CoulombRupSetConfig;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SlipRatePlots;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.util.SubSectionBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;

import com.google.common.base.Preconditions;

import scratch.kevin.nshm23.dmCovarianceTests.BvalAndSegConnectivityCorrelationSampler;
import scratch.kevin.nshm23.dmCovarianceTests.SectionCovarianceSampler;
import scratch.kevin.nshm23.dmCovarianceTests.SlipRateCovarianceSampler;

public class SlipSamplingFigures {

	public static void main(String[] args) throws IOException {
		NSHM23_FaultModels fm = NSHM23_FaultModels.WUS_FM_v2;
		NSHM23_DeformationModels dm = NSHM23_DeformationModels.GEOLOGIC;
		NSHM23_DeformationModels.HARDCODED_FRACTIONAL_STD_DEV_UPPER_BOUND = 0d;
		
		int interpSkip = 0;
		double corrDist = 200d;
		double zeroDistCoeff = 0.95;
		double negativeCorrMaxDist = 30d;
		
		File outputDir = new File("/home/kevin/Documents/papers/2024_nshm23_uncert_correlation/figures/slip_sampling");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		FaultSystemRupSet rupSet = buildExampleRupSet();
		List<? extends FaultSection> subSects = rupSet.getFaultSectionDataList();
		SectionDistanceAzimuthCalculator distAzCalc = new SectionDistanceAzimuthCalculator(subSects);
		File cacheDir = null;
		int firstSect1ID = -1;
		int lastSect1ID = -1;
		for (FaultSection sect : subSects) {
			if (sect.getParentSectionId() == 1) {
				if (firstSect1ID < 0)
					firstSect1ID = sect.getSectionId();
				lastSect1ID = sect.getSectionId();
			}
		}
		int corrSubsectID = firstSect1ID + (int)(0.5*(1+lastSect1ID-firstSect1ID));
		Region plotReg = new Region(new Location(0.35, -0.2), new Location(1d, 0.2));
		double maxSlipRate = 1.5d;
		boolean axisVisible = false;

//		File cacheDir = new File("/home/kevin/markdown/nshm23-misc/dm_slip_sampling/"
//				+ "full_GEOLOGIC_conn_corr_b0.5_MidSeg_fromMFDs_dist200.0km_zeroCoeff0.95_negCorrDist30.0km/cache/");
//		
//		List<? extends FaultSection> subSects = dm.build(fm);
//		SectionDistanceAzimuthCalculator distAzCalc = new SectionDistanceAzimuthCalculator(subSects);
//		
//		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
//		factory.setCacheDir(new File("/home/kevin/OpenSHA/nshm23/rup_sets/cache"));
//		FaultSystemRupSet rupSet = factory.buildRuptureSet(NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT, FaultSysTools.defaultNumThreads());
//		Preconditions.checkState(subSects.size() == rupSet.getNumSections());
//		
////		int corrSubsectID = 110;
////		Region plotReg = new Region(new Location(38, -120.25), new Location(39.5, -119));
////		double maxSlipRate = 0.5;
////		boolean axisVisible = true;
//		
////		int corrSubsectID = 2623;
////		Region plotReg = new Region(new Location(35.8, -114), new Location(37.3, -112.75));
////		double maxSlipRate = 0.5;
////		boolean axisVisible = true;
//		
////		int corrSubsectID = 4509;
////		Region plotReg = new Region(new Location(47.2, -123.1), new Location(47.8, -121.5));
////		double maxSlipRate = 0.5;
////		boolean axisVisible = true;
//		
//		int corrSubsectID = 2022;
//		Region plotReg = new Region(new Location(34, -107.2), new Location(35.5, -106));
//		double maxSlipRate = 0.3;
//		boolean axisVisible = true;
		
		SectionCovarianceSampler sampler = new BvalAndSegConnectivityCorrelationSampler(
				subSects, rupSet, distAzCalc, corrDist, zeroDistCoeff, negativeCorrMaxDist, 0.5d, NSHM23_SegmentationModels.MID, true);
		
		RandomGenerator rng = new Well19937c(1234l);
		boolean center = false;
		
		sampler.setDebug(true);
		
		Region calcReg = new Region(new Location(plotReg.getMinLat()-1, plotReg.getMinLon()-1),
				new Location(plotReg.getMaxLat()+1, plotReg.getMaxLon()+1));
		GeographicMapMaker mapMaker = new GeographicMapMaker(plotReg);
		mapMaker.setAxisLabelsVisible(axisVisible);
		mapMaker.setAxisTicksVisible(axisVisible);
		mapMaker.setFaultSections(subSects);
		mapMaker.setScalarThickness(5f);
		mapMaker.setDefaultPlotWidth(650);
		
		FaultSection corrSect = subSects.get(corrSubsectID);
		double[] cors = new double[subSects.size()];
		for (int i=0; i<subSects.size(); i++) {
			if (i == corrSubsectID) {
				cors[i] = 1d;
			}else {
				FaultSection oSect = subSects.get(i);
				Location traceStart = oSect.getFaultTrace().first();
				Location traceEnd = oSect.getFaultTrace().last();
				Location middle = new Location(0.5d*(traceStart.getLatitude()+traceEnd.getLatitude()),
						0.5*(traceStart.getLongitude()+traceEnd.getLongitude()));
				if (calcReg.contains(middle))
					cors[i] = sampler.getCorrelationCoefficient(corrSect, subSects.get(i));
			}
		}
		CPT corrCPT = GMT_CPT_Files.DIVERGING_BROC_UNIFORM.instance().rescale(-1d, 1d);
		mapMaker.setSectHighlights(List.of(corrSect), new PlotCurveCharacterstics(PlotLineType.SOLID, 13f, Color.BLACK));
		String corrLabel;
		if (axisVisible)
			// real data
			corrLabel = corrSect.getName()+" Correleation Coeff.";
		else
			corrLabel = "Example Fault Correlation Coeff.";
		mapMaker.plotSectScalars(cors, corrCPT, corrLabel);
		mapMaker.plot(outputDir, "corr_example", " ");
		mapMaker.clearSectHighlights();
		
//		CPT linearSlipCPT = SlipRatePlots.linearSlipCPT(1d);
		CPT linearSlipCPT = GMT_CPT_Files.SEQUENTIAL_NAVIA_UNIFORM.instance().rescale(0d, 1d).trim(0d, 0.9d).rescale(0d, maxSlipRate);
		CPT covCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, 1d);
		double[] slipRates = new double[subSects.size()];
		double[] slipRateCOVs = new double[subSects.size()];
		for (int i=0; i<subSects.size(); i++) {
			slipRates[i] = subSects.get(i).getOrigAveSlipRate();
			slipRateCOVs[i] = subSects.get(i).getOrigSlipRateStdDev()/slipRates[i];
		}
		
		mapMaker.plotSectScalars(slipRates, linearSlipCPT, "Preferred Slip Rate (mm/yr)");
		mapMaker.plot(outputDir, "slip_rates", " ");
		mapMaker.plotSectScalars(slipRateCOVs, covCPT, "Slip Rate COV");
		mapMaker.plot(outputDir, "slip_rate_covs", " ");
		
		int numSampleSets = 10;
		int numSamplesEach = 3;
		
		int numSamples = numSampleSets*numSamplesEach;
		
		if (cacheDir != null)
			sampler = sampler.loadCached(cacheDir, interpSkip).get();
		
		SlipRateCovarianceSampler slipSampler = new SlipRateCovarianceSampler(sampler);
		List<List<FaultSection>> samples = slipSampler.buildSamples(numSamples, rng, center, interpSkip);
		
		int sample = 0;
		for (int set=0; set<numSampleSets; set++) {
			for (int setN=0; setN<numSamplesEach; setN++) {
				slipRates = new double[subSects.size()];
				for (int i=0; i<subSects.size(); i++)
					slipRates[i] = samples.get(sample).get(i).getOrigAveSlipRate();
				mapMaker.plotSectScalars(slipRates, linearSlipCPT, "Slip Rate Sample "+(setN+1)+" (mm/yr)");
				mapMaker.plot(outputDir, "slip_rates_sample_"+set+"_"+setN, " ");
				sample++;
			}
		}
	}
	
	private static FaultSystemRupSet buildExampleRupSet() throws IOException {
		/*
		 * 		|		3
		 * 		|
		 * 	   / \
		 *    |   |		1  2
		 *     \ /
		 * 		|
		 * 		|		0
		 */
		double len0 = 50d;
		double jumpDist = 2d;
		double len12 = 50d;
		double width12 = 5d;
		double len3 = len0;
		List<Location[]> locs = new ArrayList<>();
		Location l01 = new Location(0d, 0d);
		Location l02 = LocationUtils.location(l01, 0d, len0);
		locs.add(new Location[] {
				l01,
				l02
		});
		
		Location l11 = LocationUtils.location(l02, 3d*Math.PI/2d, jumpDist);
		Location l21 = LocationUtils.location(l02, Math.PI/2d, jumpDist);
		
		double len12each = len12/3d;
		Location l12 = LocationUtils.location(l11, 0d, len12each); // move up
		l12 = LocationUtils.location(l12, 3d*Math.PI/2d, width12); // move left
		Location l22 = LocationUtils.location(l21, 0d, len12each); // move up
		l22 = LocationUtils.location(l22, Math.PI/2d, width12); // move right
		
		Location l13 = LocationUtils.location(l12, 0d, len12each); // move up;
		Location l23 = LocationUtils.location(l22, 0d, len12each); // move up;
		
		Location l14 = LocationUtils.location(l13, 0d, len12each); // move up
		l14 = LocationUtils.location(l14, Math.PI/2d, width12); // move right
		Location l24 = LocationUtils.location(l23, 0d, len12each); // move up
		l24 = LocationUtils.location(l24, 3*Math.PI/2d, width12); // move left
		locs.add(new Location[] {
				l11,
				l12,
				l13,
				l14
		});
		
		locs.add(new Location[] {
				l21,
				l22,
				l23,
				l24
		});
		
		Location l31 = LocationUtils.location(l02, 0d, len12);
		Location l32 = LocationUtils.location(l31, 0d, len3);
		locs.add(new Location[] {
				l31,
				l32
		});
		
		double rake = 180d;
		double dip = 90d;
		double lowDepth = 10d;
		double upDepth = 0d;
		double slipRate0 = 1d;
		double slipRate12 = 0.5d;
		double slipRate3 = 1d;
		double slipRateSD0 = 0.5d;
		double slipRateSD12 = 1d;
		double slipRateSD3 = 0.5d;
		
		List<GeoJSONFaultSection> sects = new ArrayList<>();
		
		for (int i=0; i<locs.size(); i++) {
			Location[] lls = locs.get(i);
			
			double slipRate;
			double slipRateSD;
			switch (i) {
			case 0:
				slipRate = slipRate0;
				slipRateSD = slipRateSD0;
				break;
			case 1:
				slipRate = slipRate12;
				slipRateSD = slipRateSD12;
				break;
			case 2:
				slipRate = slipRate12;
				slipRateSD = slipRateSD12;
				break;
			case 3:
				slipRate = slipRate3;
				slipRateSD = slipRateSD3;
				break;

			default:
				throw new IllegalStateException();
			}
			
			String sectJSON ="    {\n"+
							"      \"type\": \"Feature\",\n"+
							"      \"id\": "+i+",\n"+
							"      \"properties\": {\n"+
							"        \"FaultID\": "+i+",\n"+
							"        \"FaultName\": \"Example Fault "+(i+1)+"\",\n"+
							"        \"DipDeg\": "+(float)dip+",\n"+
							"        \"Rake\": "+(float)rake+",\n"+
							"        \"LowDepth\": "+(float)lowDepth+",\n"+
							"        \"UpDepth\": "+(float)upDepth+",\n"+
							"        \"SlipRate\": "+(float)slipRate+",\n"+
							"        \"SlipRateStdDev\": "+(float)slipRateSD+"\n"+
							"      },\n"+
							"      \"geometry\": {\n"+
							"        \"type\": \"LineString\",\n"+
							"        \"coordinates\": [\n";
			for (int l=0; l<lls.length; l++) {
				Location loc = lls[l];
				sectJSON +=	"          [\n"+
							"            "+loc.getLongitude()+",\n"+
							"            "+loc.getLatitude()+"\n"+
							"          ]";
				if (l < lls.length-1)
					sectJSON += ",";
				sectJSON += "\n";
			}
			sectJSON +=		"        ]\n"+
							"      }\n"+
							"    }";
			System.out.println(sectJSON);
			Feature feature = Feature.fromJSON(sectJSON);
			GeoJSONFaultSection sect = GeoJSONFaultSection.fromFeature(feature);
			sects.add(sect);
		}
		
		List<FaultSection> subSects = SubSectionBuilder.buildSubSects(sects);
		
		CoulombRupSetConfig rsConfig = new CoulombRupSetConfig(subSects, null, NSHM23_ScalingRelationships.LOGA_C4p2);
		rsConfig.setMaxJumpDist(6d);
		
		return rsConfig.build(FaultSysTools.defaultNumThreads());
	}

}
