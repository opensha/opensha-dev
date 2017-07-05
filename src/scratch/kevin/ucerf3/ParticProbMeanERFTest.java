package scratch.kevin.ucerf3;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSetMath;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.earthquake.param.ApplyGardnerKnopoffAftershockFilterParam;
import org.opensha.sha.earthquake.param.BPTAveragingTypeOptions;
import org.opensha.sha.earthquake.param.BPTAveragingTypeParam;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.analysis.CompoundFSSPlots.MapPlotData;
import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.analysis.CompoundFSSPlots.MapBasedPlot;
import scratch.UCERF3.analysis.CompoundFSSPlots.TimeDepGriddedParticipationProbPlot;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.mean.MeanUCERF3;

public class ParticProbMeanERFTest {

	public static void main(String[] args) throws GMT_MapException, IOException, DocumentException {
		MeanUCERF3 erf = new MeanUCERF3();
		erf.setMeanParams(0, false, 0, MeanUCERF3.RAKE_BASIS_NONE);
		erf.getParameter(ProbabilityModelParam.NAME).setValue(
				ProbabilityModelOptions.U3_PREF_BLEND);
		erf.setParameter(HistoricOpenIntervalParam.NAME,
				(double)(FaultSystemSolutionERF.START_TIME_DEFAULT-1875));
		erf.setParameter(BPTAveragingTypeParam.NAME,
				BPTAveragingTypeOptions.AVE_RI_AVE_NORM_TIME_SINCE);
		erf.setParameter(ApplyGardnerKnopoffAftershockFilterParam.NAME, false);
		
		double duration = 30d;
		erf.getTimeSpan().setDuration(duration);
		
		double minMag = 6.7;
		double maxMag = 10;
		
		erf.updateForecast();
		
		GriddedRegion griddedRegion = new CaliforniaRegions.RELM_TESTING_GRIDDED();

		File writeDir = new File("/tmp/timeDep");
		if (!writeDir.exists())
			writeDir.mkdir();
		
		// can do sub seis once as it doesn't vary
		GriddedGeoDataSet data = TimeDepGriddedParticipationProbPlot.getAsProbs(
				ERF_Calculator.getParticipationRatesInRegion(erf, griddedRegion, minMag, maxMag, null), duration);
		GriddedGeoDataSet.writeXYZFile(data, new File(writeDir, "mean_erf_ucerf3_data.txt"));
		
		data.log10();
		
		String plotPrefix = "meanUCERF3_data";
		CPT particCPT = FaultBasedMapGen.getParticipationCPT()
				.rescale(-5, 0);
		FaultBasedMapGen.plotMap(writeDir, plotPrefix, false,
					FaultBasedMapGen.buildMap(particCPT, null, null,
							data, griddedRegion.getSpacing(), griddedRegion,
							false, "MeanUCERF3"));
		System.out.println("DONE.");
		
		data.exp();
		
		List<MapPlotData> meanData = MapBasedPlot.loadPlotData(
				new File("/tmp/timeDep/gridded_participation_prob_plots.xml"));
		
		String compareName = "30_timedep_gridded_partic_prob_6.7+";
		MapPlotData comparePlot = null;
		for (MapPlotData d : meanData) {
			if (d.getFileName().equals(compareName)) {
				comparePlot = d;
				break;
			}
		}
		Preconditions.checkState(comparePlot  != null);
		
		CPT ratioCPT;
		try {
			ratioCPT = FaultSysSolutionERF_Calc.getScaledLinearRatioCPT(0.02);
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		GeoDataSet compGridded = comparePlot.getGriddedData();
		compGridded.exp();
		
		GriddedGeoDataSet.writeXYZFile(compGridded, new File(writeDir, "mean_of_branches_ucerf3_data.txt"));
		
		GeoDataSet ratioData = GeoDataSetMath.divide(data, compGridded);
		FaultBasedMapGen.plotMap(writeDir, "meanUCERF3_ratio", false,
				FaultBasedMapGen.buildMap(ratioCPT, null, null,
						ratioData, griddedRegion.getSpacing(), griddedRegion,
						false, "MeanUCER3/true mean"));
	}

}
