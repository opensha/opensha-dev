package scratch.kevin.pointSources.paperFigs2026;

import static scratch.kevin.pointSources.paperFigs2026.ConstantsAndSettings.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.HashSet;
import java.util.Set;

import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.ClassUtils;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Joiner;

import gov.usgs.earthquake.nshmp.model.HazardModel;
import gov.usgs.earthquake.nshmp.model.NshmErf;
import scratch.kevin.latex.LaTeXUtils;

public class NSHM23GriddedStatsTexWriter {

	public static void main(String[] args) throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(ORIG_SOL_FILE);
		
		double[] magThresholds = {5d, 6d, 7d};
		String[] magNames = {"Five", "Six", "Seven"};
		
		double griddedMomentRate = 0;
		double[] griddedRates = new double[magThresholds.length];
		
		Region wusRegion = NSHM23_RegionLoader.loadFullConterminousWUS();
		
		GridSourceList gridList = sol.requireModule(GridSourceList.class);
		
		int gridNodesInRegion = 0;
		for (int l=0; l<gridList.getNumLocations(); l++) {
			if (!wusRegion.contains(gridList.getLocation(l)))
				continue;
			gridNodesInRegion++;
			for (TectonicRegionType trt : gridList.getTectonicRegionTypes()) {
				for (GriddedRupture rup : gridList.getRuptures(trt, l)) {
					double mo = MagUtils.magToMoment(rup.properties.magnitude);
					griddedMomentRate += mo*rup.rate;
					for (int m=0; m<magThresholds.length; m++)
						if (rup.properties.magnitude >= magThresholds[m])
							griddedRates[m] += rup.rate;
				}
			}
		}
		
		double faultMomentRate = sol.getTotalFaultSolutionMomentRate();
		double[] faultRates = new double[magThresholds.length];
		FaultSystemRupSet rupSet = sol.getRupSet();
		double[] faultRupsInRegion = rupSet.getFractRupsInsideRegion(wusRegion, false);
		for (int rupIndex=0; rupIndex<rupSet.getNumRuptures(); rupIndex++) {
			double mag = rupSet.getMagForRup(rupIndex);
			double rate = sol.getRateForRup(rupIndex) * faultRupsInRegion[rupIndex];
			for (int m=0; m<magThresholds.length; m++)
				if (mag >= magThresholds[m])
					faultRates[m] += rate;
		}
		
		double[] totalRates = new double[magThresholds.length];
		for (int m=0; m<magThresholds.length; m++)
			totalRates[m] = griddedRates[m] + faultRates[m];
		double totalMomentRate = faultMomentRate + griddedMomentRate;
		
		FileWriter texFW = new FileWriter(new File(FIGURES_DIR, "nshm23_fault_grid_stats.tex"));
		
		DecimalFormat pDF = new DecimalFormat("0%");
		texFW.write("% NSHM23-WUS stats\n");
		texFW.write(LaTeXUtils.defineValueCommand("WUSFaultCount",
				LaTeXUtils.groupedIntNumber(NSHM23_FaultModels.WUS_FM_v3.getFaultSections().size()))+"\n");
		texFW.write(LaTeXUtils.defineValueCommand("WUSGridCount",
				LaTeXUtils.groupedIntNumber(gridList.getNumLocations()))+"\n");
		texFW.write(LaTeXUtils.defineValueCommand("WUSGridInRegCount",
				LaTeXUtils.groupedIntNumber(gridNodesInRegion))+"\n");
		texFW.write(LaTeXUtils.defineValueCommand("WUSFaultMoPercent", pDF.format(faultMomentRate/totalMomentRate))+"\n");
		texFW.write(LaTeXUtils.defineValueCommand("WUSGridMoPercent", pDF.format(griddedMomentRate/totalMomentRate))+"\n");
		
		for (int m=0; m<magThresholds.length; m++) {
			texFW.write(LaTeXUtils.defineValueCommand("WUSFaultM"+magNames[m]+"Percent", pDF.format(faultRates[m]/totalRates[m]))+"\n");
			texFW.write(LaTeXUtils.defineValueCommand("WUSGridM"+magNames[m]+"Percent", pDF.format(griddedRates[m]/totalRates[m]))+"\n");
		}
		texFW.flush();
		
		HazardModel model = HazardModel.load(Path.of("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-conus-6.1.3"));
		Region ceus = NSHM23_RegionLoader.SeismicityRegions.CONUS_EAST.load();
		
		for (IncludeBackgroundOption bgOption : IncludeBackgroundOption.values()) {
			if (bgOption == IncludeBackgroundOption.INCLUDE)
				continue;
			NshmErf nshmERF = new NshmErf(model, Set.of(TectonicRegionType.STABLE_SHALLOW), bgOption);
			nshmERF.getTimeSpan().setDuration(1d);
			
			double moRate = 0d;
			double[] magRates = new double[magThresholds.length];
			nshmERF.updateForecast();
			
			System.out.println("Doing NSHM23-CEUS for "+bgOption);
			
			HashSet<String> sourceClasses = new HashSet<>();
			HashSet<String> surfClasses = new HashSet<>();
			
			int sourceCount = 0;
			int rupCount = 0;
			for (ProbEqkSource source : nshmERF) {
				double rateInside = 0d;
				double mMaxInside = 0d;
				for (ProbEqkRupture rup : source) {
					boolean inside = false;
					RuptureSurface surf = rup.getRuptureSurface();
					for (Location loc : surf.getEvenlyDiscritizedListOfLocsOnSurface()) {
						if (ceus.contains(loc)) {
							inside = true;
							break;
						}
					}
					if (inside) {
						rupCount++;
						surfClasses.add(surf.getClass().getName());
						double rate = rup.getMeanAnnualRate(1d);
						double mag = rup.getMag();
						moRate += MagUtils.magToMoment(mag)*rate;
						for (int m=0; m<magThresholds.length; m++)
							if (mag >= magThresholds[m])
								magRates[m] += rate;
						rateInside += rate;
						mMaxInside = Math.max(mMaxInside, mag);
					}
				}
				if (rateInside > 0d) {
//					System.out.println(source.getName()+": "+rateInside+" with Mmax="+mMaxInside);
					sourceCount++;
					sourceClasses.add(source.getClass().getName());
				}
			}
			
			System.out.println("Included "+rupCount+" ruptures from "+sourceCount+" sources");
			System.out.println("\tSource classes:"+Joiner.on("\n\t").join(sourceClasses));
			System.out.println("\tSurface classes:"+Joiner.on("\n\t").join(surfClasses));
			
			if (bgOption == IncludeBackgroundOption.EXCLUDE) {
				faultRates = magRates;
				faultMomentRate = moRate;
			} else if (bgOption == IncludeBackgroundOption.ONLY) {
				griddedRates = magRates;
				griddedMomentRate = moRate;
			}
		}
		
		for (int m=0; m<magThresholds.length; m++)
			totalRates[m] = griddedRates[m] + faultRates[m];
		totalMomentRate = faultMomentRate + griddedMomentRate;
		
		texFW.write("\n% NSHM23-CEUS stats\n");
		texFW.write(LaTeXUtils.defineValueCommand("CEUSFaultMoPercent", pDF.format(faultMomentRate/totalMomentRate))+"\n");
		texFW.write(LaTeXUtils.defineValueCommand("CEUSGridMoPercent", pDF.format(griddedMomentRate/totalMomentRate))+"\n");
		
		for (int m=0; m<magThresholds.length; m++) {
			texFW.write(LaTeXUtils.defineValueCommand("CEUSFaultM"+magNames[m]+"Percent", pDF.format(faultRates[m]/totalRates[m]))+"\n");
			texFW.write(LaTeXUtils.defineValueCommand("CEUSGridM"+magNames[m]+"Percent", pDF.format(griddedRates[m]/totalRates[m]))+"\n");
		}
		
		texFW.close();
	}

}
