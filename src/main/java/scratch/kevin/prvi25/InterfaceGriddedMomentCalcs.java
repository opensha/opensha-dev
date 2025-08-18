package scratch.kevin.prvi25;

import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SubSeisMoRateReductions;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.PRVI25_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTree;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionBValues;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionScalingRelationships;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import scratch.kevin.prvi25.figures.PRVI_SubductionSubSectPlots;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

public class InterfaceGriddedMomentCalcs {

	public static void main(String[] args) throws IOException {
		FaultSystemSolution largeSol = FaultSystemSolution.load(SUBDUCTION_SOL_LARGE);
		FaultSystemSolution smallSol = FaultSystemSolution.load(SUBDUCTION_SOL_LARGE);
		
		PRVI25_SubductionScalingRelationships scale = PRVI25_SubductionScalingRelationships.AVERAGE;
		
		DecimalFormat pDF = new DecimalFormat("0.00%");
		
//		double maxGridMagToConsider = 7.5d;
		double maxGridMagToConsider = 10d;
		
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = List.of(
				PRVI25_LogicTree.SUB_SUPRA_B,
				NSHM23_LogicTreeBranch.SUB_SEIS_MO);
		LogicTreeBranch<LogicTreeNode> subRedBranch = new LogicTreeBranch<>(levels);
		subRedBranch.setValue(PRVI25_SubductionBValues.AVERAGE);
		subRedBranch.setValue(SubSeisMoRateReductions.SUB_B_1);
		
		File outputDir = new File("/tmp/dm_accounting");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		for (boolean large : new boolean[] { true, false }) {
			FaultSystemSolution sol = large ? largeSol : smallSol;
			FaultSystemRupSet rupSet = sol.getRupSet();
			
			GeographicMapMaker mapMaker = new GeographicMapMaker(PRVI_SubductionSubSectPlots.plotReg);
			mapMaker.setWriteGeoJSON(false);
			mapMaker.setFaultSections(rupSet.getFaultSectionDataList());
			mapMaker.setFillSurfaces(true);
			CPT slipCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, 4d);
			slipCPT.setNanColor(Color.GRAY);
			CPT ratioCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(0d, 2d);
			
			SectSlipRates origSlips = rupSet.getSectSlipRates();
			// try a sub-seis b=1 moment reduction
			new PRVI25_InvConfigFactory().buildInversionConfig(rupSet, subRedBranch, 32);
			SectSlipRates subRedSlips = rupSet.getSectSlipRates();
			Preconditions.checkState(subRedSlips != origSlips);
			
			GridSourceList gridSources = sol.requireModule(GridSourceList.class);
			
			for (boolean momentFromMag : new boolean[] { true, false }) {
				System.out.println((large ? "Large" : "Small")+" fault model, "+(momentFromMag ? "from magnitude" : "from slip and area"));
				
				double[] mappedSectGriddedMomentRates = new double[rupSet.getNumSections()];
				double sumGriddedMomentRate = 0d;
				for (int l=0; l<gridSources.getNumLocations(); l++) {
					for (GriddedRupture rup : gridSources.getRuptures(TectonicRegionType.SUBDUCTION_INTERFACE, l)) {
						if (rup.properties.magnitude > maxGridMagToConsider)
							continue;
						Preconditions.checkState(rup.associatedSections != null && rup.associatedSections.length == 1);
						int sectIndex = rup.associatedSections[0];
						double moment;
						if (momentFromMag) {
							moment = MagUtils.magToMoment(rup.properties.magnitude);
						} else {
							double ddw = rup.properties.getDownDipWidth();
							if (ddw == 0d || rup.properties.length == 0d)
								// too small, point source
								continue;
							double area = rup.properties.length * ddw * 1e6; // convert to m^2
							Preconditions.checkState(area > 0, "Bad area=%s for len=%s and ddw=%s",
									(float)area, (float)rup.properties.length, (float)ddw);
							double slip = scale.getAveSlip(area, rup.properties.length*1e3, ddw*1e3, ddw*1e3, Double.NaN);
							Preconditions.checkState(slip > 0, "Bad slip=%s", (float)slip);
							moment = FaultMomentCalc.getMoment(area, slip);
						}
						double momentRate = moment*rup.rate;
						sumGriddedMomentRate += momentRate;
						mappedSectGriddedMomentRates[sectIndex] += momentRate;
					}
				}
				
				double minFract = Double.POSITIVE_INFINITY;
				double maxFract = 0d;
				int count = 0;
				double sumFract = 0d;
				double sumSectMomentRate = 0d;
				int numExceed = 0;
				
				List<Double> gridEquivSlipRates = new ArrayList<>();
				List<Double> origSlipRates = new ArrayList<>();
				List<Double> reducedSlipRates = new ArrayList<>();
				List<Double> gridToFaultRatios = new ArrayList<>();
				
				for (int s=0; s<rupSet.getNumSections(); s++) {
					FaultSection sect = rupSet.getFaultSectionData(s);
					
					double origSlipRate = sect.getOrigAveSlipRate();
					double origMomentRate = sect.calcMomentRate(false);
					sumSectMomentRate += origMomentRate;
					double fract = mappedSectGriddedMomentRates[s] / origMomentRate;
					if (fract >= 1)
						numExceed++;
					double gridSlipRate = fract*origSlipRate;
					double reducedSlipRate = Math.max(0, origSlipRate - gridSlipRate);
					System.out.println("\t"+sect.getSectionName()+"\tred="+pDF.format(fract)+"\torigMo="+(float)origMomentRate
							+"\tgridMo="+(float)mappedSectGriddedMomentRates[s]
							+"\torigSlip="+(float)sect.getOrigAveSlipRate()+"\tredSlip="+(float)reducedSlipRate);
					
					minFract = Math.min(minFract, fract);
					maxFract = Math.max(maxFract, fract);
					
					sumFract += fract;
					count++;
					
					gridEquivSlipRates.add(gridSlipRate);
					origSlipRates.add(origSlipRate);
					reducedSlipRates.add(reducedSlipRate == 0d ? Double.NaN : reducedSlipRate);
					gridToFaultRatios.add(fract);
				}
				
				String prefix = "sub_dm_accounting";
				if (large)
					prefix += "_large";
				else
					prefix += "_small";
				if (momentFromMag)
					prefix += "_fromMag";
				else
					prefix += "_fromScale";
				
				mapMaker.plotSectScalars(origSlipRates, slipCPT, "DM Slip Rates (mm/yr)");
				mapMaker.plot(outputDir, prefix+"_dm_slips", " ");
				mapMaker.plotSectScalars(gridEquivSlipRates, slipCPT, "Gridded Equivalent Slip Rates (mm/yr)");
				mapMaker.plot(outputDir, prefix+"_grid_equiv_slips", " ");
				mapMaker.plotSectScalars(reducedSlipRates, slipCPT, "Reduced (by gridded) Slip Rates (mm/yr)");
				mapMaker.plot(outputDir, prefix+"_reduced_slips", " ");
				mapMaker.plotSectScalars(gridToFaultRatios, ratioCPT, "Grided Moment / DM Moment");
				mapMaker.plot(outputDir, prefix+"_moment_ratio", " ");
				
				System.out.println("\tRange: ["+pDF.format(minFract)+", "+pDF.format(maxFract)+"]");
				double aveFract = sumFract / count;
				System.out.println("\tAverage: "+pDF.format(aveFract));
				System.out.println("\t# exceedances: "+numExceed+" ("+pDF.format((float)numExceed/(float)rupSet.getNumSections())+")");
				System.out.println("\tOf total: "+pDF.format(sumGriddedMomentRate/sumSectMomentRate));
			}
			System.out.println((large ? "Large" : "Small")+" fault model, sub-seis b=1 moment reductions");
			double minFract = Double.POSITIVE_INFINITY;
			double maxFract = 0d;
			int count = 0;
			double sumFract = 0d;
			double sumSectMomentRate = 0d;
			double sumSectSubSeisMomentRate = 0d;
			for (int s=0; s<rupSet.getNumSections(); s++) {
				FaultSection sect = rupSet.getFaultSectionData(s);
				
				double origMomentRate = origSlips.calcMomentRate(s);
				double reducedMomentRate = subRedSlips.calcMomentRate(s);
				double subSeisMomentRate = origMomentRate - reducedMomentRate;
				sumSectMomentRate += origMomentRate;
				sumSectSubSeisMomentRate += subSeisMomentRate;
				double fract = subSeisMomentRate / origMomentRate;
				System.out.println("\t"+sect.getSectionName()+"\tred="+pDF.format(fract)+"\torigMo="+(float)origMomentRate
						+"\tsubSeisMo="+(float)subSeisMomentRate+"\torigSlip="+(float)(origSlips.getSlipRate(s)*1e3)
						+"\tredSlip="+(float)(subRedSlips.getSlipRate(s)*1e3));
				
				minFract = Math.min(minFract, fract);
				maxFract = Math.max(maxFract, fract);
				
				sumFract += fract;
				count++;
			}
			
			System.out.println("\tRange: ["+pDF.format(minFract)+", "+pDF.format(maxFract)+"]");
			double aveFract = sumFract / count;
			System.out.println("\tAverage: "+pDF.format(aveFract));
			System.out.println("\tOf total: "+pDF.format(sumSectSubSeisMomentRate/sumSectMomentRate));
		}
	}

}
