package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.modules.AverageableModule.AveragingAccumulator;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.PRVI25_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.PRVI25_GridSourceBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTree;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionCaribbeanSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionMuertosSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.collect.ImmutableList;

public class SubductionDefModAverageTest {

	public static void main(String[] args) throws IOException {
//		Location loc = new Location(18.90000, -65.90000);
//		Location loc = new Location(18.2, -62.6);
		Location loc = new Location(18.0, -68.2);
		
		List<LogicTreeLevel<? extends LogicTreeNode>> allLevels = new ArrayList<>();
		allLevels.addAll(PRVI25_LogicTree.levelsSubduction);
		allLevels.addAll(PRVI25_LogicTree.levelsSubductionGridded);
		LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(allLevels);
		for (LogicTreeNode node : PRVI25_LogicTree.DEFAULT_SUBDUCTION_INTERFACE)
			branch.setValue(node);
		branch.setValue(PRVI25_SubductionCaribbeanSeismicityRate.PREFFERRED);
		branch.setValue(PRVI25_SubductionMuertosSeismicityRate.PREFFERRED);
		branch.setValue(PRVI25_DeclusteringAlgorithms.AVERAGE);
		branch.setValue(PRVI25_SeisSmoothingAlgorithms.AVERAGE);
		
//		PRVI25_SeismicityRegions seisReg = PRVI25_SeismicityRegions.CAR_INTERFACE;
		PRVI25_SeismicityRegions seisReg = PRVI25_SeismicityRegions.MUE_INTERFACE;
		
		PRVI25_InvConfigFactory factory = new PRVI25_InvConfigFactory();
		
		List<PRVI25_SubductionScalingRelationships> scales = new ArrayList<>();
		List<GridSourceList> gridProvs = new ArrayList<>();
		for (PRVI25_SubductionScalingRelationships scale : PRVI25_SubductionScalingRelationships.values()) {
			if (scale.getNodeWeight(branch) == 0d)
				continue;
			branch.setValue(scale);
			
			FaultSystemRupSet rupSet = factory.buildRuptureSet(branch, 20);
			
			FaultSystemSolution sol = new FaultSystemSolution(rupSet, new double[rupSet.getNumRuptures()]);
			
			GridSourceList gridProv = PRVI25_GridSourceBuilder.buildInterfaceGridSourceList(
					sol, branch, seisReg);
			
			scales.add(scale);
			gridProvs.add(gridProv);
		}
		
		branch.setValue(PRVI25_SubductionScalingRelationships.AVERAGE);
		FaultSystemRupSet rupSet = factory.buildRuptureSet(branch, 20);
		FaultSystemSolution sol = new FaultSystemSolution(rupSet, new double[rupSet.getNumRuptures()]);
		
		System.out.println("Averaging");
		AveragingAccumulator<GridSourceProvider> accumulator = gridProvs.get(0).averagingAccumulator();
		double sumWeight = 0d;
		for (int i=0; i<scales.size(); i++) {
			double weight = scales.get(i).getNodeWeight(branch);
			sumWeight += weight;
			accumulator.process(gridProvs.get(i), weight);
		}
		
		GridSourceList avgGridProv = (GridSourceList) accumulator.getAverage();
		sol.setGridSourceProvider(avgGridProv);
		sol.write(new File("/tmp/sol_with_"+seisReg.name()+"_gridded.zip"));
		
		int locIndex = avgGridProv.getLocationIndex(loc);
		System.out.println("LocIndex="+locIndex);
		
		if (locIndex < 0)
			return;
		
		ImmutableList<GriddedRupture> avgRups = avgGridProv.getRuptures(TectonicRegionType.SUBDUCTION_INTERFACE, locIndex);
		List<ImmutableList<GriddedRupture>> origRups = new ArrayList<>();
		for (GridSourceList gridProv : gridProvs)
			origRups.add(gridProv.getRuptures(TectonicRegionType.SUBDUCTION_INTERFACE, locIndex));
		
		for (int r=0; r<avgRups.size(); r++) {
			GriddedRupture rup = avgRups.get(r);
			if (r > 0 && (float)rup.properties.magnitude == (float)avgRups.get(r-1).properties.magnitude)
				continue;
//			if ((float)rup.properties.magnitude != 7.85f)
//				continue;
			System.out.println("M"+(float)rup.properties.magnitude);
			List<Double> lengths = new ArrayList<>();
			List<String> strs = new ArrayList<>();
			double sumNonAvg = 0d;
			for (int i=-1; i<gridProvs.size(); i++) {
				String prefix = i < 0 ? "Average" : scales.get(i).getShortName();
				ImmutableList<GriddedRupture> rups = i < 0 ? avgRups : origRups.get(i);
				double sumRate = 0d;
				for (GriddedRupture oRup : rups) {
					if ((float)rup.properties.magnitude == (float)oRup.properties.magnitude) {
						String str = prefix+":\tlen="
								+(float)oRup.properties.length
								+"\tupper="+(float)oRup.properties.upperDepth
								+"\tlower="+(float)oRup.properties.lowerDepth
								+"\trawHypoDepth="+(float)oRup.properties.hypocentralDepth
								+"\thypoDepth="+(float)oRup.properties.getHypocentralDepth()
								+"\trate="+(float)oRup.rate;
						lengths.add(oRup.properties.length);
						strs.add(str);
						sumRate += oRup.rate;
					}
				}
				System.out.println("\t"+prefix+" sumRate:\t"+(float)sumRate);
				if (i >= 0)
					sumNonAvg += scales.get(i).getNodeWeight(branch)*sumRate;
			}
			sumNonAvg /= sumWeight;
			System.out.println("\tCalculated average:\t"+(float)sumNonAvg);
			System.out.println("\tRuptures:");
			for (String str : ComparablePairing.getSortedData(lengths, strs))
				System.out.println("\t\t"+str);
		}
	}

}
