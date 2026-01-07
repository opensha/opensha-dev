package scratch.kevin.prvi25;

import java.io.IOException;
import java.util.List;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTree;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionCouplingModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionFaultModels;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;

public class InterfaceAverageSlipRateCalc {
	
	private enum OriginalInterfaceSection {
		HISPANIOLA("Hispaniola") {
			@Override
			public boolean isMember(FaultSection subSect) {
				return subSect.getParentSectionId() == 7501;
			}
		},
		PRVI("PRVI") {
			@Override
			public boolean isMember(FaultSection subSect) {
				return subSect.getParentSectionId() == 7500 && traceCenterLon(subSect) < -63.06423033396625d;
			}
		},
		LESER_ANTILLES("Lesser Antilles") {
			@Override
			public boolean isMember(FaultSection subSect) {
				return subSect.getParentSectionId() == 7500 && traceCenterLon(subSect) >= -63.06423033396625d;
			}
		},
		MUERTOS("Muertos") {
			@Override
			public boolean isMember(FaultSection subSect) {
				return subSect.getParentSectionId() == 7550;
			}
		};
		
		private String name;

		private OriginalInterfaceSection(String name) {
			this.name = name;
		}
		
		@Override
		public String toString() {
			return name;
		}
		
		public abstract boolean isMember(FaultSection subSect);
		
		private static double traceCenterLon(FaultSection sect) {
			FaultTrace trace = sect.getFaultTrace();
			return 0.5*(trace.first().lon + trace.last().lon);
		}
		
		public static OriginalInterfaceSection detect(FaultSection subSect) {
			for (OriginalInterfaceSection sect : values())
				if (sect.isMember(subSect))
					return sect;
			throw new IllegalStateException("Couldn't detect original section for "+subSect.getSectionName());
		}
	}

	public static void main(String[] args) throws IOException {
		LogicTreeBranch<LogicTreeNode> branch = PRVI25_LogicTree.DEFAULT_SUBDUCTION_INTERFACE;
		OriginalInterfaceSection[] origSects = OriginalInterfaceSection.values();
		AreaWeightedSlipRateAverage[] averages = new AreaWeightedSlipRateAverage[origSects.length];
		for (int i=0; i<averages.length; i++)
			averages[i] = new AreaWeightedSlipRateAverage();
		for (PRVI25_SubductionFaultModels fm : PRVI25_SubductionFaultModels.values()) {
			branch.setValue(fm);
			List<? extends FaultSection> subSects = fm.buildSubSects(fm);
			for (PRVI25_SubductionDeformationModels dm : PRVI25_SubductionDeformationModels.values()) {
//				if (dm == PRVI25_SubductionDeformationModels.PARTIAL)
//					continue;
				branch.setValue(dm);
				for (PRVI25_SubductionCouplingModels coupling : PRVI25_SubductionCouplingModels.values()) {
//					if (coupling != PRVI25_SubductionCouplingModels.PREFERRED)
//						continue;
					branch.setValue(coupling);
					List<? extends FaultSection> dmSects = dm.apply(fm, branch, subSects);
					for (FaultSection subSect : dmSects) {
						OriginalInterfaceSection sect = OriginalInterfaceSection.detect(subSect);
						averages[sect.ordinal()].add(subSect.getArea(false), subSect.getOrigAveSlipRate());
					}
				}
			}
		}
		
		for (int i=0; i<averages.length; i++) {
			System.out.println(origSects[i].name+":\t"+(float)averages[i].getAverageSlipRate());
		}
	}
	
	private static class AreaWeightedSlipRateAverage {
		private double sumArea;
		private double sumAreaTimesSlip;
		
		public void add(double area, double slipRate) {
			sumArea += area;
			sumAreaTimesSlip += slipRate*area;
		}
		
		public double getAverageSlipRate() {
			return sumAreaTimesSlip/sumArea;
		}
	}

}