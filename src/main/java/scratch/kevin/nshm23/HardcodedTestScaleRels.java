package scratch.kevin.nshm23;

import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.logicTree.Affects;
import org.opensha.commons.logicTree.DoesNotAffect;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.RupSetScalingRelationship;

@DoesNotAffect(FaultSystemRupSet.SECTS_FILE_NAME)
@DoesNotAffect(FaultSystemRupSet.RUP_SECTS_FILE_NAME)
@Affects(FaultSystemRupSet.RUP_PROPS_FILE_NAME)
@Affects(FaultSystemSolution.RATES_FILE_NAME)
public enum HardcodedTestScaleRels implements RupSetScalingRelationship {
	ELLB_C_4p2 {
		@Override
		public double getMag(double area, double length, double width, double origWidth, double aveRake) {
			return  4.2 + Math.log10(area*1e-6);
		}
	},
	ELLB_C_4p1 {
		@Override
		public double getMag(double area, double length, double width, double origWidth, double aveRake) {
			return  4.1 + Math.log10(area*1e-6);
		}
	};

	@Override
	public double getNodeWeight(LogicTreeBranch<?> fullBranch) {
		return 1d;
	}

	@Override
	public String getFilePrefix() {
		return name();
	}

	@Override
	public String getShortName() {
		return name();
	}

	@Override
	public String getName() {
		return name();
	}

	@Override
	public double getAveSlip(double area, double length, double width, double origWidth, double aveRake) {
		double mag = getMag(area, length, width, origWidth, aveRake);
		double moment = MagUtils.magToMoment(mag);
		return FaultMomentCalc.getSlip(area, moment);
	}

	@Override
	public abstract double getMag(double area, double length, double width, double origWidth, double aveRake);

}
