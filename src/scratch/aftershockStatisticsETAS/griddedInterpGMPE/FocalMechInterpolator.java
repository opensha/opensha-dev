package scratch.aftershockStatisticsETAS.griddedInterpGMPE;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.opensha.commons.data.Site;
import org.opensha.nshmp2.util.FocalMech;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.imr.ScalarIMR;

import com.google.common.collect.Lists;

import scratch.aftershockStatisticsETAS.griddedInterpGMPE.AbstractGMPEInterpolation.Discrete;

//public class FocalMechInterpolator extends Discrete<FocalMech> {
//	
//	public FocalMechInterpolator() {
//		this(FocalMech.values());
//	}
//
//	public FocalMechInterpolator(FocalMech... values) {
//		super("Focal Mechanism", Lists.newArrayList(values));
//	}
//
//	@Override
//	public void setGMPE_Params(ScalarIMR gmpe, ProbEqkSource source, int index) {}
//
//	@Override
//	public Set<EqkRupture> getViableRuptures(Set<EqkRupture> ruptures, int index) {
//		FocalMech mech = getValue(index);
//		HashSet<EqkRupture> viableRups = new HashSet<>();
//		
//		for (EqkRupture rup : ruptures)
//			if ((float)rup.getAveRake() == (float)mech.rake())
//				viableRups.add(rup);
//		
//		return viableRups;
//	}
//
//	@Override
//	public FocalMech detectCurrentVal(ScalarIMR gmpe, Site site) {
//		return forRake(gmpe.getEqkRupture().getAveRake());
//	}
//	
//	public static FocalMech forRake(double rake) {
//		for (FocalMech mech : FocalMech.values())
//			if ((float)rake == (float)mech.rake())
//				return mech;
//		throw new IllegalStateException("Unknown rake detected: "+rake);
//	}
//
//}
