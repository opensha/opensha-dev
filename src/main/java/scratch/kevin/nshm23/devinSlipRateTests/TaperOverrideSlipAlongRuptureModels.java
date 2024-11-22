package scratch.kevin.nshm23.devinSlipRateTests;

import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.modules.SlipAlongRuptureModel;
import org.opensha.sha.earthquake.faultSysSolution.util.SlipAlongRuptureModelBranchNode;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

import scratch.kevin.nshm23.devinSlipRateTests.DevinModDeformationModels.SlipOverride;

public enum TaperOverrideSlipAlongRuptureModels implements SlipAlongRuptureModelBranchNode {
	UNIFORM("Uniform", "Uniform") {
		@Override
		public SlipAlongRuptureModel getModel() {
			return new SlipAlongRuptureModel.Uniform();
		}
	},
	TAPER_OVERRIDE_INDIVIDUAL("Taper override (select faults), individual-taper on matching parents", "TaperOverrideIndividual") {
		@Override
		public SlipAlongRuptureModel getModel() {
			return new InvidualFaultTaperOverrideSlipAlongRuptureModel(getOverrideParents(), true);
		}
	},
	TAPER_OVERRIDE_COMBINED("Taper override (select faults), single-taper across contiguous matching parents", "TaperOverrideCombined") {
		@Override
		public SlipAlongRuptureModel getModel() {
			return new InvidualFaultTaperOverrideSlipAlongRuptureModel(getOverrideParents(), false);
		}
	};
	
	private String name;
	private String shortName;
	
	private static Set<Integer> parentIDs;

	private TaperOverrideSlipAlongRuptureModels(String name, String shortName) {
		this.name = name;
		this.shortName = shortName;
	}
	
	private static Set<Integer> getOverrideParents() {
		if (parentIDs == null) {
			synchronized (TaperOverrideSlipAlongRuptureModels.class) {
				Set<Integer> parents;
				try {
					List<? extends FaultSection> subSects = DevinModDeformationModels.GEO_FROM_DEVIN.build(NSHM23_FaultModels.WUS_FM_v3);
					Map<Integer, SlipOverride> overrides = DevinModDeformationModels.GEO_FROM_DEVIN.getSlipRateOverrides(subSects);
					
					parents = new HashSet<>();
					for (SlipOverride override : overrides.values()) {
						Preconditions.checkState(subSects.get(override.id).getName().equals(override.name));
						int parentID = subSects.get(override.id).getParentSectionId();
						parents.add(parentID);
					}
				} catch (IOException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
				TaperOverrideSlipAlongRuptureModels.parentIDs = parents;
			}
		}
		return parentIDs;
	}

	@Override
	public double getNodeWeight(LogicTreeBranch<?> fullBranch) {
		return 1;
	}

	@Override
	public String getFilePrefix() {
		return name();
	}

	@Override
	public String getShortName() {
		return shortName;
	}

	@Override
	public String getName() {
		return name;
	}
	
	public abstract SlipAlongRuptureModel getModel();

}
